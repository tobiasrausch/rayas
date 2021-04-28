#ifndef CALL_H
#define CALL_H

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

namespace ted
{

  struct CallConfig {
    uint16_t minMapQual;
    uint16_t minClip;
    uint16_t minSplit;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
    boost::filesystem::path tumor;
    boost::filesystem::path control;
  };



  template<typename TConfig, typename TVector>
  inline bool
  parseChr(TConfig& c, samFile* samfile, hts_idx_t* idx, bam_hdr_t* hdr, int32_t refIndex, std::string const& str, TVector& left, TVector& right, TVector& cov) {
    typedef typename TVector::value_type TValue;
    // Any data?
    bool nodata = true;
    std::string suffix("cram");
    if ((str.size() >= suffix.size()) && (str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0)) nodata = false;
    uint64_t mapped = 0;
    uint64_t unmapped = 0;
    hts_idx_get_stat(idx, refIndex, &mapped, &unmapped);
    if (mapped) nodata = false;
    if (nodata) return false;

    // Max value
    TValue maxval = std::numeric_limits<TValue>::max();

    // Read alignments
    hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
    bam1_t* rec = bam_init1();
    while (sam_itr_next(samfile, iter, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;
      std::size_t seed = hash_string(bam_get_qname(rec));
	
      // SV detection using single-end read
      uint32_t rp = rec->core.pos; // reference pointer
      uint32_t sp = 0; // sequence pointer
      
      // Parse the CIGAR
      uint32_t* cigar = bam_get_cigar(rec);
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k)
	    if (cov[rp] < maxval) ++cov[rp];
	  sp += bam_cigar_oplen(cigar[i]);
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
	  if (bam_cigar_oplen(cigar[i]) >= c.minClip) {
	    if (sp == 0) {
	      if (left[rp] < maxval) left[rp] += 1;
	    } else {
	      if (right[rp] < maxval) right[rp] += 1;
	    }
	    //if (rec->core.flag & BAM_FREAD1) read1.push_back(std::make_pair(seed, rp));
	    //else read2.push_back(std::make_pair(seed, rp));
	  }
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  rp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Warning: Unknown Cigar operation!" << std::endl;
	}
      }
    }
    bam_destroy1(rec);
    hts_itr_destroy(iter);
    return true;
  }


  template<typename TBitSet, typename TVector>
  inline void
  covParams(TBitSet const& nrun, TVector const& cov, uint32_t const seedwin, uint32_t& avgcov, uint32_t& sdcov) {
    boost::accumulators::accumulator_set<double, boost::accumulators::features<boost::accumulators::tag::mean, boost::accumulators::tag::variance> > acc;
    for(uint32_t i = seedwin; i < nrun.size(); i = i + seedwin) {
      uint32_t lcov = 0;
      bool ncontent = false;
      for(uint32_t k = i - seedwin; k < i; ++k) {
	if (nrun[k]) {
	  ncontent = true;
	  break;
	}
	lcov += cov[k];
      }
      if (!ncontent) acc(lcov);
    }
    sdcov = sqrt(boost::accumulators::variance(acc));
    avgcov = boost::accumulators::mean(acc);
  }

  template<typename TBitSet, typename TVector>
  inline bool
  getcov(TBitSet const& nrun, TVector const& cov, uint32_t const start, uint32_t const end, uint32_t& lcov) {
    lcov = 0;
    for(uint32_t k = start; k < end; ++k) {
      if (nrun[k]) return false;
      lcov += cov[k];
    }
    return true;
  }
  
  
  template<typename TConfig>
  inline int32_t
  runCall(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("ted.prof");
#endif

    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    samFile* cfile = sam_open(c.control.string().c_str(), "r");
    hts_set_fai_filename(cfile, c.genome.string().c_str());
    hts_idx_t* cidx = sam_index_load(cfile, c.control.string().c_str());
    faidx_t* fai = fai_load(c.genome.string().c_str());
    char* seq = NULL;
    
    // Parse genome, process chromosome by chromosome
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      if (hdr->target_len[refIndex] > 1000000) {
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();	  
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parsing " << hdr->target_name[refIndex] << std::endl;
      }
      
      // Tumor
      std::vector<uint16_t> left(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> right(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> cov(hdr->target_len[refIndex], 0);
      if (!parseChr(c, samfile, idx, hdr, refIndex, c.tumor.string(), left, right, cov)) continue;
      
      // Control
      std::vector<uint16_t> cleft(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> cright(hdr->target_len[refIndex], 0);
      std::vector<uint16_t> ccov(hdr->target_len[refIndex], 0);
      if (!parseChr(c, cfile, cidx, hdr, refIndex, c.control.string(), cleft, cright, ccov)) continue;

      // Load sequence
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, hdr->target_name[refIndex], 0, hdr->target_len[refIndex], &seqlen);
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr->target_len[refIndex], 0);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }
      if (seq != NULL) free(seq);
      
      // Get background coverage
      uint32_t seedwin = 200;
      if (2 * seedwin < hdr->target_len[refIndex]) {
	uint32_t sdcov = 0;
	uint32_t avgcov = 0;
	covParams(nrun, cov, seedwin, avgcov, sdcov);
	//std::cout << "Tumor avg. coverage and SD coverage " << avgcov << "," << sdcov << std::endl;
	uint32_t csdcov = 0;
	uint32_t cavgcov = 0;
	covParams(nrun, ccov, seedwin, cavgcov, csdcov);
	//std::cout << "Control avg. coverage and SD coverage " << cavgcov << "," << csdcov << std::endl;

	// Erase noise
	for(uint32_t i = seedwin; i < hdr->target_len[refIndex] - seedwin; ++i) {
	  if (left[i] >= c.minSplit) {
	    uint32_t lcov = 0;
	    if (!getcov(nrun, cov, i, i+seedwin, lcov)) continue;
	    if (lcov > avgcov + 3 * sdcov) std::cerr << hdr->target_name[refIndex] << ',' << i << ',' << left[i] << ",left," << lcov/seedwin << std::endl;
	  }
	  if (right[i] >= c.minSplit) {
	    uint32_t rcov = 0;
	    if (!getcov(nrun, cov, i - seedwin, i, rcov)) continue;
	    if (rcov > avgcov + 3 * sdcov) std::cerr << hdr->target_name[refIndex] << ',' << i << ',' << right[i] << ",right," << rcov/seedwin << std::endl;
	  }
	}
      }
    }
    
    // Clean-up
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    hts_idx_destroy(idx);
    sam_close(samfile);
    hts_idx_destroy(cidx);
    sam_close(cfile);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int call(int argc, char** argv) {
    CallConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("qual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(1), "min. mapping quality")
      ("clip,c", boost::program_options::value<uint16_t>(&c.minClip)->default_value(25), "min. clipping length")
      ("split,s", boost::program_options::value<uint16_t>(&c.minSplit)->default_value(3), "min. split-read support")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("matched,m", boost::program_options::value<boost::filesystem::path>(&c.control), "matched control BAM")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("sv.bcf"), "SV BCF output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.tumor), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("matched"))) {
      std::cout << "Usage: ted " << argv[0] << " [OPTIONS] -g <ref.fa> -m <control.bam> <tumor.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "ted ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runCall(c);
  }

}

#endif
