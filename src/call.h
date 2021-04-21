#ifndef CALL_H
#define CALL_H

#include <fstream>
#include <iomanip>

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

namespace ted
{

  struct CallConfig {
    boost::filesystem::path gtfFile;
    boost::filesystem::path annofile;
    boost::filesystem::path db;
    boost::filesystem::path matchfile;
    boost::filesystem::path infile;
  };


  template<typename TConfig>
  inline int32_t
  runCall(TConfig& c) {
    
#ifdef PROFILE
    ProfilerStart("ted.prof");
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
      ("anno,a", boost::program_options::value<boost::filesystem::path>(&c.annofile)->default_value("anno.bcf"), "output annotation VCF/BCF file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "query VCF/BCF file")
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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: ted " << argv[0] << " [OPTIONS] input.bcf" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // SV database
    if (!vm.count("db")) {
      // Set input SV file as DB to fill chr array
      c.db = c.infile;
    }
    
    // Match SV types
    if (vm.count("notype")) c.matchSvType = false;
    else c.matchSvType = true;

    // Report no matches
    if (vm.count("nomatch")) c.reportNoMatch = true;
    else c.reportNoMatch = false;

    // Report contained genes
    if (vm.count("contained")) c.containedGenes = true;
    else c.containedGenes = false;
    
    // Check size ratio
    if (c.sizediff < 0) c.sizediff = 0;
    else if (c.sizediff > 1) c.sizediff = 1;

    // Check strategy
    if (strategy == "all") c.bestMatch = false;
    else c.bestMatch = true;

    // Check output directory
    if (!_outfileValid(c.matchfile)) return 1;
    if (!_outfileValid(c.annofile)) return 1;

    // GTF/GFF3/BED
    if (vm.count("gtf")) {
      if (is_gff3(c.gtfFile)) c.gtfFileFormat = 2; // GFF3
      else if (is_gtf(c.gtfFile)) c.gtfFileFormat = 0; // GTF/GFF2
      else c.gtfFileFormat = 1;  // BED
    } else c.gtfFileFormat = -1;

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
