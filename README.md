[![C/C++ CI](https://github.com/tobiasrausch/rayas/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/rayas/actions)
[![Docker CI](https://github.com/tobiasrausch/rayas/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/rayas/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/rayas/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/rayas.svg)](https://github.com/tobiasrausch/rayas/releases)

# Rayas

Discovery of templated insertion threads and somatic retrocopy insertions

## Installing rayas

Rayas is available as a [statically linked binary](https://github.com/tobiasrausch/rayas/releases/), a [singularity container (SIF file)](https://github.com/tobiasrausch/rayas/releases/) or as a [docker container](https://hub.docker.com/r/trausch/rayas/).
You can also build Rayas from source using a recursive clone and make.
Rayas depends on [HTSlib](https://github.com/samtools/htslib) and [Boost](https://www.boost.org/).

`git clone --recursive https://github.com/tobiasrausch/rayas.git`

`cd rayas/`

`make all`


## Running rayas

Rayas was designed for paired-end or single-end matched tumor-normal data. It can be used with short reads, for long reads please use [lorax](https://github.com/tobiasrausch/lorax).

`rayas call -g <genome.fa> -m <control.bam> <tumor.bam>`

## Simple graph visualization

You can convert the output into a dot graph. Each component represents one templated insertion cluster. Nodes are genomic segments and edges represent the cancer genome structure with edge weights equalling the sequencing read support.

`cut -f 4,9 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`

## Somatic retrocopy insertions

For somatic retrocopy insertions, the distance between exons tends to be smaller than the default cutoff of 10kbp. To detect clusters involving retrocopies you have to lower the segment distance threshold:

`rayas call -e 1000 -g <genome.fa> -m <control.bam> <tumor.bam>`

To catch smaller exons, you may also want to adjust the minimal size:

`rayas call -d 5 -i 50 -e 1000 -g <genome.fa> -m <control.bam> <tumor.bam>`


## Citation

Tobias Rausch, Rene Snajder, Adrien Leger, Milena Simovic, Mădălina Giurgiu, Laura Villacorta, Anton G. Henssen, Stefan Fröhling, Oliver Stegle, Ewan Birney, Marc Jan Bonder, Aurelie Ernst, Jan O. Korbel     
Long-read sequencing of diagnosis and post-therapy medulloblastoma reveals complex rearrangement patterns and epigenetic signatures     
Cell Genomics, 2023, 100281, [DOI: 10.1016/j.xgen.2023.100281](https://doi.org/10.1016/j.xgen.2023.100281)     

License
-------
Rayas is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/rayas/blob/master/LICENSE) file for more details.
