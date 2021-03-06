[![C/C++ CI](https://github.com/tobiasrausch/rayas/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/rayas/actions)
[![Docker CI](https://github.com/tobiasrausch/rayas/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/rayas/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/rayas/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/rayas.svg)](https://github.com/tobiasrausch/rayas/releases)

# Rayas

Discovery of templated insertion threads

## Running rayas

Rayas was designed for paired-end or single-end matched tumor-normal data. It can be used with short reads, for long reads please use [lorax](https://github.com/tobiasrausch/lorax).

`rayas call -g <genome.fa> -m <control.bam> <tumor.bam>`

## Simple graph visualization

You can convert the output into a dot graph. Each component represents one templated insertion cluster. Nodes are genomic segments and edges represent the cancer genome structure with edge weights equalling the sequencing read support.

`cut -f 4,9 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`

## Citation

Tobias Rausch, Rene Snajder, Adrien Leger, Milena Simovic, Oliver Stegle, Ewan Birney, Marc Jan Bonder, Aurelie Ernst, Jan O. Korbel.          
Long-read sequencing of diagnosis and post-therapy medulloblastoma reveals complex rearrangement patterns and epigenetic signatures.            
[bioRxiv 2022.02.20.480758](https://doi.org/10.1101/2022.02.20.480758)

License
-------
Rayas is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/rayas/blob/master/LICENSE) file for more details.
