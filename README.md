[![C/C++ CI](https://github.com/tobiasrausch/rayas/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/rayas/actions)
[![Docker CI](https://github.com/tobiasrausch/rayas/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/rayas/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/rayas/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/rayas.svg)](https://github.com/tobiasrausch/rayas/releases)

# Rayas

Discovery of templated insertion threads

# Running rayas

Rayas was designed for paired-end or single-end matched tumor-normal data. It can be used with short read or long read data.

`rayas call -g <genome.fa> -m <control.bam> <tumor.bam>`

# Simple graph visualization

You can convert the output into a dot graph. Each component represents one templated insertion cluster. Nodes are genomic segments and edges represent the cancer genome structure with edge weights equalling the sequencing read support.

`cut -f 4,9 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`
