# Rayas

Templated insertion discovery

# Running rayas

Rayas was designed for paired-end matched tumor-normal data.

`rayas call -g <genome.fa> -m <control.bam> <tumor.bam>`

# Simple graph visualization

You can convert the output into a dot graph. Each component represents one templated insertion cluster. Nodes are genomic segments and edges represent the cancer genome structure with edge weights equalling the sequencing read support.

`cut -f 4,9 out.bed | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tpdf out.dot -o out.pdf`
