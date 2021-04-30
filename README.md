# Rayas

Templated insertion discovery

# Running rayas

Rayas was designed for paired tumor-normal data.

`rayas call -g <genome.fa> -m <control.bam> <tumor.bam>`

# Simple graph visualization

You can convert the output into a dot graph. For instance, component 0 can be visualized using

`awk '$6==0' out.bed | cut -f 7 | sed -e '1s/^/graph {\n/' | sed -e '$a}' > out.dot`

`dot -Tps out.dot -o out.pdf`
