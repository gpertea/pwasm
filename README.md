# pafreport
This program analyzes PAF alignments with `cs` tags (generated by `minimap2` with `--cs` option) producing a report about 
the sequence differences found (indels and substitusions) and their sequence context (e.g. if specific motifs
or homopolymers were found at that location).

The primary use case is when a known bacterial gene sequence (CDS) is aligned to one or more genomes (e.g. nanopore assemblies).
In this case the report will include an entry for each assembly detailing the gene sequence coverage and the context for each
putative mutation event, including the downstream consequences of the mutation (aminoacid changes, premature stop codon).

#Building

This project depends on my gclib source. The build steps from the github repository are like this:
    
    git clone https://github.com/gpertea/gclib
    git clone https://github.com/gpertea/pwasm
    cd pwasm
    make release
    
This should build the `pafreport` binary.
