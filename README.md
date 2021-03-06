# pafreport
This program analyzes PAF alignments with `cs` tags (generated by `minimap2` with `--cs` option) producing a report about 
the sequence differences found (indels and substitutions) and their sequence context (e.g. if specific motifs
or homopolymers were found at that location).

The primary use case is when a known bacterial gene sequence (CDS) is aligned to one or more genomes (e.g. nanopore assemblies).
In this case the report will include an entry for each assembly detailing the gene sequence coverage and the context for each
putative mutation event, including the downstream consequences of the mutation (aminoacid changes, premature stop codon).

# Building

This project depends on my gclib source. The build steps from the github repository are like this:
    
    git clone https://github.com/gpertea/gclib
    git clone https://github.com/gpertea/pwasm
    cd pwasm
    make release
    
This should build the `pafreport` binary.


# Usage
For aligning a reference CDS (bacterial gene sequence, cds.fa) vs. a set of Nanopore assemblies (asms.fa):
    
    minimap2 -x map-ont -d asms.mmi asms.fa #build the minimizer index (optional)
    minimap2 -c --cs -r1000 -P -x map-ont asms.mmi cds.fa > mappings.paf

Then pafreport can be run on the resulting PAF file:

    pafreport -r cds.fa mappings.paf -o report.txt

