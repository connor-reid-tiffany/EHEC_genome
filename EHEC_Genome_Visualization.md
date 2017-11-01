### Visualize and make a comparison plot with K12. check for prophage and tRNA sequences

Considering Circos is the best option, but has a steep learning curve, I opted to do this in Artemis for now (probably better for people in the class as well). Constructing a circular genome representation in R could also be challenging for a beginner. The software they used in the paper is unfortunately not open source.

First, now that we have our EHEC genome, we need the ecoli K12 genome. Downloading it from either the BLAST web app or command line is easy. I downloaded the genbank file for use in Artemis, and the FASTA file for use in LAST

After playing around with viewing the genome in Artemis (not many gaps! yay!), I used a web application called PHASTER to look for prophage sequences in the genome (they did this in the paper)

Following some prophage gene analysis, I built a genome comparison dotplot using the command line tool LAST. I would recommend that anyone who wants to use LAST should use anaconda to install it. When I compiled it using clang on mac the last-dotplot command was completely busted.

```
conda install last
```
make a true FASTA file by concatenating the contig file of the EHEC genome using a BASH script (though i would also recommend doing an alignment with the contigs fa file, as this will randomly select some and have the effect of giving you a zoomed in view of shared alignment)

from the directory with your final_contigs.fa file
```
grep -v "^>" myecoli.fna | awk 'BEGIN { ORS=""; print ">Sequence_name\n" } { print }' > EHEC.fasta
```

make a directory for our alignments
```
cd
mkdir last_alignments
lastdb -cR01 EHEC_db_2 EHEC.fasta

lastal EHEC_db_2 ~/ecoli_k12.fasta > myalns.maf

last-dotplot --lengths1 --lengths2 -y myalns_2.maf dotplot.png


```
