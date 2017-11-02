# Reassembly of the Escherichia coli O157:H7 genome

## Recreating the results of a classic paper

Here, I will be going through the assembly and annotation of Escherichia coli O157:H7,
and comparing the results to the first complete assembled genome of an EHEC strain
[(you can find the paper here!)][0]

### Step 1, download the SRA toolkit and raw reads

First we will need the toolkit
```
wget "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
```

then unzip
```
tar -xzf sratoolkit.current-centos_linux64.tar.gz
```
then make it useable from any directory
```
export PATH=$PATH:$PWD/sratoolkit.2.8.2-1-mac64/bin
```

Now that the toolkit is installed, we can download and split our raw reads using the appropriate accession number and the `fastq-dump` command with the split files flag

First make a new directory from your home directory:

```
mkdir EHEC_genome
```
then download and split the files
```
fastq-dump --split-files SRR4301589
```
check to make sure the files are there using `ls`

### Step 2, Assess the quality of your reads

Note: I had trouble figuring out how to run FastQC from any directory, and was unable to custom set a path (probably security junk since i did this on an MMI machine), so I just ran the fastqc function from inside the FastQC folder and used the full path for my files in the arguments (clunky, I know)
```
./fastqc ~/EHEC_genome/seqs/SRR4301589_1.fastq ~/EHEC_genome/seqs/SRR4301589_2
```
then I moved into the directory with my QC info files and opened them.

```
cd ~/EHEC_genome/seqs
open SRR4301589_1_fastqc.html
open SRR4301589_2_fastqc.html

```
Looking at the readout, I could see per base quality was low and adapters were present, specifically Nextera transposase sequences


`spoiler alert`: It seems there are some Nextera Transposase Adapters contaminating our sequences, but there is no mention of which sequences were used in this experiment on the SRA....so we need to figure out what the adapter sequences are, trim them, and then reassess the quality of our reads.

### Step 3: Trim adapters and bad stuff
```
TrimmomaticPE SRR1976948_1.fastq.gz \
                 SRR1976948_2.fastq.gz \
        SRR4301589_1.fastq  s1_se \
        SRR4301589_2.fastq  s2_se \
        ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
        LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 \
        MINLEN:25
```
Trimmomatic didn't actually trim any adapters (or sequence for that matter). So lets try another QC trimmer!

Here, I use a cutadapt, along with the wrapper trim_galore, and use the --paired flag (because they are paired end reads!) and the --nextera flag (because from the QC analysis showed the adapter contamination is from nextera transposase sequences)


```
./trim_galore --paired --nextera --quality 25 --length 150 -o EHEC_genome/phred_25_threshold_150min/  SRR4301589_1.fastq SRR4301589_2.fastq

```

it worked! you can now see that the nextera transposase sequences are no longer there after re running fastqc on the trimmed files.

```
./fastqc ~/EHEC_genome/seqs/SRR4301589_1_val_1.fq ~/EHEC_genome/seqs/SRR4301589_2_val_2.fq

open SRR4301589_1_val_1_fastqc.html
open SRR4301589_2_val_2_fastqc.html
```

you can find the QC figures [here][1]
### Step 4: Assemble the genome

For Assembly, I decided to go with Megahit because word on the street was that it was a great assembler for Illumina reads if computational power is limited, and both of those criteria apply here!

Note: A common theme I have been seeing is that anyone on a Mac should really just use anaconda to install tools when possible, as Xcodes garbage g++ imitator sucks at compiling, and it does a great job at setting paths for you.

So 3 cores and 6GB of RAM is pretty self explanatory, since i have 4 cores and 8GB of RAM on this laptop (guessing -m stands for memory and -t stands for threads). Letting ppl crash their machines could be both a learning experience and a common computational biology occupational hazard though.

```
megahit -1 SRR4301589_1_val_1.fq  -2 SRR4301589_2_val_2.fq  -m 0.75  -t 3  -o megahit_result
```

Okay so now I have contigs, and after a nice `head | less` (hah) combo, I can see its still DNA, but who knows if this assembly is garbage or not. Maybe (?) QUAST knows? Even though Daniel recently told me its kind of garbage i dont know of a reasonable alternative (if there even is one).

after running QUAST:

```
cd megahit_result/
~/quast/quast.py final.contigs.fa -o ecoli_report
```


I got this output table. The contigs and total genome length seem a bit high for an e.coli (maybe like half a million bases over)

```
Assembly                    final.contigs
# contigs (>= 0 bp)         1392         
# contigs (>= 1000 bp)      200          
# contigs (>= 5000 bp)      84           
# contigs (>= 10000 bp)     67           
# contigs (>= 25000 bp)     51           
# contigs (>= 50000 bp)     28           
Total length (>= 0 bp)      5703429      
Total length (>= 1000 bp)   5284784      
Total length (>= 5000 bp)   5037663      
Total length (>= 10000 bp)  4918477      
Total length (>= 25000 bp)  4659567      
Total length (>= 50000 bp)  3866983      
# contigs                   327          
Largest contig              443851       
Total length                5370548      
GC (%)                      50.30        
N50                         108305       
N75                         43069        
L50                         14           
L75                         32           
# N's per 100 kbp           0.00
```
So maybe I should have trimmed more? Or Less?

Upon consulting with the lovely article linked in the ANGUS trimming tutorial, I decided to go back and trim with an additional `-q 4` flag to eliminate bases of <phred 4. It seems trim_galore defaults at 20, so maybe less is more in this case and i will get a less weird assembly with a lower QC threshold, or maybe the assembly is actually fine.

Okay so doing assemblies at <25 and <4 gave comparable results to <20, and on the SRA it says the genome they got was 5.3 MB (same as mine), so assembly is as good as it gets. heres the trim code and the quast readout from the assembly


```
megahit -1 SRR4301589_1_val_1.fq  -2 SRR4301589_2_val_2.fq  -m 0.75  -t 3  -o megahit_result

~/quast/quast.py final.contigs.fa -o ecoli_report_original
```

```
All statistics are based on contigs of size >= 500 bp, unless otherwise noted (e.g., "# contigs (>= 0 bp)" and "Total length (>= 0 bp)" include all contigs).

Assembly                    final.contigs
# contigs (>= 0 bp)         882          
# contigs (>= 1000 bp)      203          
# contigs (>= 5000 bp)      86           
# contigs (>= 10000 bp)     68           
# contigs (>= 25000 bp)     53           
# contigs (>= 50000 bp)     30           
Total length (>= 0 bp)      5538988      
Total length (>= 1000 bp)   5297885      
Total length (>= 5000 bp)   5045853      
Total length (>= 10000 bp)  4918550      
Total length (>= 25000 bp)  4682936      
Total length (>= 50000 bp)  3879205      
# contigs                   293          
Largest contig              443875       
Total length                5357863      
GC (%)                      50.29        
N50                         98587        
N75                         43069        
L50                         15           
L75                         34           
# N's per 100 kbp           0.00
```
### Annotate with prokka
First inside of the EHEC_genome directory, make a new directory for the annotation data, then link this directory to the directory with the assembly.

```
mkdir EHEC_Annotation/
cd EHEC_Annotation/
ln -fs ~/EHEC_genome/EHEC_Annotation/megahit_result/final.contigs.fa
```

Now run prokka on the assembled contigs
```
prokka final.contigs.fa --outdir prokka_annotation --prefix myecoli
```
You should get several output files. Now the genome is annotated! The next step is to  take these output files and start doing cool things like genome visualization and prophage visualization.

### Visualizing the genome

Here I go from using mostly command line tools to mostly GUI tools in the interest of time and sanity (Circos is cool but the learning curve is painful).

First I used [Artemis][2] to visualize the genome. By itself and comparatively with the ecoli_k12 genome (using GFF and GBK files respectively)

<img src="https://github.com/recursive-deletion/EHEC_genome/blob/master/figures/EHEC_circular_fig.png" style="width:60%" />

Then, I viewed the EHEC genome alongside the ecoli_k12 genome. You can see the sequence difference is about 700Kb, why is that? (hint, its the prophage DNA!)

<img src="https://github.com/recursive-deletion/EHEC_genome/blob/master/figures/EHEC_vs_K12.png" style="width:60%" />

Speaking of prophage DNA, I wanted to get a closer look at where my potential phage sequences were on the genome. Fortunately theres a web application called [PHASTER][3] that does a decent job.

<img
src="https://github.com/recursive-deletion/EHEC_genome/blob/master/figures/circularview_prophage.png" style="width:60%" />

### Using LAST to do multiple genome alignments

To create a visual showing overall conservation between EHEC and K12, I used a command line tool called LAST to generate an alignment and a subsequent dotplot. (A more practical use of this program would probably be to align concatenated protein files for tree building)

unfortunately, compiling it using clang left me with a broken tool, so I had to make a virtual environment using anaconda, and then install and compile LAST in there along with all the dependencies.

```
conda create last_stuff
source activate last_stuff

conda install last
```

I then did alignments using the final.contigs.fa and a EHEC.fasta file I created to the ecoli_k12 genome FASTA file

Create an EHEC.fasta file
```
grep -v "^>" myecoli.fna | awk 'BEGIN { ORS=""; print ">Sequence_name\n" } { print }' > EHEC.fasta
```
Perform the alignments

```
cd
mkdir last_alignments
lastdb -cR01 EHEC_db_2 EHEC.fasta

lastal EHEC_db_2 ~/ecoli_k12.fasta > myalns.maf
```
in the first step, i made a new directory and then made my EHEC database for the alignments, then suqsequently aligned the ecoli_k12 genome to them.

Then take the .maf file and create some dotplots

```
last-dotplot --lengths1 --lengths2 -y myalns.maf dotplot.png
```

this should give you a figure like this if you used EHEC.fasta

<img
src="https://github.com/recursive-deletion/EHEC_genome/blob/master/figures/dotplot2.png" style="width:30%" />

or this if you used final_contigs.fa

<img
src="https://github.com/recursive-deletion/EHEC_genome/blob/master/figures/dotplot.png" style="width:60%" />


### Future directions

The paper had more detailed information on genomic features, mostly ORFs and prophage related stuff, as well as a nice codon usage table. The best way to approach this would probably involve some programming in python or perl.


[0]:https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/dnaresearch/8/1/10.1093_dnares_8.1.11/3/dnares_8_1_11.pdf?Expires=1509649205&Signature=Vbc9g22UA95lDR6HCiIN2EYy8k~DyvuJwbWJ4UI8swOaJwA~heS2Y20ti24RKqocf0V3Ziv1cP~5Baz0Q5oktT~mcs4TOz6GbiI-1gA97SpCBmPRW0GA8BM9pWarMRqTBt~jIMzeW-Dyiqe7c8UKy~xh~Q5qqmQvpE261SyoHvkKUVASh99rVDTcBeEwjHdeFLVBNSVVnplFQsfbZGUiNqBL2BMfAS5jZFiTnQ712F0NC~CzNHcJUmoKkInBY4LUxHXEUrZ0AUbB6RO9wF47BC4Nitb7ENPVVtP8u~4CuiU3BvBu~O9cblKTP2lCWU4~Ck-Bxnm7BIuGwXHa1q4PdQ__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q
[1]:https://github.com/recursive-deletion/EHEC_genome/tree/master/QC_Analysis
[2]:http://www.sanger.ac.uk/science/tools/artemis
[3]:http://phaster.ca/
