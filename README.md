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

[placeholder for QC score images]
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

Here, I use a cutadapt, along with the wrapper trim_galore, and use the --paired flag (we have paired end read! but maybe neglect to tell the students which flag to use?) and the nextera flag (because we have nextera adapter contamination, again maybe let them figure this out on their own =])


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

<img src="https://github.com/dib-lab/sourmash/blob/master/doc/_static/cmp.matrix.png?raw=true" style="width:60%" />
