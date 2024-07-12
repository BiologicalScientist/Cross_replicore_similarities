# Cross replicore similarities
These scripts can be used to identify regions of increased similarity that are either quasisymmetrically - defined as + or - 30kb (Cross_genome_searcher.py) or dispersed throughout the genome (Cross_genome_searcher_random_connect_window_mode.py). 

Two python scripts are provided, one for either analysis.
Python>=3.9 is required along with installation of [cd-hit](https://github.com/weizhongli/cdhit), [seqkit](https://github.com/shenwei356/seqkit), and [Circlator](https://github.com/sanger-pathogens/circlator).


help message and Flags are as follows
```
usage: Cross_genome_searcher.py [-h] -g file.fasta [file.fasta ...] -o OUTPUT_FOLDER [-l LENGTH] [-i]

script to assess whether genomes have repeat regions symmetrically (+-30kb) across the genome

options:
  -h, --help            show this help message and exit
  -g file.fasta [file.fasta ...], --genomes file.fasta [file.fasta ...]
                        a path to input genomes
  -o OUTPUT_FOLDER, --output OUTPUT_FOLDER
                        Path to output folder
  -l LENGTH, --length LENGTH
                        length of the sliding window (in bp) - default 3000bp
  -i, --inverted_repeats
                        restricts search to ONLY consider inverted repeats
```
## optional flags
### -i / --inverted_repeats
If only wishing to consider inverted repeats in your search use the -i flag (default is to consider both inverted and direct repeats due to the dynamic nature of the genome and the potential for inversions within a replichore to change strands of sequences.)
### -l / --length
controls the size of the sliding window in the seqkit sliding window step. Default is 3000bp though smaller values may be desired in some use cases. Note: as the size decreases the number of comparisons increases so expect slower run times for both scripts but particularly asymmetric comparisons.


## Troubleshooting:
Cross_genome_searcher_random_connect_window_mode.py is considering a much larger search space than the symmetric plot as it features an "all vs all" style comparison. As a result this runs much slower than the symmetric script. 


## post-processing visualisation

The R script provided (Shadow_plot_script.R) outlines how to create the shadow plot (see figure 1)

<!--- Add in image and caption --->
<p>
<p align = "center">
<img src = "illustrations/S_pyogenes_shadow_plot.png">
</p>

<p>
<p align = "center">
Figure 1: Shadow plot for for 249 complete <i>Streptococcus pyogenes</i> genomes. Left is the shadow plot for quasisymmetrically placed regions. Right is the shadow plot for regions dispursed thorughout the genome with connection between reigon from opposing halfs of the genome. Produced using Shadow_plot_script.R and the python script given beneth each shadow plot. 
</p>
<!--- _______________________ --->

An example of an aoutput file used to produce the left shadown plot of figure 1 is given (Streptococcus_pyogenes_relative_genome_connections_example.tsv).