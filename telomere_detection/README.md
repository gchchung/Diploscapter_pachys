#Telomere detection algorithm using TideHunter
Below you will find a brief description on how to use the telomere detection algorithm. 

##How it works
Much of our algorithm depends on the TideHunter program written by Yan Gao, Bo Liu, Yadong Wang, and Yi Xing (Gao & al. 2019 Bioinformatics). Our algorithm first takes a collection of long reads (from Oxford Nanopore or PacBio sequencing platforms) and partitions each read into two - the first 1000 bps and the last 1000 bps. TideHunter is run on these partitions to detect any tandem repeats of length k. Next, our algorithm reads the TideHunter output and ranks the repeated sequences together with their reverse complements based on occupancy (saved as *_TideHunter_parser_OUTPUT_CONDENSED.txt). This ranking step also merges repeat sequences that are circularly permuted (repeats of ATCG are considered the same class as repeats of TCGA, CGAT and GATC).
To graph the most frequently occurring repeats, our algorithm runs TideHunter again on full-length reads and parses through the output. Only the top n repeat patterns represented in the first and last 1000 bps of reads are considered. The user provides a length window l for which occupancies will be graphed - the value of l should be determined by trial and error to fully cover the lengths of telomeres. Finally, our algorithm captures the occupancy patterns in the first l bps, last l bps and the middle of all reads longer than 2l, then graphs these patterns.

##Requirements
The script is written in Python 3. Required packages include:
-BioPython (I/O of sequence files)
-PIL (for graph png output) 
-svgutils (for graph svg output)
-plotly (for graphing)
-kaleido
In addition, TideHunter v1.4.2 should be installed on the $PATH. Later versions of TideHunter may produce an ouptut whose column orders are not the same.
-TideHunter v1.4.2


