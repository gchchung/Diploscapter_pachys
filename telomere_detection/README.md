# Novel telomere detection algorithm using TideHunter
Below you will find a brief description on how to use the telomere detection algorithm.

## Usage

## Background
The motivation for this algorithm was initially to identify the telomeric repeats of _Diploscapter pachys_ and _Diploscapter coronatus_. Previously published assemblies of these two organisms did not contain contigs that terminate in repeated clusters of TTAGGC (3' end) or begin with clusters GCCTAA (5'end), the canonical nematode telomeric repeats. These two species were subsequently re-sequenced using Nanopore and PacBio. We reasoned that if the genomic reads were not intentionally sheared, telomeres could be captured at the 5' and the 3' ends of reads. Furthermore, if Diploscapter had conventional telomeres maintained by a functional telomerase, these telomeric repeats would satisfy at least three conditions usually observed for conventional telomeres:
1. They would consist of tandemly repeating sequences
2. Clusters of telomeric repeats and their reverse complement would be found at the 3’ and 5’ ends of reads, respectively, and
3. The occupancy pattern of the repeats at the ends of reads would resemble an inverted log-normal cumulative curve, which reflects a naturally log-normal distribution of telomere lengths.

Thus, the telomeric repeat pattern can be found by first identifying the most frequent repeat patterns at the ends of the reads, followed by the elimination of patterns that do not satisfy the 3 conditions above.

## How the algorithm works
Much of our algorithm depends on the TideHunter program written by Yan Gao, Bo Liu, Yadong Wang, and Yi Xing (Gao & al. 2019 Bioinformatics). Our algorithm first takes a collection of long reads (from Oxford Nanopore or PacBio sequencing platforms) and partitions each read into two - the first 1000 bps and the last 1000 bps. TideHunter is run on these partitions to detect any tandem repeats of length ```k```. Next, our algorithm reads the TideHunter output and ranks the repeated sequences together with their reverse complements based on occupancy (saved as *_TideHunter_parser_OUTPUT_CONDENSED.txt). This ranking step also merges repeat sequences that are circularly permuted (repeats of ATCG are considered the same class as repeats of TCGA, CGAT and GATC).

To graph the most frequently occurring repeats, our algorithm runs TideHunter again on full-length reads and parses through the output. Only the top n repeat patterns represented in the first and last 1000 bps of reads are considered. The user provides a length window ```l``` for which occupancies will be graphed - the value of ```l``` should be determined by trial and error to fully cover the lengths of telomeres: we found that 2000 bps were sufficient for _Diploscapter_ telomeres, while _Caenorhabditis_ telomeres required at least 6000 bps. Finally, our algorithm captures the occupancy patterns in the first ```l``` bps, last ```l``` bps and the middle of all reads longer than ```2l```, then graphs these patterns.

Finally, the user should scan the graphed occupancy patterns to determine if the three conditions above have been satisfied. Any _stranded_ occupancy - that is, a repeat pattern occuring only on one end of the reads and its reverse complement only on the other end - should be noted and tested further.

## Requirements
The script is written in Python 3. Required packages include:
- BioPython (I/O of sequence files)
- PIL (for graph png output) 
- svgutils (for graph svg output)
- plotly (for graphing)
- kaleido

In addition, TideHunter v1.4.2 should be installed in the $PATH. Later versions of TideHunter may produce an ouptut whose column orders are not the same.
- TideHunter v1.4.2

## Usage

