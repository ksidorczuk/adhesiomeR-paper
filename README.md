# adhesiomeR-paper

### 01-Filter_RefSeq.R

This script was used to filter downloaded genomes based on associated metadata.

Data files necessary to run the scripts are available [here](https://www.dropbox.com/scl/fo/gc2mz0n4d66c65gxgl0ym/h?rlkey=h1446yb8hopv47axag0ankv8o&dl=0).

### 02-Adhesin_sequence_comparison.R

Analysis and generation of plots from blast all to all comparison of adhesin sequences. 

Files with adhesin information necessary to run the script: [Adhesins.xlsx](https://www.dropbox.com/scl/fi/mac3wzohoxu0wpjloj2wy/Adhesins.xlsx?rlkey=dbjjmofbtvw8h9hjxmixc1mhq&dl=0)
Directory with BLAST results and plots we obtained: [BLAST](https://www.dropbox.com/scl/fo/s2lma468wqlcdmgrjjxmm/h?rlkey=hhd385dlpmjy6t3wpxevr7897&dl=0)

### 03-Initial_pathotypes_analysis.R

Initial analysis of pathotyped genomes with older version of adhesiomeR (commit 0f675a9), which should correspond to current relaxed version with thresholds of identity and coverage set to 80%. 
You can download result files [here](https://www.dropbox.com/scl/fo/fiqua9rhrnsg1wve1y0ww/h?rlkey=hxvgu9lnu5ookbl8jfxwze3fm&dl=0)

### 04-Fimbrial_adhesins_analysis.R

This script describes analysis of fimbrial adhesins and their gene co-localization to select the reference set for calculation of bit score thresholds.
It uses results linked in the previous step. 
You can download all results obtained by us [here](https://www.dropbox.com/scl/fo/vxyiivnv000wvb0g6d6r2/h?rlkey=1t18zyn72k5tsu4zdq67vf814&dl=0). Note that this directory contains a lot of large files. 

### 04-Intimin_analysis.R

Analysis of intimin variants to select one representative sequence for the database. This script runs BLAST on pathotyped genome collection. You can download our results and source sequences for intimin [here](https://www.dropbox.com/scl/fo/cz3hqlp6muvlcaan6byur/h?rlkey=s06odj1z1ylt7qieiro1047ws&dl=0).

### 04-Nonfimbrial_adhesins_analysis.R

Analysis of genomic context of nonfimbrial adhesins to select sequences for the reference set. The analyses were run multiple times, you can download all our results [here](https://www.dropbox.com/scl/fo/1jq9jjx5zhe8oqev7dowx/h?rlkey=vblhi6q2yuiszxquhiv3wrh9o&dl=0). Note that this directory contains a lot of large files. 

### 05-Bitscore_thresholds.R

Filter results from step 4 to obtain reference sequences and calculate bit score thresholds (Table S5). Files necessary to run this script are linked in previous steps. 

### 06-Pathotypes_analysis.R

Analysis of pathotyped genomes with final version of adhesiomeR using both types of search settings.
You can download our results for [strict version](https://www.dropbox.com/scl/fo/t4ebnz5y04q970q26x32l/h?rlkey=t2xiyijnyuukipraqs6vk3j2x&dl=0) and [relaxed version](https://www.dropbox.com/scl/fo/2l2km8ibpcaj4d5bgul6w/h?rlkey=e90e7ixxbzkyl9nway1f3w1wx&dl=0).

### 07-Selecting_k_clara.R

Script running selection of the optimal number of clusters according to gap statistic for all adhesins. 

### 07-Selecting_k_clara_fimbrial.R

Script running selection of the optimal number of clusters according to gap statistic for fimbrial adhesins. 

### 07-Selecting_k_clara_nonfimbrial.R

Script running selection of the optimal number of clusters according to gap statistic for nonfimbrial adhesins. 

### 08-Profile_analysis_and_initial_clustering.R

This script generates adhesin profiles and performs initial clustering according to selected in a previous step optimal numbers of clusters. It generates HTML reports with overview of each clustering (based on **Clara_clustering.Rmd**). 
Clustering data generated in this script is available for download [here](https://www.dropbox.com/scl/fi/41kynbkilc2sv2smafb8r/clustering.RData?rlkey=efghp64murz1yxptvrlkum4jr&dl=0).

### 09-Clustering.R

Performs final clustering, assigns names to clusters, calculates gene importance for each clustering and generates clustering plots with fractions of pathotypes and numbers of profiles and genomes assigned to each cluster. 
File with clustering data used in this script is obtained and linked in a previous step. 

### 10-Benchmark.R

Performs comparison of adhesiomeR results with experimental annotations from Von Mentzer et al. on ETEC strains and generates figures. 
Files necessary to run this analysis: [ETEC metadata](https://www.dropbox.com/scl/fo/53trher9qkozmys4sk3oe/h?rlkey=7j94uxib0ds1enlb0wl1hxw9x&dl=0),
[adhesiomeR results](https://www.dropbox.com/scl/fo/iiw93tqacriv7u1qsgfqz/h?rlkey=j9ovagvsxaeawxwzizq4qktgi&dl=0).
