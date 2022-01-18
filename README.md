# DRUID
A pipeline for reproducible transcriptome-wide measurements of mRNA stability

## Example usage

>The code below has been simplified as much as possible to focus only on the 
required steps to properly run DRUID. The required input data is described 
verbally, rather than with specific commands, as different users may have 
varying means of generating the same data.

```
source("DRUID.R") # loads function like the DRUID() wrapper function and Filter()

# Gathering required data (for details see the wiki that is linked below) ----

# 1. Table of exonic read counts.
# Row names are genes and column names are  the different samples (time points).
data.reads <- read.csv() # dummy code, needs actual path to the data, of course

# 2. Table of spike-in read counts (not required for DRUID), similar in
#   appearance to the previous table.
spike.reads <- read.csv()

# 3. Vector of library sizes e.g. can be obtained from STAR log files.
library.sizes <- read.csv()

# Read preDruid results = intron counts ---------------------------------------
# hopefully the only csvs present in the directory
csvs <- dir("path/to/intersect", pattern = "csv$", full.names = TRUE)  

# Reorder the files so that they are in increasing time order, there are four
#   files per timepoint, positive and negative strand for both exons and
#   introns. The order of these four files should be consistent for each time-
#   point.
csvs <- csvs[c(5:8, 13:24, 1:4, 9:12, 25:28)] 

# Read in and merge the preDRUID results for introns
introns.negative.strand <- ExtractAndMerge(csvs, 4, 2)
introns.positive.strand <- ExtractAndMerge(csvs, 4, 3)
introns.all <- rbind(introns.negative.strand, introns.positive.strand)

# Calculate half-lives with DRUID
times <- c(1, 2, 4, 8, 12, 24)  # timepoints at which samples were collected
data <- Filter(data.reads) # exonic counts; all entries with zeros removed because nls() cannot handle them. An alternative might be nlsLM()

# preDRUID results will contain information about reads mapping to all 
#   organisms / spike-ins. We need to specify that we want only those
#   corresponding to the organism under study.
introns.human <- introns.all[introns.all$organism == "h", ]

# run the wrapper functiom
half.lives.DRUID <- DRUID(
  x = times,
  introns = introns.human,
  data = data,
  lib.sizes = library.sizes,
  path = "DRUID.replicate.1",
  kmeans.seed = 42
)

# Downstream analysis ...
```

## Full Documentation

See the [Wiki](https://github.com/risslandlab/DRUID/wiki) for full 
documentation, examples, and other information.

## Installation

DRUID is available as a set of scripts and does not need to be compiled or
installed. DRUID does however have some prerequisites.

The main pipeline requires [samtools](https://github.com/samtools/samtools) and 
[bedtools](http://bedtools.readthedocs.io/en/latest/index.html), [python](http://www.python.org/) with [NumPy](http://www.numpy.org/), and a working 
[R](https://cran.r-project.org/) environment with 
[gplots](https://cran.r-project.org/web/packages/gplots/index.html) for optional
plots. 

In addition to the main pipeline, the helper R script 
[processGTF.R](https://github.com/risslandlab/DRUID/blob/master/processGTF.R)
requires the [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html), [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html), and [plyr](https://cran.r-project.org/web/packages/plyr/index.html) libraries.

We have strived to make DRUID independent of any specific RNA-seq pipeline,
however, the [Wiki](https://github.com/risslandlab/DRUID/wiki) outlines our 
standard pipeline based around [STAR](https://github.com/alexdobin/STAR) and 
[HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html).

## Reference

* Lugowski A, Nicholson B, Rissland OS. DRUID: a pipeline for transcriptome-wide measurements of mRNA stability. RNA. 2018 May;24(5):623-632. doi: 10.1261/rna.062877.117. Epub 2018 Feb 8. PMID: 29438994; PMCID: PMC5900561.
