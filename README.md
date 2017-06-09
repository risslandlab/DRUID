# DRUID
A pipeline for reproducible transcriptome-wide measurements of mRNA stability

## Full Documentation

See the [Wiki](https://github.com/risslandlab/DRUID/wiki) for full 
documentation, examples, operational details and other information.

## Installation

DRUID is available as a set of scripts and does not need to be compiled or
installed. DRUID does however have some prerequisites.

The main pipeline requires [samtools](https://github.com/samtools/samtools) and 
[bedtools](http://bedtools.readthedocs.io/en/latest/index.html), [python](http://www.python.org/) environment with [NumPy](http://numpy.scipy.org/), a working 
[R](https://cran.r-project.org/) environment with 
[gplots](https://cran.r-project.org/web/packages/gplots/index.html) for optional
plots. 

In addition to the main pipeline, the helper R script 
[processGTF.R](https://github.com/risslandlab/DRUID/blob/master/processGTF.R)
requires the [GenomicFeatures](http://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html), [rtracklayer](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html), and [plyr](https://cran.r-project.org/web/packages/plyr/index.html) libraries.

We have strived to make DRUID independent of any specific RNA-seq pipeline,
however, the [Wiki](https://github.com/risslandlab/DRUID/wiki) outlines our 
standard pipeline based around [STAR](https://github.com/alexdobin/STAR) and 
[HTseq](http://www-huber.embl.de/HTSeq/doc/overview.html).
