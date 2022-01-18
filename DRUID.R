################################################################################
#                              Global settings                                 #
################################################################################
require("gplots")  # heatmap.2

# To ensure consistent file ordering between linux and windows.
Sys.setlocale("LC_COLLATE", "C") 

################################################################################
#                                  FUNCTIONS                                   #
################################################################################

Filter <- function(gene.matrix, all = 1, any = 5) {
  # Removes rows from a numeric matrix of gene abundances failing to meet
  #   minimum expression requirements.
  #
  # Args:
  #   gene.matrix: A numeric matrix of gene abudances.
  #   all: the minimum gene abundance across all conditions. Default is 1.
  #   any: a minimum abudance that must be observed in at least one condition.
  #     Default is 5.
  #
  # Returns:
  #   The the filtered matrix as a data frame.
  gene.matrix <- 
    gene.matrix[
      apply(
        gene.matrix,
        1,
        function(row) {
          all(row >= 1) & any(row >= 5) 
        }
      )
      ,
    ];
  return( as.data.frame(gene.matrix))
}

HalfLives4SU <- function(x, y, doubling.time = 24 ) {
  # Computes the half-life of a gene in an approach-to-equilibrium experiment.
  #
  # Args:
  #   x: Vector of experimental time points.
  #   y: Vector of normalized gene abundances.
  #   doubling.time: The cell cycle length of the system under study, used for
  #     correcting determined half-life. No cell cycle correction is performed
  #     if set to 0. Default is 0.
  #
  # Returns:
  #   The determined half-life in the same units used for the time vector x. NA
  #     if the algorithm fails to converge or the determined half-life is
  #     negative.
  suppressWarnings(
    fit <- nls(
      y ~ C * (1 - exp(k * x)),
      start = c(C = max(y), k = -0.25),
      algorithm = "port",
      lower = c(C = 0,k = -Inf), upper = c(C = Inf, k = 0),
      control = list(warnOnly = TRUE),
      weights = 1/y
    )
  )
  
  # check to see if it converged, return if it didn't
  if (!fit$convInfo$isConv) {
    return(NA)
  }
  
  k <- coef(fit)['k']
  C <- coef(fit)['C']  
  hl <- NA

  # calculate half-life performing cell cycle correction if required
  if (doubling.time) {
    hl <- log(2) / -(k + log(2) / doubling.time) 
  } else {
    hl <- log(2) / -(k) 
  }
  
  if (hl < 0) {  # we don't want any negative half-lives
    hl <- NA
  }
  
  return(hl)
}


HalfLivesBlock <- function(x, y, doubling.time = 0) {
  # Computes the half-life of a gene in a transcription inhibition experiment.
  #
  # Args:
  #   x: Vector of experimental time points.
  #   y: Vector of normalized gene abundances.
  #   doubling.time: The cell cycle length of the system under study, used for
  #     correcting determined half-life. No cell cycle correction is performed
  #     if set to 0. Default is 0.
  #
  # Returns:
  #   The determined half-life in the same units used for the time vector x. NA
  #     if the algorithm fails to converge or the determined half-life is
  #     negative.
  suppressWarnings(
    fit <- nls(
      y ~ C * (exp(k * x)),
      start = c(C = max(y), k = -0.25),
      algorithm = "port",
      lower = c(C = 0, k = -Inf), upper = c(C = Inf, k = 0),
      control = list(warnOnly = TRUE)
    )
  )
  
  # check to see if it converged, return if it didn't
  if( !fit$convInfo$isConv ) {
    return( NA )
  }
  
  k <- coef(fit)['k']
  C <- coef(fit)['C']
  hl <- NA

  # calculate half-life performing cell cycle correction if required
  if (doubling.time) {
    hl <- log(2) / -(k + log(2) / doubling.time) 
  } else {
    hl <- log(2) / -(k) 
  }

  if (hl < 0) {  # we don't want any negative half-lives
    hl <- NA
  }
  
  return(hl)
}

HalfLives <- function(x, data, norm, func, doubling.time = 24) {
  # Computes the half lives of multiple genes.
  #
  # Args:
  #   x: Vector of experimental time points.
  #   data: matrix of unnormalized gene abundances.
  #   norm: vector of  
  #   func: the function to use for calculating individual half-lives. Currently
  #     choices are HalfLives4SU for approach-to-equilibirum experiments, and 
  #     HalfLivesBlock for transcription inhibition experiments.
  #   doubling.time: The cell cycle length of the system under study, used for
  #     correcting determined half-lives. No cell cycle correction is performed
  #     if set to 0. Default is 0.
  #
  # Returns:
  #   Vector of determined half-life in the same units used for the time vector.
  half.lives <- apply(
    data, 1,
    function(y, x) {
      func(x, y / norm, doubling.time)
    },
    x
  )
  return(half.lives)
}

LeaveOneOut <- function(x, data, norm, func, doubling.time = 24) {
  # Computes the half lives of multiple genes using all timepoints and again 
  #   leaving out individual time points.
  #
  # Args:
  #   x: Vector of experimental time points.
  #   data: Matrix of unnormalized gene abundances.
  #   norm: Vector of values for normalizing gene abundance.
  #   func: The function to use for calculating individual half-lives. Currently
  #     choices are HalfLives4SU for approach-to-equilibirum experiments, and 
  #     HalfLivesBlock for transcription inhibition experiments.
  #   doubling.time: The cell cycle length of the system under study, used for
  #     correcting determined half-lives. No cell cycle correction is performed
  #     if set to 0. Default is 0.
  #
  # Returns:
  #   Matrix of determined half-life in the same units used for the time vector.
  hl <- HalfLives(x, data, norm, func, doubling.time)
  loo <- cbind(
    hl,
    sapply(
      1:length(x),
      function( leave.out ) {
        HalfLives(
          x[-leave.out],
          data[, -leave.out],
          norm[-leave.out],
          func,
          doubling.time
        )
      }
    )
  )
  row.names(loo) <- row.names(data)
  colnames(loo) <- c(
    "std",
    sapply(
      x,
      function(x) {
        paste("no_", x, sep = "")
      }
    )
  )
  return(loo)
}

ExtractAndMerge <- function(
  csvs, jump = 4, offset = 0,
  suffixes=c(".1", ".2", ".4", ".8", ".12", ".24", ".SS") 
) {
  # Merges intersection tables generated by preDRUID.sh script into a single 
  #   data frame.
  #
  # Args:
  #   csvs: A list of paths to the intersection files. Should be sufficient to
  #     list the files in the intersect directory created by preDruid.sh and
  #     potentialy reorder if default naming order is insufficient.  
  #   jump: The number of files for every data point Default is 4 to accomodate
  #     positive and negative strand for both exons and introns.
  #   offset: Sets the specific file type to merge, e.g. positive strand exons. 
  #     This function assumes that the four possible file types will be
  #     consistently ordered for all data points.
  #   suffixes: Labels that will be appended to shared column names.
  #
  # Returns:
  #   Data frame of merged intersection tables usually for a single subtype of 
  #     tables, e.g positive strand exonic coverage. Reported properties include
  #     mean base coverage and associated standard deviation, feature size, and 
  #     proportion of bases coverage by at least one read.
  files <- lapply(csvs, read.csv)
  table.names <- colnames(files[[1]])
  a <- files[[1 + offset]]
  colnames(a)[-1:-5] <- sapply(
    colnames(a)[-1:-5], paste, suffixes[1], sep = "", collapse = ""
  )
  for (i in 1:(length(suffixes) - 1)) {
    index <- i * jump + offset + 1
    b <- files[[index]]
    colnames(b)[-1:-5] <- sapply(
      colnames(b)[-1:-5], paste, suffixes[i + 1], sep="", collapse=""
    )
    a <- merge(a, b, by = c(table.names[1:5]), all = TRUE)
  }
  return(a)
}

DRUID <- function(x, introns, data, lib.sizes, path, filterIntrons = TRUE,
  plot = TRUE, minCutoff = 0.5, kmeans.seed = 0) {
  # Determines half-lives for an approach-to-equilibium experiment by 
  #   normalizing to a subset of well-behaved intronic reads. Performs
  #   leave-one-out analysis.
  #
  # Args:
  #   x: Vector of experimental time points.
  #   introns: Data frame of intronic coverage as generated by ExtractAndMerge
  #     from preDruid.sh results.
  #   data: Matrix of unnormalized gene abundances, generally the same table
  #     that would be used for normal half-life calculations.
  #   lib.sizes: Vector of library sizes.
  #   path: Path prefix for saving results including introns used for
  #     normalization as well as heatmap of introns passing abundance cutoff.
  #   filterIntrons: If TRUE, only well-behaved introns will be used for 
  #     normalization. Default is TRUE.
  #   plot: If TRUE, a heatmap of introns passing abundance cutoff will be
  #     plotted. Default is TRUE.
  #   minCutoff: The minimum mean base coverage of an intron to be included in
  #     analysis. Default is 0.5.
  #   kmeans.seed: Set the kmeans.seed if results need to be reproducible, if
  #     set to zero, default seed will be used. Default is 0.
  #
  # Returns:
  #   Matrix of determined half-life in the same units used for the time vector.
  indices <- 1:length(x)
  
  # Isolate and filter intron mean coverage.
  introns.mean <- introns[, grep("mean", colnames(introns))]
  introns.mean <- introns.mean[, indices]
  row.names(introns.mean) <- paste(introns$gene, '.', introns$exon, sep = "")
  introns.mean <- introns.mean[complete.cases(introns.mean), ]
  introns.mean.filt <- introns.mean[
    apply(introns.mean, 1, function(x) {
      all(x >= minCutoff)
    }),
  ]
  
  # Normalize to library size.
  lib.sizes.matrix <- matrix(
    rep(lib.sizes, nrow(introns.mean.filt)),
    byrow = TRUE,
    ncol = length(lib.sizes)
  )
  introns.mean.filt.norm <- introns.mean.filt / lib.sizes.matrix
  maxes <- apply(introns.mean.filt.norm, 1, max)
  maxes.matrix <- matrix(rep(maxes, length(maxes)), ncol = length(maxes))
  introns.mean.filt.norm <- introns.mean.filt.norm / maxes.matrix
  
  # Plot heatmap.
  if (plot) {
    pdf(paste(path, ".heatmap.DRUID.pdf", sep = ""), pointsize = 8)
    par(lty = 1)
    heatmap.2(
      as.matrix(introns.mean.filt.norm),
      trace = "none",
      dendrogram = "row", labRow = NA,
      Colv = NA,
      notecol = "black", notecex = 0.75,
      density.info = 'none',
      breaks = seq(0.0, 1.0, length.out = 101)
    )
    dev.off()
  }
  
  # Cluster
  if (kmeans.seed) {
    set.seed(kmeans.seed)
  }
  fit <- kmeans(introns.mean.filt.norm, 4)
  
  # Pick the right cluster (max average first timepoint).
  # A more robust method could definetly be implemented, this approach will
  #   likely fail if very short labeling periods (<5 min) are included and even
  #   well-behaved introns have not yet reached their steady-state levels.
  mean.first.timepoint <- aggregate(
    introns.mean.filt.norm[, 1],
    by = list(fit$cluster),
    mean
  )
  desired.cluster <- mean.first.timepoint$Group.1[
    mean.first.timepoint$x == max(mean.first.timepoint$x)
  ]
  stopifnot(length(desired.cluster) == 1)
  
  # Get the correct intron coverages and sizes.
  desired.exons <- row.names(introns.mean.filt)[fit$cluster == desired.cluster]
  introns.size <- introns[, grep("size", colnames(introns))]
  introns.size <- introns.size[, indices]
  row.names(introns.size) <- paste(introns$gene, ".", introns$exon, sep = "")
  introns.size.filt <- introns.size[desired.exons, ]
  introns.mean.filt <- introns.mean.filt[desired.exons, ]
  
  # Compute normalization vector.
  norm <- colSums(introns.mean.filt * introns.size.filt)
  if (filterIntrons) {
    write.csv(
      introns.mean.filt * introns.size.filt,
      file = paste(path, ".filterIntrons.norm.csv", sep = "")
    )
  } else {
    norm <- colSums(introns.mean * introns.size, na.rm = TRUE)
  }
  
  # Now determine half-lives as per usual, note cell-cycle compensation is not
  #   performed.
  DRUID.HL <- LeaveOneOut(
    x, data, norm, HalfLives4SU, 0
  )
}

################################################################################
#                                   Example                                    #
################################################################################

# The code below has been simplified as much as possible to focus only on the 
#   required steps to properly run DRUID. The required input data is described 
#   verbally, rather than with specific commands, as different users may have 
#   varying means of generating the same data.


# Gathering required data.

# Table of exonic read counts. Row names are genes and column names are
#   the different samples (time points).
data.reads <- read.csv()

# Table of spike-in read counts (not required for DRUID), similar in
#   appearance to the previous table.
spike.reads <- read.csv()

# Vector of library sizes e.g. can be obtained from STAR log files.
library.sizes <- read.csv()

# Read preDruid results, hopefully the only csvs present in the directory
csvs <- dir("path/to/intersect", pattern = "csv$", full.names = TRUE)  
# Reorder the files so that they are in increasing time order, there are four
#   files per timepoint, positive and negative strand for both exons and
#   introns. The order of these four files should be consistent for each time-
#   point.
csvs <- csvs[c(5:8, 13:24, 1:4, 9:12, 25:28)] 

# We will read in preDRUID results for exonic features for completeness, but we
#   wont actually be using them.
exons.negative.strand <- ExtractAndMerge(csvs, 4, 0)
exons.positive.strand <- ExtractAndMerge(csvs, 4, 1)
exons.all <- rbind(exons.negative.strand, exons.positive.strand)

# Read in and merge the preDRUID results for introns
introns.negative.strand <- ExtractAndMerge(csvs, 4, 2)
introns.positive.strand <- ExtractAndMerge(csvs, 4, 3)
introns.all <- rbind(introns.negative.strand, introns.positive.strand)

# Calculate half-lives

# Standard spike-in method
times <- c(1, 2, 4, 8, 12, 24)  # timepoints at which samples were collected
data <- data.reads
norm <- colSums(spike.reads)

half.lives.spikeins <- LeaveOneOut(times, data, norm, HalfLives4SU)

# DRUID
times <- c(1, 2, 4, 8, 12, 24)  # timepoints at which samples were collected
data <- data.reads
# preDRUID results will contain information about reads mapping to all 
#   organisms / spike-ins. We need to specify that we want only those
#   corresponding to the organism under study.
introns.human <- introns.all[introns.all$organism == "h", ]
half.lives.DRUID <- DRUID(
  x = times,
  introns = introns.human,
  data = data,
  lib.sizes = library.sizes,
  path = "DRUID.replicate.1",
  kmeans.seed = 42
)

# Downstream analysis ...
