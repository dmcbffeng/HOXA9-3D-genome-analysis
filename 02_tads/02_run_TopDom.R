# TopDom (v0.0.2) TAD caller driver for the CellLinesAnd4943 pipeline.
#
# Reads samples.tsv (sample_name<TAB>hic_path) and runs TopDom on every
# temp_matrix/<sample>_<chr>.matrix file produced by 01_process_matrix.py.
# Per-chromosome outputs are written to results/<sample>_<chr>.{bed,binSignal,domain}.
#
# Usage:
#   Rscript 02_run_TopDom.R                      # all samples / autosomes + chrX
#   Rscript 02_run_TopDom.R --samples 4943 MOLM13
#   Rscript 02_run_TopDom.R --window-size 5
#
# Original TopDom code attribution:
#   @author : Hanjun Shin(shanjun@usc.edu)
#   @credit : Harris Lazaris (NYU), Dr. Gangqing Hu (NIH)

# ---------------- TopDom v0.0.2 (unchanged from previous workflow) -----------

TopDom <- function( matrix.file, window.size, outFile=NULL, statFilter=T)
{
  window.size = as.numeric(window.size)
  matdf = read.table(matrix.file, header=F)

  if( ncol(matdf) - nrow(matdf) == 3) {
    colnames(matdf) <- c("chr", "from.coord", "to.coord")
  } else if( ncol(matdf) - nrow(matdf) ==4 ) {
    colnames(matdf) <- c("id", "chr", "from.coord", "to.coord")
  } else {
    print("Unknown Type of matrix file"); return(0)
  }
  n_bins = nrow(matdf)
  mean.cf <- rep(0, n_bins)
  pvalue <- rep(1, n_bins)
  local.ext = rep(-0.5, n_bins)

  bins <- data.frame(id=1:n_bins, chr=matdf[, "chr"],
                     from.coord=matdf[, "from.coord"],
                     to.coord=matdf[, "to.coord"])
  matrix.data <- as.matrix( matdf[, (ncol(matdf) - nrow(matdf)+1 ):ncol(matdf)] )

  for(i in 1:n_bins) {
    diamond = Get.Diamond.Matrix(mat.data=matrix.data, i=i, size=window.size)
    mean.cf[i] = mean(diamond)
  }

  gap.idx = Which.Gap.Region2(matrix.data=matrix.data, window.size)
  proc.regions = Which.process.region(rmv.idx=gap.idx, n_bins=n_bins, min.size=3)

  for( i in 1:nrow(proc.regions)) {
    start = proc.regions[i, "start"]; end = proc.regions[i, "end"]
    local.ext[start:end] = Detect.Local.Extreme(x=mean.cf[start:end])
  }

  if(statFilter) {
    scale.matrix.data = matrix.data
    for( i in 1:(2*window.size) ) {
      scale.matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] =
        scale( matrix.data[ seq(1+(n_bins*i), n_bins*n_bins, 1+n_bins) ] )
    }
    for( i in 1:nrow(proc.regions)) {
      start = proc.regions[i, "start"]; end = proc.regions[i, "end"]
      pvalue[start:end] <- Get.Pvalue(matrix.data=scale.matrix.data[start:end, start:end],
                                      size=window.size, scale=1)
    }
    local.ext[intersect( union(which(local.ext==-1), which(local.ext==-1)),
                         which(pvalue<0.05))] = -2
    local.ext[which(local.ext==-1)] = 0
    local.ext[which(local.ext==-2)] = -1
  } else pvalue = 0

  domains = Convert.Bin.To.Domain.TMP(bins=bins,
                                  signal.idx=which(local.ext==-1),
                                  gap.idx=which(local.ext==-0.5),
                                  pvalues=pvalue, pvalue.cut=0.05)

  bins = cbind(bins, local.ext=local.ext, mean.cf=mean.cf, pvalue=pvalue)
  bedform = domains[, c("chr", "from.coord", "to.coord", "tag")]
  colnames(bedform) = c("chrom", "chromStart", "chromEnd", "name")

  if( !is.null(outFile) ) {
    write.table(bins,    file=paste0(outFile, ".binSignal"), quote=F, row.names=F, col.names=T, sep="\t")
    write.table(domains, file=paste0(outFile, ".domain"),    quote=F, row.names=F, col.names=T, sep="\t")
    write.table(bedform, file=paste0(outFile, ".bed"),       quote=F, row.names=F, col.names=F, sep="\t")
  }

  return(list(binSignal=bins, domain=domains, bed=bedform))
}

Get.Diamond.Matrix <- function(mat.data, i, size) {
  n_bins = nrow(mat.data); if(i==n_bins) return(NA)
  lb = max(1, i-size+1); ub = min(i+size, n_bins)
  return(mat.data[lb:i, (i+1):ub])
}

Which.process.region <- function(rmv.idx, n_bins, min.size=3) {
  gap.idx = rmv.idx
  proc.regions = data.frame(start=numeric(0), end=numeric(0))
  proc.set = setdiff(1:n_bins, gap.idx); n_proc.set = length(proc.set)
  i = 1
  while(i < n_proc.set) {
    start = proc.set[i]; j = i + 1
    while(j <= n_proc.set) {
      if(proc.set[j] - proc.set[j-1] <= 1) j = j + 1
      else { proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1])); i = j; break }
    }
    if(j >= n_proc.set) { proc.regions = rbind(proc.regions, c(start=start, end=proc.set[j-1])); break }
  }
  colnames(proc.regions) = c("start", "end")
  proc.regions[ which(abs(proc.regions[,"end"] - proc.regions[, "start"]) >= min.size), ]
}

Which.Gap.Region2 <- function(matrix.data, w) {
  n_bins = nrow(matrix.data); gap = rep(0, n_bins)
  for(i in 1:n_bins)
    if(sum(matrix.data[i, max(1, i-w):min(i+w, n_bins)]) == 0) gap[i] = -0.5
  which(gap == -0.5)
}

Detect.Local.Extreme <- function(x) {
  n_bins = length(x); ret = rep(0, n_bins); x[is.na(x)] = 0
  if(n_bins <= 3) { ret[which.min(x)] = -1; ret[which.max(x)] = 1; return(ret) }
  new.point = Data.Norm(x=1:n_bins, y=x); x = new.point$y
  cp = Change.Point(x=1:n_bins, y=x)
  if(length(cp$cp) <= 2) return(ret); if(length(cp$cp) == n_bins) return(ret)
  for(i in 2:(length(cp$cp)-1)) {
    if(x[cp$cp[i]] >= x[cp$cp[i]-1] && x[cp$cp[i]] >= x[cp$cp[i]+1]) ret[cp$cp[i]] = 1
    else if(x[cp$cp[i]] < x[cp$cp[i]-1] && x[cp$cp[i]] < x[cp$cp[i]+1]) ret[cp$cp[i]] = -1
    min.val = min(x[cp$cp[i-1]], x[cp$cp[i]]); max.val = max(x[cp$cp[i-1]], x[cp$cp[i]])
    if(min(x[cp$cp[i-1]:cp$cp[i]]) < min.val) ret[cp$cp[i-1] - 1 + which.min(x[cp$cp[i-1]:cp$cp[i]])] = -1
    if(max(x[cp$cp[i-1]:cp$cp[i]]) > max.val) ret[cp$cp[i-1] - 1 + which.max(x[cp$cp[i-1]:cp$cp[i]])] = 1
  }
  ret
}

Data.Norm <- function(x, y) {
  ret.x = rep(0, length(x)); ret.y = rep(0, length(y))
  ret.x[1] = x[1]; ret.y[1] = y[1]
  diff.x = diff(x); diff.y = diff(y)
  scale.x = 1 / mean(abs(diff(x))); scale.y = 1 / mean(abs(diff(y)))
  for(i in 2:length(x)) {
    ret.x[i] = ret.x[i-1] + (diff.x[i-1]*scale.x)
    ret.y[i] = ret.y[i-1] + (diff.y[i-1]*scale.y)
  }
  list(x=ret.x, y=ret.y)
}

Change.Point <- function(x, y) {
  if(length(x) != length(y)) { print("ERROR : length(x) != length(y)"); return(0) }
  n_bins = length(x); Fv = rep(NA, n_bins); Ev = rep(NA, n_bins); cp = 1
  i = 1; Fv[1] = 0
  while(i < n_bins) {
    j = i + 1; Fv[j] = sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2)
    while(j < n_bins) {
      j = j + 1; k = (i+1):(j-1)
      Ev[j] = (sum(abs((y[j]-y[i])*x[k] - (x[j]-x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]))) / sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2))
      Fv[j] = sqrt((x[j]-x[i])^2 + (y[j]-y[i])^2) - Ev[j]
      if(is.na(Fv[j]) || is.na(Fv[j-1])) { j = j - 1; cp = c(cp, j); break }
      if(Fv[j] < Fv[j-1]) { j = j - 1; cp = c(cp, j); break }
    }
    i = j
  }
  cp = c(cp, n_bins); list(cp=cp, objF=Fv, errF=Ev)
}

Get.Pvalue <- function(matrix.data, size, scale=1) {
  n_bins = nrow(matrix.data); pvalue = rep(1, n_bins)
  for(i in 1:(n_bins-1)) {
    dia = as.vector(Get.Diamond.Matrix2(matrix.data, i, size=size))
    ups = as.vector(Get.Upstream.Triangle(matrix.data, i, size=size))
    downs = as.vector(Get.Downstream.Triangle(matrix.data, i, size=size))
    wil.test = wilcox.test(x=dia*scale, y=c(ups, downs), alternative="less", exact=F)
    pvalue[i] = wil.test$p.value
  }
  pvalue[is.na(pvalue)] = 1; pvalue
}

Get.Upstream.Triangle <- function(mat.data, i, size) {
  lower = max(1, i-size); tmp.mat = mat.data[lower:i, lower:i]
  tmp.mat[upper.tri(tmp.mat, diag=F)]
}

Get.Downstream.Triangle <- function(mat.data, i, size) {
  n_bins = nrow(mat.data); if(i == n_bins) return(NA)
  ub = min(i+size, n_bins); tmp.mat = mat.data[(i+1):ub, (i+1):ub]
  tmp.mat[upper.tri(tmp.mat, diag=F)]
}

Get.Diamond.Matrix2 <- function(mat.data, i, size) {
  n_bins = nrow(mat.data); new.mat = matrix(rep(NA, size*size), nrow=size, ncol=size)
  for(k in 1:size) {
    if(i-(k-1) >= 1 && i < n_bins) {
      lower = min(i+1, n_bins); upper = min(i+size, n_bins)
      new.mat[size-(k-1), 1:(upper-lower+1)] = mat.data[i-(k-1), lower:upper]
    }
  }
  new.mat
}

Convert.Bin.To.Domain.TMP <- function(bins, signal.idx, gap.idx, pvalues=NULL, pvalue.cut=NULL) {
  n_bins = nrow(bins)
  ret = data.frame(chr=character(0), from.id=numeric(0), from.coord=numeric(0),
                   to.id=numeric(0), to.coord=numeric(0), tag=character(0), size=numeric(0))
  levels(x=ret[, "tag"]) = c("domain", "gap", "boundary")

  rmv.idx = setdiff(1:n_bins, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  n_procs = nrow(proc.region)
  gap = data.frame(chr=rep(bins[1,"chr"], n_procs), from.id=rep(0,n_procs),
                   from.coord=from.coord, to.id=rep(0,n_procs), to.coord=rep(0,n_procs),
                   tag=rep("gap",n_procs), size=rep(0,n_procs), stringsAsFactors=F)

  rmv.idx = union(signal.idx, gap.idx)
  proc.region = Which.process.region(rmv.idx, n_bins, min.size=0)
  n_procs = nrow(proc.region)
  from.coord = bins[proc.region[, "start"], "from.coord"]
  domain = data.frame(chr=rep(bins[1,"chr"], n_procs), from.id=rep(0,n_procs),
                      from.coord=from.coord, to.id=rep(0,n_procs), to.coord=rep(0,n_procs),
                      tag=rep("domain",n_procs), size=rep(0,n_procs), stringsAsFactors=F)

  rmv.idx = setdiff(1:n_bins, signal.idx)
  proc.region = as.data.frame(Which.process.region(rmv.idx, n_bins, min.size=1))
  n_procs = nrow(proc.region)
  if(n_procs > 0) {
    from.coord = bins[proc.region[, "start"]+1, "from.coord"]
    boundary = data.frame(chr=rep(bins[1,"chr"],n_procs), from.id=rep(0,n_procs),
                          from.coord=from.coord, to.id=rep(0,n_procs), to.coord=rep(0,n_procs),
                          tag=rep("boundary",n_procs), size=rep(0,n_procs), stringsAsFactors=F)
    ret = rbind(ret, boundary)
  }

  ret = rbind(gap, domain); ret = ret[order(ret[,3]), ]
  ret[, "to.coord"] = c(ret[2:nrow(ret), "from.coord"], bins[n_bins, "to.coord"])
  ret[, "from.id"] = match(ret[, "from.coord"], bins[, "from.coord"])
  ret[, "to.id"]   = match(ret[, "to.coord"],   bins[, "to.coord"])
  ret[, "size"]    = ret[,"to.coord"] - ret[,"from.coord"]

  if(!is.null(pvalues) && !is.null(pvalue.cut)) {
    for(i in 1:nrow(ret)) {
      if(ret[i, "tag"] == "domain") {
        idx = ret[i, "from.id"]:ret[i, "to.id"]
        if(length(idx) == length(which(pvalues[idx] < pvalue.cut))) ret[i,"tag"] = "boundary"
      }
    }
  }
  ret
}

# ---------------- driver: read samples.tsv & iterate -------------------------

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(samples_tsv  = NULL, matrix_dir = NULL, results_dir = NULL,
               window_size  = 10, samples = NULL, chroms = NULL)
  i <- 1
  while(i <= length(args)) {
    a <- args[i]
    if(a == "--samples-tsv")      { opts$samples_tsv <- args[i+1]; i <- i + 2 }
    else if(a == "--matrix-dir")  { opts$matrix_dir <- args[i+1];  i <- i + 2 }
    else if(a == "--results-dir") { opts$results_dir <- args[i+1]; i <- i + 2 }
    else if(a == "--window-size") { opts$window_size <- as.numeric(args[i+1]); i <- i + 2 }
    else if(a == "--samples")     {
      j <- i + 1
      while(j <= length(args) && !startsWith(args[j], "--")) j <- j + 1
      opts$samples <- args[(i+1):(j-1)]; i <- j
    } else if(a == "--chroms")    {
      j <- i + 1
      while(j <= length(args) && !startsWith(args[j], "--")) j <- j + 1
      opts$chroms <- args[(i+1):(j-1)]; i <- j
    } else {
      cat("[WARN] unknown arg:", a, "\n"); i <- i + 1
    }
  }
  opts
}

opts <- parse_args()

# Resolve paths relative to this script's directory.
script_dir <- tryCatch({
  ca <- commandArgs(trailingOnly = FALSE)
  fa <- ca[grep("^--file=", ca)]
  if(length(fa) == 1) dirname(normalizePath(sub("^--file=", "", fa)))
  else normalizePath(".")
}, error = function(e) normalizePath("."))

if(is.null(opts$samples_tsv)) opts$samples_tsv <- file.path(script_dir, "samples.tsv")
if(is.null(opts$matrix_dir))  opts$matrix_dir  <- file.path(script_dir, "temp_matrix")
if(is.null(opts$results_dir)) opts$results_dir <- file.path(script_dir, "results")
if(is.null(opts$chroms)) {
  opts$chroms <- c(paste0("chr", 1:22), "chrX")
}
dir.create(opts$results_dir, showWarnings = FALSE, recursive = TRUE)

samples_df <- read.table(opts$samples_tsv, sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE, comment.char = "#",
                         col.names = c("name", "hic"))
cell_names <- samples_df$name
if(!is.null(opts$samples)) cell_names <- intersect(cell_names, opts$samples)

cat("Samples:", paste(cell_names, collapse=", "), "\n")
cat("Chroms :", paste(opts$chroms, collapse=", "), "\n")
cat("Window :", opts$window_size, "\n")

for(cell_name in cell_names) {
  for(chr in opts$chroms) {
    matrix_file <- file.path(opts$matrix_dir, paste0(cell_name, "_", chr, ".matrix"))
    output_file <- file.path(opts$results_dir, paste0(cell_name, "_", chr))
    if(!file.exists(matrix_file) || file.info(matrix_file)$size == 0) {
      cat("[skip] missing or empty:", matrix_file, "\n"); next
    }
    if(file.exists(paste0(output_file, ".bed"))) {
      cat("[skip] already done:", output_file, ".bed\n"); next
    }
    cat("=== TopDom:", cell_name, chr, "===\n")
    res <- try(TopDom(matrix_file, window.size = opts$window_size,
                      outFile = output_file, statFilter = TRUE), silent = FALSE)
    if(inherits(res, "try-error")) {
      cat("[ERROR]", cell_name, chr, ":", conditionMessage(attr(res, "condition")), "\n")
    }
  }
}

cat("Done.\n")
