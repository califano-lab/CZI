## libraries
library(optparse)
library(viper)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input raw count matrix for meta cell creation (genesXsamples).'),
  make_option(c('-a', '--activity_file'), type="character", help='Matrix of protein activity for distance calculation (proteinsXsamples).'),
  make_option(c('-c', '--num_cells'), type="integer", help='Number of neighbors to use..', default=5),
  make_option(c('-s', '--subset_size'), type="integer", help='Subsample size.', default=200),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
in.dat <- readRDS(opt$input_file)
pAct.dat <- readRDS(opt$activity_file)

in.dat <- readRDS('C://Users/lvlah/linux/ac_lab/data/czi/pipelineDev/d1-lung_mergedFiltered.rds')
pAct.dat <- readRDS('C://Users/lvlah/linux/ac_lab/data/czi/pipelineDev/d1-lung_GTExActivity.rds')

## create distance matrix
pAct.dat <- pAct.dat[, colnames(in.dat)]
dist_matrix <- as.matrix(as.dist(viperSimilarity(pAct.dat)))
## comput k nearest neighbors
KNN <- function(dist.mat, k){
  dist.mat <- as.matrix(dist.mat)
  n <- nrow(dist.mat)
  neighbor.mat <- matrix(0L, nrow = n, ncol = k)
  for (i in 1:n) {
    neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
  }
  return(neighbor.mat)
}
knn.neighbors <- KNN(dist_matrix, 10)
## create imputed matrix
imputed.dat <- matrix(0, nrow = nrow(in.dat), ncol = ncol(in.dat))
rownames(imputed.dat) <- rownames(in.dat); colnames(imputed.dat) <- colnames(in.dat)
for (i in 1:ncol(in.dat)) {
  neighbor.mat <- in.dat[,c(i, knn.neighbors[i,])]
  imputed.dat[,i] <- rowSums(neighbor.mat)
}
## subset if specified
imputed.dat <- imputed.dat[,sample(colnames(imputed.dat), size = 200)]
imputed.dat <- imputed.dat[ rowSums(imputed.dat) != 0, ]
## normalize if specified
imputed.dat <- log2(t(t(imputed.dat) / (colSums(imputed.dat) / 1e6)) + 1)
save(imputed.dat, file = paste(opt$out_dir, opt$out_name, '_metaCells.rda', sep = ''))