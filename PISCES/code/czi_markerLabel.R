library(biomaRt)
library(viper)
library(ggplot2)
library(optparse)
library(umap)

## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input gene expression matrix (proteinsXsamples).'),
  make_option(c('-m', '--marker_table'), type="character", help='Marker table file.'),
  make_option(c('-u', '--umap'), type="character", help='UMAP rds.'),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))

## gene name conversion function
Ensemble2GeneName<-function(dataset2Convert) {
  require(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  names_dataset<-getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id'),filters = 'ensembl_gene_id', values = (rownames(dataset2Convert)), mart = ensembl)
  rownames(names_dataset)<-make.unique(names_dataset$ensembl_gene_id)
  dataset2Convert_2<-merge(names_dataset,dataset2Convert,by=c("row.names"))
  dim(dataset2Convert_2)
  rownames(dataset2Convert_2)<-make.unique(dataset2Convert_2$hgnc_symbol)
  GeneName_dataset<-dataset2Convert_2[,-c(1:4)]
  head(GeneName_dataset[,1:4])
  return(GeneName_dataset)
}
## rank transformation function
RankTransform <- function(mat) {
  # replace values with their ranks in the cell (column)
  rank.mat <- apply(mat, 2, rank)
  # get the median rank of each gene (row)
  median <- apply(rank.mat, 1, median)
  # get the mad of each gene (row)
  mad <- apply(rank.mat, 1, mad)
  rank.mat <- (rank.mat - median) / mad
  return(rank.mat)
}
## dummy regulon function
DummyReg <- function(dat.mat, marker.table) {
  cellTypes <- names(table(marker.table[,1]))
  dummy.reg <- list() 
  for (i in 1:length(cellTypes)) {
    subTable <- marker.table[which(marker.table[,1] == cellTypes[i]),]
    markers <- intersect(subTable[,2], rownames(dat.mat))
    if (length(markers) > 0) {
      likelihood <- subTable[which(subTable[,2] %in% markers), 3]
      tf.mode <- likelihood; names(tf.mode) <- markers
      reg <- list('tfmode' = tf.mode, 'likelihood' = abs(likelihood))
      dummy.reg[[cellTypes[i]]] <- reg
    }
  }
  class(dummy.reg) <- "regulon"
  return(dummy.reg)
}

## read in marker and starting file
gExp <- readRDS(opt$input_file)
gExp <- Ensemble2GeneName(gExp)
marker.table <- read.csv(opt$marker_table, header = TRUE)
dat.umap <- readRDS(opt$umap)
plot.dat <- data.frame('UMAP1' = dat.umap[,1], 'UMAP2' = dat.umap[,2])

### ROUND 1 ###
gExp.r1 <- RankTransform(gExp)
## viper with dummy regulon
r1.reg <- DummyReg(gExp.r1, marker.table)
markerVip.r1 <- viper(eset = gExp.r1, regulon = r1.reg, method = 'none', minsize = 1)
## plots based on the viper results
plot.pref <- paste(opt$out_dir, opt$out_name, '_allCells-markerVip_', sep = '')
p1.dat <- cbind(plot.dat, as.data.frame(t(markerVip.r1)))
colnames(p1.dat) <- gsub(' ', '.', colnames(p1.dat))
colnames(p1.dat) <- gsub('-', '.', colnames(p1.dat))
cellTypes <- colnames(p1.dat)[3:12]
for (i in 1:length(cellTypes)) {
  ggplot(p1.dat, aes(x=UMAP1, y=UMAP2)) + geom_point(aes(colour=get(cellTypes[i]))) + 
    scale_color_gradient2(low = 'blue', mid = "grey", high = 'red', name = NULL, midpoint = 0, 
                          limits = c(-4.1, 4.1), labels = waiver()) + ggtitle(cellTypes[i]) + 
    ggsave(paste(plot.pref, cellTypes[i], '.jpeg', sep = ''), width = 8, height = 8)
}
## find the non t-cells
NES.thresh <- 2
bCells <- which(markerVip.r1['B Cells',] > NES.thresh)
fibs <- which(markerVip.r1['Fibroblasts',] > NES.thresh)
macrophage <- which(markerVip.r1['Macrophage',] > NES.thresh)
mastCells <- which(markerVip.r1['Mast Cells',] > NES.thresh)
neutro <- which(markerVip.r1['Neutrophils',] > NES.thresh)
plasmaCells <- which(markerVip.r1['Plasma Cells',] > NES.thresh)
non.tCell <- unique(c(bCells, fibs, macrophage, mastCells, neutro, plasmaCells))
tCells <- colnames(gExp.r1)[setdiff(1:ncol(gExp.r1), non.tCell)]
tCell.gExp <- gExp[, tCells]
saveRDS(tCell.gExp, file = paste(opt$out_dir, opt$out_name, '_tCells-cpm.rds', sep = ''))

### ROUND 2 $###
gExp.r2 <- RankTransform(tCell.gExp)
## round 2 viper
r2.reg <- DummyReg(gExp.r2, marker.table)
markerVip.r2 <- viper(eset = gExp.r2, regulon = r2.reg, method = 'none', minsize = 1)
## round 2 plots
plot.pref <- paste(opt$out_dir, opt$out_name, '_tCells-markerVip_', sep = '')
p2.dat <- cbind(plot.dat[tCells,], as.data.frame(t(markerVip.r2)))
colnames(p2.dat) <- gsub(' ', '.', colnames(p2.dat))
colnames(p2.dat) <- gsub('-', '.', colnames(p2.dat))
cellTypes <- colnames(p2.dat)[3:12]
r <- range(markerVip.r2[c('T Cell Activation', 'CD4_T-Cells', 'CD8_T-Cells'),])
for (i in 1:length(cellTypes)) {
  ggplot(p2.dat, aes(x=UMAP1, y=UMAP2)) + geom_point(aes(colour=get(cellTypes[i]))) + 
    scale_color_gradient2(low = 'blue', mid = "grey", high = 'red', name = NULL, midpoint = 0, 
                          limits = r, labels = waiver()) + ggtitle(cellTypes[i]) + 
    ggsave(paste(plot.pref, cellTypes[i], '.jpeg', sep = ''), width = 8, height = 8)
}
## separating activated T Cells
activated <- colnames(markerVip.r2)[which(markerVip.r2['T Cell Activation',] > 0)]
cell.Types <- data.frame('Cell Type' = rep('Resting T Cell', ncol(gExp)), stringsAsFactors = FALSE); rownames(cell.Types) <- colnames(gExp)
cell.Types[activated,] <- 'Activated T Cell'
cell.Types[non.tCell,] <- 'Non-T Cell'
write.csv(cell.Types, paste(opt$out_dir, opt$out_name, '_cellTypes.csv', sep = ''), quote = FALSE)
saveRDS(gExp[, activated], file = paste(opt$out_dir, opt$out_name, '_activated-tCells-cpm.rds', sep = ''))
