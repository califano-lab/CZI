library(optparse)
## load arguments
option_list <- list(
  make_option(c('-i', '--input_file'), type="character", help='Input .rds.'),
  make_option(c('-o', '--out_file'), type="character", help='Output .tsv.')
)
opt <- parse_args(OptionParser(option_list = option_list))
### read in data, then convert
rds_file <- readRDS(opt$input_file)
rds_file <- rds_file[,sample(colnames(rds_file), 1000)]
ARACNeTable <- function(data, file) {
  sample.names <- colnames(data)
  gene.ids <- rownames(data)
  m <- data
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = file , 
               sep="\t", quote = F , row.names = F , col.names = F )
}
ARACNeTable(rds_file, opt$out_file)
