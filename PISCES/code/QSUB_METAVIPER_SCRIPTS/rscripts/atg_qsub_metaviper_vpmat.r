## Author = Aaron T. Griffin
## Date = 11/27/18
## Group = Califano Lab

## Objective = compute vpmat protein activity using previously computed eset, nullmodel, and elements of the regulon list.

# Read In Command Line Arguments using optparse
suppressMessages(library(optparse))
option_list = list(make_option(c("--test_samples"), type = "character", default = NA, help = "test RNA-Seq rawcounts matrix",metavar = "/path/to/test_samples.rda"), make_option(c("--ref_samples"), type = "character", default = NA, help = "reference RNA-Seq rawcounts matrix", metavar = "/path/to/ref_samples.rda"), make_option(c("--regulon_list"), type = "character", default = NA, help = "regulon list object", metavar = "/path/to/regulon_object.rda"), make_option(c("--per_num"), type = "numeric", default = 1000, help = "viperSignature() permutation number", metavar = "1000"), make_option(c("--core_num"), type = "numeric", default = 1, help = "core number for parallel computing", metavar = "1"), make_option(c("--output_dir"), type = "character", default = NA, help = "output directory path", metavar = "/path/to/output_dir"),make_option(c("--sub_num"), type = "integer", default = NA, help = "metaViper subprocess number", metavar = "1"))

opt_parser = OptionParser(option_list=option_list)
param = parse_args(opt_parser)

# temp local options
running.local <- FALSE
if(running.local){
	param$test_samples <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_luad_tumor_rawmat.rds"
	param$ref_samples <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_lusc_tumor_rawmat.rda"
	param$regulon_list <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/regulon_list_1.rds"
	param$per_num <- 100
	param$norm_method <- "cpm"
	param$core_num <- 1
	param$output_dir <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/temp_scripts"
	param$sub_num <- 6
}


if(anyNA(param)){
	print_help(opt_parser)
	stop("WARNING: Please supply all parameters", call.=FALSE)
}

# Load the required packages
suppressMessages(library(viper))

# Load the current_signature object
current_signature <- readRDS(paste(param$output_dir,"/current_signature.rds",sep = ""))

# load the regulon list object
if(substr(param$regulon_list, start = nchar(param$regulon_list)-3, stop = nchar(param$regulon_list)) == ".rds"){
	regulon_list <- readRDS(file = param$regulon_list)
} else if(substr(param$regulon_list, start = nchar(param$regulon_list)-3, stop = nchar(param$regulon_list)) == ".rda"){
	assign("regulon_list",value = get(load(param$regulon_list)))
} else if(TRUE){
	stop("WARNING: regulon list object must be an RDA or RDS file.")
}

if(is(regulon_list,"regulon")){
	temp <- list(first_reg = regulon_list)
	regulon_list <- temp
}

# set the current regulon
current_regulon <- regulon_list[[param$sub_num]]

# calculate protein activity vpmat
if(param$per_num == 0){
	current_signature <- as.matrix(current_signature)
	vpmat <- viper(eset = current_signature, regulon = current_regulon, method = "none", cores = param$core_num, verbose = FALSE)
} else if(param$per_num != 0){
	vpmat <- viper(current_signature, regulon = current_regulon, cores = param$core_num, verbose = FALSE)
}

# save the protein activity vpmat
vpmat_savename <- paste(param$output_dir,paste("/sub_vpmat_",param$sub_num,sep = ""),".rda",sep = "")
save(vpmat,file = vpmat_savename)
