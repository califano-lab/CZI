## Author = Aaron T. Griffin
## Date = 11/27/18
## Group = Califano Lab

## Objective = compute nullmodel permutations for metaviper calculations.

# Read In Command Line Arguments using optparse
suppressMessages(library(optparse))
option_list = list(make_option(c("--test_samples"), type = "character", default = NA, help = "test RNA-Seq rawcounts matrix",metavar = "/path/to/test_samples.rda"), make_option(c("--ref_samples"), type = "character", default = NA, help = "reference RNA-Seq rawcounts matrix", metavar = "/path/to/ref_samples.rda"), make_option(c("--regulon_list"), type = "character", default = NA, help = "regulon list object", metavar = "/path/to/regulon_object.rda"), make_option(c("--per_num"), type = "numeric", default = 1000, help = "viperSignature() permutation number", metavar = "1000"), make_option(c("--norm_method"), type = "character", default = NA, help = "normalization method",metavar = 'c("none","cpm","log2cpm","rank")'), make_option(c("--core_num"), type = "numeric", default = 1, help = "core number for parallel computing", metavar = "1"), make_option(c("--output_dir"), type = "character", default = NA, help = "output directory path", metavar = "/path/to/output_dir"),make_option(c("--sub_num"), type = "integer", default = NA, help = "nullmodel subprocess number", metavar = "1"))

opt_parser = OptionParser(option_list=option_list)
param = parse_args(opt_parser)

# temp local options
running.local <- FALSE
if(running.local){
	param$test_samples <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_luad_tumor_rawmat.rds"
	param$ref_samples <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_lusc_tumor_rawmat.rda"
	param$regulon_list <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/regulon_list_1.rda"
	param$per_num <- 100
	param$norm_method <- "cpm"
	param$core_num <- 1
	param$output_dir <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/temp_scripts"
	param$sub_num <- 10
}

if(anyNA(param)){
	print_help(opt_parser)
	stop("WARNING: Please supply all parameters", call.=FALSE)
}

# Load the required packages
suppressMessages(library(viper))

# load in test samples
if(length(grep(".rds",param$test_samples)) == 1){
	test_rawmat <- readRDS(file = param$test_samples)
} else if(length(grep(".rda",param$test_samples)) == 1){
	assign("test_rawmat",value = get(load(param$test_samples)))
} else if(TRUE){
	stop("WARNING: test sample matrix must be an RDA or RDS file.")
}

# load in reference samples
if(length(grep(".rds",param$ref_samples)) == 1){
	ref_rawmat <- readRDS(file = param$ref_samples)
} else if(length(grep(".rda",param$ref_samples)) == 1){
	assign("ref_rawmat",value = get(load(param$ref_samples)))
} else if(TRUE){
	stop("WARNING: ref sample matrix must be an RDA or RDS file.")
}

# check that the rownames of the samples are identical and stop the program if they are not
rowname_check_var <- identical(rownames(test_rawmat),rownames(ref_rawmat))
if(!rowname_check_var){
	stop("WARNING: test sample matrix and reference sample matrix must have identical rownames.n", call.=FALSE)
}

# perform normalization on the test and reference matrices
# options = c("none","cpm","log2cpm","rank")

# none normalization function
none_norm_fn <- function(input.matrix){
	input.matrix <- as.matrix(input.matrix)
	input.matrix[!is.finite(input.matrix)] <- 0
	output.matrix <- input.matrix
	return(output.matrix)
}

# cpm normalization function
cpm_norm_fn <- function(input.matrix){
	input.matrix <- as.matrix(input.matrix)
	input.matrix[!is.finite(input.matrix)] <- 0
	output.matrix <- input.matrix
	for(i in 1:ncol(output.matrix)){
		output.matrix[,i] <- 1E6*output.matrix[,i]/sum(output.matrix[,i])
	}
	return(output.matrix)
}

# log2cpm normalization function
log2cpm_norm_fn <- function(input.matrix){
	input.matrix <- as.matrix(input.matrix)
	input.matrix[!is.finite(input.matrix)] <- 0
	output.matrix <- input.matrix
	for(i in 1:ncol(output.matrix)){
		output.matrix[,i] <- log((1E6*output.matrix[,i]/sum(output.matrix[,i]) + 1), base = 2)
	}
	return(output.matrix)
}

# rank normalization function
rank_norm_fn <- function(input.matrix){
	input.matrix <- as.matrix(input.matrix)
	input.matrix[!is.finite(input.matrix)] <- 0
	output.matrix <- input.matrix
	for(i in 1:ncol(output.matrix)){
		output.matrix[,i] <- rank(output.matrix[,i], ties.method = "average")
	}
	return(output.matrix)
}

if(param$norm_method == "none"){
	test_norm <- none_norm_fn(test_rawmat)
	ref_norm <- none_norm_fn(ref_rawmat)
} else if(param$norm_method == "cpm"){
	test_norm <- cpm_norm_fn(test_rawmat)
	ref_norm <- cpm_norm_fn(ref_rawmat)
} else if(param$norm_method == "log2cpm"){
	test_norm <- log2cpm_norm_fn(test_rawmat)
	ref_norm <- log2cpm_norm_fn(ref_rawmat)
} else if(param$norm_method == "rank"){
	test_norm <- rank_norm_fn(test_rawmat)
	ref_norm <- rank_norm_fn(ref_rawmat)
} else if(TRUE){
	stop("WARNING: normalization method must be specified")
}

if(running.local){
print(summary(colSums(test_norm)))
print(summary(colSums(ref_norm)))
}

# compute 10 nullmodel permutations using ttestNull

sub_nullmodel <- ttestNull(x = test_norm, y = ref_norm, per = 10, seed = param$sub_num, cores = param$core_num, verbose = FALSE)
sub_nullmodel[is.na(sub_nullmodel)] <- 0

# save the sub_nullmodel object
sub_nullmodel_savename <- paste(param$output_dir,paste("/sub_null_",param$sub_num,".rda",sep = ""),sep = "")

save(sub_nullmodel,file = sub_nullmodel_savename)

