## Author = Aaron T. Griffin
## Date = 11/27/18
## Group = Califano Lab

## Objective = compute the gene expression signature of the test samples vs. reference samples and aggregate the null model permutations. 

# Read In Command Line Arguments using optparse
suppressMessages(library(optparse))
option_list = list(make_option(c("--test_samples"), type = "character", default = NA, help = "test RNA-Seq rawcounts matrix",metavar = "/path/to/test_samples.rda"), make_option(c("--ref_samples"), type = "character", default = NA, help = "reference RNA-Seq rawcounts matrix", metavar = "/path/to/ref_samples.rda"), make_option(c("--regulon_list"), type = "character", default = NA, help = "regulon list object", metavar = "/path/to/regulon_object.rda"), make_option(c("--per_num"), type = "numeric", default = 1000, help = "viperSignature() permutation number", metavar = "1000"), make_option(c("--norm_method"), type = "character", default = NA, help = "normalization method",metavar = 'c("none","cpm","log2cpm","rank")'), make_option(c("--core_num"), type = "numeric", default = 1, help = "core number for parallel computing", metavar = "1"), make_option(c("--output_dir"), type = "character", default = NA, help = "output directory path", metavar = "/path/to/output_dir"))

opt_parser = OptionParser(option_list=option_list)
param = parse_args(opt_parser)

# temp local options
running.local <- FALSE
if(running.local){
	param$test_samples <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_luad_tumor_rawmat.rds"
	param$ref_samples <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_lusc_tumor_rawmat.rda"
	param$regulon_list <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/regulon_list_1.rda"
	param$per_num <- 0
	param$norm_method <- "cpm"
	param$core_num <- 1
	param$output_dir <- "/Volumes/ac_lab_scratch/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/temp_scripts"
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

# compute the signature of test_norm vs. ref_norm using vipersignature if per_num > 0 or median and mad if per_num = 0
if(param$per_num == 0){
	
	current_signature <- test_norm
	for(i in 1:nrow(test_norm)){
		old_values <- test_norm[i,]
		new_values <- (old_values - median(ref_norm[i,]))/mad(ref_norm[i,])
		new_values[is.na(new_values)] <- 0
		current_signature[i,] <- new_values
	}
	
	# save the current_signature matrix
	current_signature_savename <- paste(param$output_dir,"/current_signature.rds",sep = "")
	saveRDS(current_signature,file = current_signature_savename)
	
} else if(param$per_num != 0){

	current_signature <- apply(test_norm, 2, function(x, ctrl) {
		tmp <- rowTtest(x - ctrl)
		(qnorm(tmp$p.value/2, lower.tail = FALSE)*sign(tmp$statistic))[,1]
	}, ctrl = ref_norm)
	current_signature[is.na(current_signature)] <- 0
	
	# load in previously generated nullmodel data
	temp_nullmodel <- matrix(data = 0, nrow = nrow(current_signature), ncol = 10)
	rownames(temp_nullmodel) <- rownames(current_signature)
	initial_temp_size <- ncol(temp_nullmodel)
	
	sub_nullmodel_index <- grep("sub_null_",list.files(param$output_dir))
	
	for(i in 1:length(sub_nullmodel_index)){
		load(paste(param$output_dir,"/",list.files(param$output_dir)[sub_nullmodel_index[i]],sep = ""))
		temp_nullmodel <- cbind(temp_nullmodel,sub_nullmodel)
	}
	
	temp_nullmodel <- temp_nullmodel[,-c(1:initial_temp_size)]
	colnames(temp_nullmodel) <- NULL
	
	# create a vpsig object using the current_signature and temp_nullmodel
	vpsig <- list(signature = current_signature, nullmodel = temp_nullmodel)
	class(vpsig) <- "viperSignature"
	
	# save the current_signature viperSignature object
	current_signature <- vpsig
	current_signature_savename <- paste(param$output_dir,"/current_signature.rds",sep = "")
	saveRDS(current_signature,file = current_signature_savename)
}

# load the regulon list object
if(substr(param$regulon_list, start = nchar(param$regulon_list)-3, stop = nchar(param$regulon_list)) == ".rds"){
	regulon_list <- readRDS(file = param$regulon_list)
} else if(substr(param$regulon_list, start = nchar(param$regulon_list)-3, stop = nchar(param$regulon_list)) == ".rda"){
	assign("regulon_list",value = get(load(param$regulon_list)))
} else if(TRUE){
	stop("WARNING: test sample matrix must be an RDA or RDS file.")
}

if(is(regulon_list,"regulon")){
	temp <- list(first_reg = regulon_list)
	regulon_list <- temp
}

# save the length of the regulon list in a text file in the output directory
regulon_list_length <- length(regulon_list)
write.table(regulon_list_length,file = paste(param$output_dir,"/regulon_list_length.txt",sep = ""),quote = FALSE,sep = "\t",row.names=FALSE,col.names=FALSE)

