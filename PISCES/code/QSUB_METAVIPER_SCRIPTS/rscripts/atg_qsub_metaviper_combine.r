## Author = Aaron T. Griffin
## Date = 10/20/18
## Group = Califano Lab

## Objective = combine vpmat protein activity to produce final vpmat_list object.

# Read In Command Line Arguments using optparse
suppressMessages(library(optparse))
option_list = list(make_option(c("--test_samples"), type = "character", default = NA, help = "test RNA-Seq rawcounts matrix",metavar = "/path/to/test_samples.rda"), make_option(c("--ref_samples"), type = "character", default = NA, help = "reference RNA-Seq rawcounts matrix", metavar = "/path/to/ref_samples.rda"), make_option(c("--regulon_list"), type = "character", default = NA, help = "regulon list object", metavar = "/path/to/regulon_object.rda"), make_option(c("--per_num"), type = "numeric", default = 1000, help = "viperSignature() permutation number", metavar = "1000"), make_option(c("--core_num"), type = "numeric", default = 1, help = "core number for parallel computing", metavar = "1"), make_option(c("--output_dir"), type = "character", default = NA, help = "output directory path", metavar = "/path/to/output_dir"), make_option(c("--combine_method"), type = "character", default = NA, help = "metaVIPER combination method", metavar = 'c("Stouffer","AbsMax")'),make_option(c("--keep"), type = "character", default = NA, help = "specify which files to keep", metavar = 'c("vpmat","vpmat_list","all")'))

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
	param$combine_method <- "Stouffer"
	param$keep <- "all"
}


if(anyNA(param)){
	print_help(opt_parser)
	stop("WARNING: Please supply all parameters", call.=FALSE)
}

# Load the required packages
suppressMessages(library(viper))

# define Stouffer integration function for combining NES values
metaVIPER_Stouffer_combine <- function(vpmat){
	final_activity_vector <- rep(NA,times = ncol(vpmat))
	for(i in 1:ncol(vpmat)){
		temp_vec <- vpmat[,i]
		temp_vec <- temp_vec[is.finite(temp_vec)]
		temp_res <- sum(temp_vec)/sqrt(length(temp_vec))
		final_activity_vector[i] <- temp_res
	}
	return(final_activity_vector)
}

# define AbsMax integration function for combining NES values
metaVIPER_AbsMax_combine <- function(vpmat){
	final_activity_vector <- rep(NA,times = ncol(vpmat))
	for(i in 1:ncol(vpmat)){
		temp_vec <- vpmat[,i]
		temp_vec <- temp_vec[is.finite(temp_vec)]
		temp_res <- temp_vec[(abs(temp_vec) == max(abs(temp_vec)))]
		final_activity_vector[i] <- temp_res
	}
	return(final_activity_vector)
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

# create the vpmat_list object
vpmat_list <- list()

# load each vpmat and add it to the list

for(i in 1:length(regulon_list)){
	load(paste(param$output_dir,paste("/sub_vpmat_",i,sep = ""),".rda",sep = ""))
	vpmat_list[[i]] <- vpmat
}
if(!is.null(names(regulon_list)[1])){
	names(vpmat_list) <- paste(names(regulon_list),rep("_vpmat",times = length(regulon_list)),sep = "")
}

# remove the regulon_list object from the environment to free up memory
remove(regulon_list)

# compute final protein activity by taking the average protein activity from all vpmat objects in the list 
final_proteins <- NULL
for(i in 1:length(vpmat_list)){
	final_proteins <- c(final_proteins,rownames(vpmat_list[[i]]))
}
final_proteins <- unique(final_proteins)

sample_names <- colnames(vpmat_list[[1]])

final_vpmat <- matrix(data = 0, nrow = length(final_proteins), ncol = length(sample_names), dimnames = list(final_proteins,sample_names))

for(i in 1:nrow(final_vpmat)){
	query_protein <- rownames(final_vpmat)[i]
	temp_activity_mat <- matrix(data = NA, nrow = length(vpmat_list),ncol = ncol(final_vpmat))
	colnames(temp_activity_mat) <- colnames(final_vpmat)
	rownames(temp_activity_mat) <- names(vpmat_list)
	for(k in 1:length(vpmat_list)){
		sub_vpmat_index <- grep(paste("_",query_protein,"_",sep = ""),paste("_",rownames(vpmat_list[[k]]),"_",sep = ""), fixed = TRUE)
		if(length(sub_vpmat_index) == 1){
			temp_activity_mat[k,] <- vpmat_list[[k]][sub_vpmat_index,]
		} 
	}
	if(param$combine_method == "Stouffer"){
		final_vpmat[i,] <- metaVIPER_Stouffer_combine(temp_activity_mat)
	} else if(param$combine_method == "AbsMax"){
		final_vpmat[i,] <- metaVIPER_AbsMax_combine(temp_activity_mat)
	} else if(TRUE){
		stop("WARNING: normalization method must be Stouffer or AbsMax")
	}
}

if(param$keep == "vpmat"){
	final_vpmat_savename <- paste(param$output_dir,"/final_vpmat.rds",sep = "")
	saveRDS(final_vpmat,file = final_vpmat_savename)
} else if(param$keep != "vpmat"){
	vpmat_list$final_vpmat <- final_vpmat
	final_vpmat_savename <- paste(param$output_dir,"/final_vpmat_list.rds",sep = "")
	saveRDS(vpmat_list,file = final_vpmat_savename)
}

# clean up temporary files
if(param$keep != "all"){
	command <- paste("rm ",param$output_dir,"/*.o*",sep = "")
	system(command = command)
	command <- paste("rm ",param$output_dir,"/metaviper_nullmodel_sub_*",sep = "")
	system(command = command)
	command <- paste("rm ",param$output_dir,"/sub_null_*",sep = "")
	system(command = command)
	command <- paste("rm ",param$output_dir,"/metaviper_*.sh",sep = "")
	system(command = command)
	command <- paste("rm ",param$output_dir,"/sub_vpmat_*.rda",sep = "")
	system(command = command)
	command <- paste("rm ",param$output_dir,"/regulon_list_length.txt",sep = "")
	system(command = command)
	command <- paste("rm ",param$output_dir,"/sub_vpmat_*.rda",sep = "")
	system(command = command)
}






