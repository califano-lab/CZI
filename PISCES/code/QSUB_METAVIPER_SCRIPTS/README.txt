[ Project Title ]

QSUB METAVIPER

[ Motivation ]

With the development of VIPER (Alvarez et. al. Nature Genetics 2016), protein activity can be inferred from gene expression data using gene regulatory networks (interactomes) created by ARACNe-AP (Lachmann et. al. Bioinformatics 2016). VIPER runs efficiently on a small number of samples in a local R environment with one interactome. However, the transformation of VIPER to metaVIPER by integrating protein activity inferred across multiple interactomes (Hongxu et. al. Nature Communications 2018) has changed the standard way in which protein activity is inferred. Now it is common to infer protein activity using a robust null model and numerous interactomes, which becomes computationally intensive in a local R environment.

In order to speed up protein activity inference, this script was created to distribute the individual aspects of metaVIPER across multiple jobs submitted to the ARCS isilon Linux server. By massively parallelizing metaVIPER, protein inference can be accomplished an order of magnitude faster than in a local R environment, thus decreasing the time from hypothesis to results for computational biologists in the Califano Lab. 

[ Build Status ]

Version 1.2 (12/04/2018)

[ Code Style ]

The atg_qusb_metaviper_wrapper.sh script is written in BASH.

The atg_qsub_metaviper_nullmodel.r, atg_qsub_metaviper_signature.r, atg_qsub_metaviper_vpmat.r, and atg_qsub_metaviper_combine.r scripts are written in R. 

[ Features and Options]

--test_samples=/path/to/test_samples.rda
	(test_samples.rda should be an R matrix object containing rawcounts for the test samples whose protein activity will be inferred; RDA and RDS file formats accepted)

--ref_samples=/path/to/ref_samples.rda
	(ref_samples.rda should be an R matrix object containing rawcounts for the reference samples which will be used to generate the gene expression signature of the test samples; RDA and RDS file formats accepted)

--regulon_list=/path/to/regulon_list.rda
	(regulon_list.rda should be an R list object containing regulons which will be used to perform metaVIPER protein activity inference; RDA and RDS file formats accepted)

--output_dir=/path/to/output_directory
	(output_directory should be an existing empty directory in the isilon ac_lab scratch directory where results will be output)

--per_num=1000
	(integer number of nullmodel permutations, must be a multiple of 10; recommended 1000)

--norm_method=none
	(method for normalizing columns of test and reference sample matrices; should be none, rank, cpm, or log2cpm)

--combine_method=Stouffer
	(method for combining protein activity inferred from different interactomes; should be Stouffer or AbsMax)

--keep=vpmat
	(variable to determine which output file are kept; should be vpmat, vpmat_list, or all)

The atg_qsub_metaviper_nullmodel.r script distributes nullmodel computation across multiple isilon jobs (10 nullmodel permutations per job). This script uses the VIPER function ttestNull() and each job is created with a different random seed. 

The atg_qsub_metaviper_signature.r script computes a gene expression signature of test samples vs. reference samples using the VIPER function rowTtest() if per_num > 0 or calculates a gene expression signature by scaling test sample matrix rows by the median and mad of the reference sample matrix rows. The script then aggregates all nullmodel files to create a final nullmodel if per_num > 0, and saves the gene expression signature and nullmodel together as a viperSignature object (equivalent to running viperSignature() locally with method = "ttest"). 

The atg_qsub_metaviper_vpmat.r script distributes protein activity inference by VIPER across multiple isilon jobs (one job per interactome in the regulon_list.rda file). The VIPER function viper() is used to compute protein activity (method = "none" if per_num = 0).

The atg_qsub_metaviper_combine.r script combines protein activity inferred from each regulon into a list of protein activity matrices. The final protein activity matrix is created with the number of columns equal to number of test samples and number of rows equal to the union of all proteins in the individual viper protein activity matrices. Final activity for a given protein is calculated by Stouffer integration of protein activity NES values inferred for that protein in each of the viper protein matrices in which it is represented (if combine_method = Stouffer) or by selection of the NES value with the greatest magnitude (if combine_method = AbsMax).

[ Installation ]

Prior to running this script the user must install two R packages on your ARCS isilon account: "viper" and "optparse". Both can be installed from bioconductor with the follow commands:

qrsh -l mem=8G,time=1::
R
source("https://bioconductor.org/biocLite.R")
biocLite("viper")
biocLite("optparse")
q()

[ Tests ]

Test data has been included in the folder /ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data :
- luad_tumor_rawmat.rda (515 TCGA LUAD primary tumor samples)
- lusc_tumor_rawmat.rda (501 TCGA LUSC primary tumor samples)
- sub_luad_tumor_rawmat.rda (50 TCGA LUAD primary tumor samples)
- sub_lusc_tumor_rawmat.rda (50 TCGA LUSC primary tumor samples)
- regulon_list_1.rda (CNV corrected TCGA interactomes from the DarwinHealth n1platform built with an ARACNe-AP bootstrap bug for BLCA, BRCA, CESC, COAD, ESCA)
- regulon_list_2.rda (CNV corrected TCGA interactomes from the DarwinHealth n1platform built with an ARACNe-AP bootstrap bug for LIHC, LUAD, LUSC, PRAD, READ)
- n1platform_regulon_list.rda (28 CNV corrected TCGA interactomes obtained from the DarwinHealth n1platform built with an ARACNe-AP bootstrap bug pruned to 50 top targets per regulator)
- tcga_regulon_list.rda (CNV corrected TCGA interactomes constructed by Hyunjin Kim from the Aristidis Floratos Lab using high-qualty primary tumor samples without the ARACNe-AP bootstrap bug pruned to 50 top targets per regulator)

To run metaVIPER with a 1000 permutation null model, Stouffer combination, cpm normalization, and keeping only the final protein activity matrix using the test data, execute the following command on an isilon login node:

my_output_directory=/path/to/output_directory
my_test_samples=/ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_luad_tumor_rawmat.rda
my_reference_samples=/ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/sub_lusc_tumor_rawmat.rda
my_regulon_list=/ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/test_data/regulon_list_1.rda

bash /ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/atg_qusb_metaviper_wrapper.sh --test_samples=${my_test_samples} --ref_samples=${my_reference_samples} --regulon_list=${my_regulon_list} --output_dir=${my_output_directory} --norm_method=cpm --per_num=1000 --combine_method=Stouffer --keep=vpmat


[ Contribute ]

Please direct comments, bug reports, and questions to Aaron T. Griffin (atg2142@cumc.columbia.edu).

[ Credit ]

Credit to Pasquale Laise for comments and guidance.

Credit to Alessandro Vasciaveo for comments and guidance. 

Credit to Lukas Vlahos for comments, troubleshooting, bug fixes, and guidance. 

Credit to Ding Hongxu for the inspiration to massively parallelize time-intensive Califano Lab programs.

Credit to Alvarez et. al. for the development of the VIPER algorithm and implementation in R.

Credit to Hongxu et. al. for the theoretical development of metaVIPER.

Credit to members of the Califano Lab for constructing interactomes from The Cancer Genome Atlas (TCGA) and Genotype-Tissue Expression (GTEx) data sets.

