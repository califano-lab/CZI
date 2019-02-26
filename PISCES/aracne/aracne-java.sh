## load arguments
RUN_ID=$1
BASE_DIR=$2
REG_FILE=$3
EXP_FILE=$4
## code paths
CONSOLIDATE_CODE=/ifs/scratch/c2b2/ac_lab/CZI/aracne-1.11/aracneConsolidate.r
ARACNE_BIN=/ifs/scratch/c2b2/ac_lab/CZI/aracne-1.11/aracne-ap.jar
java=/nfs/apps/java/1.8.0_91/bin/java
P_THRESH=1E-8
## create directory
RUN_DIR=${BASE_DIR}/${RUN_ID}
mkdir -p ${RUN_DIR}
## ARACNe Threshold
qsub -N aThresh-${RUN_ID} -wd ${RUN_DIR} -l mem=18G -l s_rt=2:0:0 -l h_rt=2:0:0 -b y ${java} -Xmx5G -jar ${ARACNE_BIN} \
-e ${EXP_FILE} -o ${RUN_DIR} --tfs ${REG_FILE} --pvalue ${P_THRESH} --seed 666 --calculateThreshold
## bootstraps
for i in {1..200}
do
	qsub -N ABS-${RUN_ID} -hold_jid aThresh-${RUN_ID} -wd ${RUN_DIR} \
	-l mem=24G -l s_rt=2:0:0 -l h_rt=2:0:0 -b y ${java} -Xmx12G -jar ${ARACNE_BIN} -e ${EXP_FILE} \
	-o ${RUN_DIR} --tfs ${REG_FILE} --pvalue ${P_THRESH} --seed $i
done
## consolidate
qsub -N aCon-${RUN_ID} -hold_jid ABS-${RUN_ID} -wd ${RUN_DIR} -l mem=16G -l s_rt=16:0:0 -l h_rt=16:0:0 -b y \
/nfs/apps/R/3.4.3/bin/Rscript ${CONSOLIDATE_CODE} ${RUN_DIR} ${EXP_FILE} ${REG_FILE} bonferroni 0.01
