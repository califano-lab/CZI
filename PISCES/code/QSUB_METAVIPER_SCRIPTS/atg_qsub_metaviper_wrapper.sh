#!/bin/bash
## Author = Aaron T. Griffin
## Date = 11/27/18
## Group = Califano Lab

## Objective = Perform parallelized metaVIPER using the isilon cluster's qsub capabilities.

# set internal parameters for core usage, RAM usage per qsub job, time limit per qsub job, and time in between verbose messages
CORE_NUM=1
MEMORY_LIMIT=16G
TIME_LIMIT=3::
SLEEP_TIME=60

# define the path to the rscripts which will be called by the bash wrapper for qsub metaVIPER
qsub_metaviper_rscripts_directory="/ifs/scratch/c2b2/ac_lab/atg2142/Fall_2018/QSUB_METAVIPER_SCRIPTS/rscripts"

# read in parameters from the command line
for i in "$@"
do
case $i in
    --test_samples=*)
    TEST_SAMPLES="${i#*=}"
    shift # past argument=value
    ;;
    --ref_samples=*)
    REF_SAMPLES="${i#*=}"
    shift # past argument=value
    ;;
    --regulon_list=*)
    REGULON_LIST="${i#*=}"
    shift # past argument=value
    ;;
    --per_num=*)
    PER_NUM="${i#*=}"
    shift # past argument=value
    ;;
    --norm_method=*)
    NORM_METHOD="${i#*=}"
    shift # past argument=value
    ;;
    --combine_method=*)
    COMBINE_METHOD="${i#*=}"
    shift # past argument=value
    ;;
    --output_dir=*)
    OUTPUT_DIR="${i#*=}"
    shift # past argument=value
    ;;
    --keep=*)
    KEEP="${i#*=}"
    shift # past argument=value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument with no value
    ;;
    *)
          # unknown option
    ;;
esac
done

# echo input parameters to the terminal

echo ""
echo ""
echo "*** CLUSTER METAVIPER PARAMETERS ***"
echo "- TEST SAMPLES = ${TEST_SAMPLES} "
echo "- REFERENCE SAMPLES = ${REF_SAMPLES} " 
echo "- REGULON LIST = ${REGULON_LIST} "
echo "- NORMALIZATION METHOD = ${NORM_METHOD} "
echo "- COMBINATION METHOD = ${COMBINE_METHOD} "
echo "- PERMUTATIONS = ${PER_NUM} "
echo "- OUTPUT DIRECTORY = ${OUTPUT_DIR} "

# check if the input files exist and stop the script if they do not
echo ""
echo "*** CHECKING FOR INPUT FILES ***"
continue_var="TRUE"

if [ -e ${TEST_SAMPLES} ]
then
    echo "- test samples exist "
else
    echo "- test samples do not exist "
    continue_var="FALSE"
fi

if [ -e ${REF_SAMPLES} ]
then
    echo "- reference samples exist "
else
    echo "- reference samples do not exist "
    continue_var="FALSE"
fi

if [ -e ${REGULON_LIST} ]
then
    echo "- regulon list exists "
else
    echo "- regulon list does not exist "
    continue_var="FALSE"
fi

if [ -e ${OUTPUT_DIR} ]
then
    echo "- output directory exists "
else
    echo "- output directory does not exist "
    continue_var="FALSE"
fi

if [ ${continue_var} = "FALSE" ]
then
    echo ""
    echo "*** STOPPING COMPUTATION *** "
    exit 1
fi

# if permutation number is greater than 0, compute null model permutations in parallel using individual qsubs (10 permutations per qsub)

if [ ${PER_NUM} -gt 0 ]
then

    echo ""
    echo "*** COMPUTING NULLMODEL (${CORE_NUM} CORE(S)) (${PER_NUM} permutations) ***"
    current_time=$(date)

    echo "- Nullmodel start time = ${current_time}"

    nullmodel_job_number=$((${PER_NUM}/10))

    for i in `seq 1 $nullmodel_job_number` ; do
        echo "- initializing nullmodel calculation $i/$nullmodel_job_number "
        command="Rscript --vanilla "${qsub_metaviper_rscripts_directory}"/atg_qsub_metaviper_nullmodel.r --test_samples=${TEST_SAMPLES} --ref_samples=${REF_SAMPLES} --regulon_list=${REGULON_LIST} --per_num=${PER_NUM} --core_num=${CORE_NUM} --norm_method=${NORM_METHOD} --output_dir=${OUTPUT_DIR} --sub_num=$i"

        echo '#!/bin/bash' > ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#$ -cwd' >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#$ -j y' >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#$ -o' ${OUTPUT_DIR} >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#$ -l mem='${MEMORY_LIMIT}',time='${TIME_LIMIT} >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#$ -S /bin/bash' >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#$ -N' metaviper_nullmodel_sub_${i} >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo '#' >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        echo ${command} >> ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
        qsub ${OUTPUT_DIR}/metaviper_nullmodel_sub_${i}.sh
    done

    # check whether or not the nullmodel calculations have finished 
    nullmodel_complete_var=0

    while [ ${nullmodel_complete_var} = 0 ]
    do
        echo "- nullmodel calculation ongoing ..."
        sleep ${SLEEP_TIME}
        sub_check_var=1
        for i in `seq 1 $nullmodel_job_number` ; do
            if [ -e ${OUTPUT_DIR}/sub_null_${i}.rda ] ; then
                sub_check_var=$((sub_check_var*1))
            else
                sub_check_var=$((sub_check_var*0))
            fi
        done
        nullmodel_complete_var=${sub_check_var}
    done

fi

# compute the gene expression signature 
echo ""
echo "*** COMPUTING GENE EXPRESSION SIGNATURE  (${CORE_NUM} CORES) ***"
current_time=$(date)
echo "- gene expression signature start time = ${current_time}"

command="Rscript --vanilla "${qsub_metaviper_rscripts_directory}"/atg_qsub_metaviper_signature.r --test_samples=${TEST_SAMPLES} --ref_samples=${REF_SAMPLES} --regulon_list=${REGULON_LIST} --per_num=${PER_NUM} --core_num=${CORE_NUM} --norm_method=${NORM_METHOD} --output_dir=${OUTPUT_DIR}"

echo '#!/bin/bash' > ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#$ -cwd' >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#$ -j y' >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#$ -o' ${OUTPUT_DIR} >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#$ -l mem='${MEMORY_LIMIT}',time='${TIME_LIMIT} >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#$ -S /bin/bash' >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#$ -N' metaviper_vpsig_job >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo '#' >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
echo ${command} >> ${OUTPUT_DIR}/metaviper_vpsig_script.sh
qsub ${OUTPUT_DIR}/metaviper_vpsig_script.sh

# check whether or not the gene expression signature calculation has finished 
vpsig_complete_var="FALSE"

while [ ${vpsig_complete_var} = "FALSE" ]
do
    echo "- gene expression signature calculation ongoing ..."
    sleep ${SLEEP_TIME}
    if [ -e ${OUTPUT_DIR}/regulon_list_length.txt ]
    then
        vpsig_complete_var=TRUE
    fi
done

# compute protein activity for each member of the regulon using individual qsubs
echo ""
echo "*** COMPUTING VIPER PROTEIN ACTIVITY (${CORE_NUM} CORES) ***"
current_time=$(date)

echo "- Viper start time = ${current_time}"

while IFS="" read -r line || [[ -n "$line" ]]; do
    number_of_regulons=$line
done < "${OUTPUT_DIR}/regulon_list_length.txt"

for i in `seq 1 $number_of_regulons` ; do
    echo "- initializing Viper calculation $i/$number_of_regulons "
    command="Rscript --vanilla "${qsub_metaviper_rscripts_directory}"/atg_qsub_metaviper_vpmat.r --test_samples=${TEST_SAMPLES} --ref_samples=${REF_SAMPLES} --regulon_list=${REGULON_LIST} --per_num=${PER_NUM} --core_num=${CORE_NUM} --output_dir=${OUTPUT_DIR} --sub_num=$i"

    echo '#!/bin/bash' > ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#$ -cwd' >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#$ -j y' >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#$ -o' ${OUTPUT_DIR} >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#$ -l mem='${MEMORY_LIMIT}',time='${TIME_LIMIT} >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#$ -S /bin/bash' >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#$ -N' metaviper_vpmat_sub_${i} >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo '#' >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    echo ${command} >> ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
    qsub ${OUTPUT_DIR}/metaviper_vpmat_sub_${i}.sh
done

# check whether or not the vpmat calculations have finished 
vpmat_complete_var=0

while [ ${vpmat_complete_var} = 0 ]
do
    echo "- vpmat calculation ongoing ..."
    sleep ${SLEEP_TIME}
    sub_check_var=1
    for i in `seq 1 $number_of_regulons` ; do
        if [ -e ${OUTPUT_DIR}/sub_vpmat_${i}.rda ] ; then
            sub_check_var=$((sub_check_var*1))
        else
            sub_check_var=$((sub_check_var*0))
        fi
    done
    vpmat_complete_var=${sub_check_var}
done

# combine the protein activity and save the final vpmat list
echo ""
echo "*** COMBINING PROTEIN ACTIVITY ***"
current_time=$(date)
echo "- protein combination start time = ${current_time}"

command="Rscript --vanilla "${qsub_metaviper_rscripts_directory}"/atg_qsub_metaviper_combine.r --test_samples=${TEST_SAMPLES} --ref_samples=${REF_SAMPLES} --regulon_list=${REGULON_LIST} --per_num=${PER_NUM} --core_num=${CORE_NUM} --combine_method=${COMBINE_METHOD} --output_dir=${OUTPUT_DIR} --keep=${KEEP}"

echo '#!/bin/bash' > ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#$ -cwd' >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#$ -j y' >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#$ -o' ${OUTPUT_DIR} >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#$ -l mem='${MEMORY_LIMIT}',time='${TIME_LIMIT} >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#$ -S /bin/bash' >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#$ -N' metaviper_combine_job >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo '#' >> ${OUTPUT_DIR}/metaviper_combine_script.sh
echo ${command} >> ${OUTPUT_DIR}/metaviper_combine_script.sh
qsub ${OUTPUT_DIR}/metaviper_combine_script.sh

# check every 30 seconds whether or not the vpmat combination calculation has finished
combine_complete_var="FALSE"

while [ ${combine_complete_var} = "FALSE" ]
do
    echo "- protein activity combination ongoing ..."
    sleep ${SLEEP_TIME}
    if [ -e ${OUTPUT_DIR}/final_vpmat_list.rds ]
    then
        combine_complete_var=TRUE
    fi

    if [ -e ${OUTPUT_DIR}/final_vpmat.rds ]
    then
        combine_complete_var=TRUE
    fi

done

# finished
echo ""
echo "*** METAVIPER COMPUTATION COMPLETE ***"
current_time=$(date)
echo "- metaviper end time = ${current_time}"


