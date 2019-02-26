Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP001swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP002swap.filtered.matrix.txt.bz2 --out_name=d1-lung --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP003swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP004swap.filtered.matrix.txt.bz2 --out_name=d1-boneMarrow --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP005swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP006swap.filtered.matrix.txt.bz2 --out_name=d1-lymphNode --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP009swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP010swap.filtered.matrix.txt.bz2 --out_name=d2-lung --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP011swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP012swap.filtered.matrix.txt.bz2 --out_name=d2-boneMarrow --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP013swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP014swap.filtered.matrix.txt.bz2 --out_name=d2-lymphNode --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP015swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP016swap.filtered.matrix.txt.bz2 --out_name=d2-spleen --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP017swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP018swap.filtered.matrix.txt.bz2 --out_name=d3-blood --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/preProcess.R --rest_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP019swap.filtered.matrix.txt.bz2 --act_file=/ifs/scratch/c2b2/ac_lab/CZI/peter_data/PP020swap.filtered.matrix.txt.bz2 --out_name=d4-blood --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d1-lung_mergedRaw.rds --out_name=d1-lung --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d1-boneMarrow_mergedRaw.rds --out_name=d1-boneMarrow --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d1-lymphNode_mergedRaw.rds --out_name=d1-lymphNode --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d2-lung_mergedRaw.rds --out_name=d2-lung --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d2-boneMarrow_mergedRaw.rds --out_name=d2-boneMarrow --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d2-lymphNode_mergedRaw.rds --out_name=d2-lymphNode --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d2-spleen_mergedRaw.rds --out_name=d2-spleen --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d3-blood_mergedRaw.rds --out_name=d3-blood --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl
Rscript /ifs/scratch/c2b2/ac_lab/CZI/lv_pipeline/qualityControl.R --input_file=/ifs/scratch/c2b2/ac_lab/CZI/preProcessed/d4-blood_mergedRaw.rds --out_name=d4-blood --out_dir=/ifs/scratch/c2b2/ac_lab/CZI/qualityControl