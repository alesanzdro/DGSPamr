#!/bin/bash

# Example of how to set the environment variable in the bash shell. Remember this is only temporary, if you want it set every time you log in you need to add this line to for example your .bashrc file.
export CGE_RESFINDER_RESGENE_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/db_cge/resfinder"
export CGE_RESFINDER_RESPOINT_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/db_cge/pointfinder"
export CGE_DISINFINDER_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/db_cge/disinfinder"
export CGE_MLST_DB="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/db_cge/mlst"

#export PERL5LIB=/software/miniconda3/envs/dgsp_amr_mlst_detection/lib/perl5/5.32


INPUT_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW"
RUN="LSPV_001_22_M07580"
SAMPLESHEET="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/samplesheet.csv"
IR=${INPUT_PATH}"/"${RUN}
OUTPUT_PATH="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/analysis/"${RUN}
CONDAPATH="/software/miniconda3/envs"
DATABASE_KMERFINDER="/ALMEIDA/PROJECTS/BACTERIAS/DGSP/resources/kmerfinder/databases/bacteria/bacteria"
THREADS=16
SPADESMEM=80

# AMR MLST env
AMR_MLST_ENV=${CONDAPATH}/dgsp_amr_mlst_detection



##########################################################
#
# 01 Create OUTPUT dirs
#
##########################################################

if [[ ! -e ${OUTPUT_PATH} ]]; then
    mkdir -p ${OUTPUT_PATH}"/qc/fastqc/raw"
    mkdir -p ${OUTPUT_PATH}"/qc/fastqc/trim"
    mkdir -p ${OUTPUT_PATH}"/qc/kmerfinder"
    mkdir -p ${OUTPUT_PATH}"/qc/quast"
    mkdir -p ${OUTPUT_PATH}"/qc/stats/trimmomatic"
    mkdir -p ${OUTPUT_PATH}"/out/0_fastq"
    mkdir -p ${OUTPUT_PATH}"/out/1_assembly"
    mkdir -p ${OUTPUT_PATH}"/out/2_amr"
    mkdir -p ${OUTPUT_PATH}"/out/3_annotation"
    mkdir -p ${OUTPUT_PATH}"/out/4_more_analysis/mlst"
    mkdir -p ${OUTPUT_PATH}"/out/4_more_analysis/serotype"
    mkdir -p ${OUTPUT_PATH}"/log/trimmomatic"
    mkdir -p ${OUTPUT_PATH}"/log/bbmap"
    mkdir -p ${OUTPUT_PATH}"/log/kmerfinder"
    mkdir -p ${OUTPUT_PATH}"/log/spades"
    mkdir -p ${OUTPUT_PATH}"/log/quast"
    mkdir -p ${OUTPUT_PATH}"/log/prokka"
    mkdir -p ${OUTPUT_PATH}"/log/resfinder"
    mkdir -p ${OUTPUT_PATH}"/log/mlst_cge"
    mkdir -p ${OUTPUT_PATH}"/log/mlst_torsten"
    mkdir -p ${OUTPUT_PATH}"/log/serotypefinder"
    mkdir -p ${OUTPUT_PATH}"/report"
    elif [[ ! -d ${OUTPUT_PATH} ]]; then
    echo "${OUTPUT_PATH} already exists but is not a directory" 1>&2
fi


        
##########################################################
#
# ANALYSIS
#
##########################################################

# De momento generamos sample sheet
#ls /ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/*fastq.gz | awk -F "/" '{print $8}' | grep "_R1_" | awk -F "_R1_" '{print $1}'
#ls /ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/*fastq.gz | awk -F "/" '{print $8}' | grep "_R1_"
#ls /ALMEIDA/PROJECTS/BACTERIAS/DGSP/RAW/LSPV_001_22_M07580/*fastq.gz | awk -F "/" '{print $8}' | grep "_R2_"

# Con el extra read, eliminamos la primera línea, en teoría los valores son reales
{
    read
    while IFS=, read -r sample fastq_1 fastq_2
    do
        echo "$sample FASTQ1: $fastq_1 FASTQ2: $fastq_2"
        


        ##########################################################
        #
        # 02  Assessment of the genomic sequence quality
        #
        ##########################################################
        sample="STEC_00757_1_17_S20"
        fastq_1="STEC_00757_1_17_S20_1.fastq.gz"
        fastq_2="STEC_00757_1_17_S20_2.fastq.gz "
        
        # RAW FASTQC
        ${CONDAPATH}/dgsp_amr_qc/bin/fastqc --quiet --threads ${THREADS} --outdir ${OUTPUT_PATH}"/qc/fastqc/raw" ${IR}"/"${fastq_1} ${IR}"/"${fastq_2}
        
        # TRIM ADAPTERS
        ${CONDAPATH}/dgsp_amr_qc/bin/trimmomatic PE ${IR}"/"${fastq_1} ${IR}"/"${fastq_2} -threads ${THREADS} \
        ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_1.fastq.gz ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_unpaired_1.fastq.gz  \
        ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_2.fastq.gz ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_unpaired_2.fastq.gz \
        ILLUMINACLIP:/software/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:50 \
        -summary ${OUTPUT_PATH}"/qc/stats/trimmomatic/"${sample}.txt 1> ${OUTPUT_PATH}"/log/trimmomatic/"${sample}.out 2> ${OUTPUT_PATH}"/log/trimmomatic/"${sample}.err
        
        # BBDUK TRIM START / END
        ${CONDAPATH}/dgsp_amr_qc/bin/bbduk.sh -Xmx4g in=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_1.fastq.gz out=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_1.fastq.gz forcetrimleft=14 forcetrimright2=5 minlength=50 1> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_1.out 2> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_1.err
        ${CONDAPATH}/dgsp_amr_qc/bin/bbduk.sh -Xmx4g in=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_2.fastq.gz out=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_2.fastq.gz forcetrimleft=14 forcetrimright2=5 minlength=50 1> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_2.out 2> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_2.err
        ${CONDAPATH}/dgsp_amr_qc/bin/bbduk.sh -Xmx4g in=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_unpaired_1.fastq.gz out=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_unpaired_1.fastq.gz forcetrimleft=14 forcetrimright2=5 minlength=50 1> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_unpaired_1.out 2> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_unpaired_1.err
        ${CONDAPATH}/dgsp_amr_qc/bin/bbduk.sh -Xmx4g in=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_unpaired_2.fastq.gz out=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_unpaired_2.fastq.gz forcetrimleft=14 forcetrimright2=5 minlength=50 1> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_unpaired_2.out 2> ${OUTPUT_PATH}"/log/bbmap/bbduk_"${sample}_unpaired_2.err
        
        # BBDUK Repairing disordered dual files:
        ${CONDAPATH}/dgsp_amr_qc/bin/repair.sh \
        in1=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_1.fastq.gz \
        in2=${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_2.fastq.gz \
        out1=${OUTPUT_PATH}"/out/0_fastq/"${sample}_1.fastq.gz \
        out2=${OUTPUT_PATH}"/out/0_fastq/"${sample}_2.fastq.gz \
        outs=${OUTPUT_PATH}"/out/0_fastq/"${sample}_singletons.fastq.gz repair overwrite=true \
        1> ${OUTPUT_PATH}"/log/bbmap/repair_"${sample}_pair.out 2> ${OUTPUT_PATH}"/log/bbmap/repair_"${sample}_pair.err
        
        # Merge singletons
        cat ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_unpaired_1.fastq.gz ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre_bbmap_unpaired_2.fastq.gz ${OUTPUT_PATH}"/out/0_fastq/"${sample}_singletons.fastq.gz > ${OUTPUT_PATH}"/out/0_fastq/"${sample}_unpaired.fastq.gz
        rm ${OUTPUT_PATH}"/out/0_fastq/"${sample}_pre* ${OUTPUT_PATH}"/out/0_fastq/"${sample}_singletons*
        
        # TRIM FASTQC
        ${CONDAPATH}/dgsp_amr_qc/bin/fastqc \
        --quiet \
        --threads ${THREADS} \
        --outdir ${OUTPUT_PATH}"/qc/fastqc/trim" \
        ${OUTPUT_PATH}"/out/0_fastq/"${sample}_1.fastq.gz ${OUTPUT_PATH}"/out/0_fastq/"${sample}_2.fastq.gz
        
        # KMERFINDER
        #https://bitbucket.org/genomicepidemiology/kmerfinder/src/master/
        ${CONDAPATH}/dgsp_amr_qc/bin/kmerfinder.py \
        -i ${OUTPUT_PATH}"/out/0_fastq/"${sample}_1.fastq.gz ${OUTPUT_PATH}"/out/0_fastq/"${sample}_2.fastq.gz \
        -o ${OUTPUT_PATH}"/qc/kmerfinder/"${sample} \
        -db ${DATABASE_KMERFINDER}.ATG \
        -tax ${DATABASE_KMERFINDER}.tax -x \
        -kp  ${CONDAPATH}/dgsp_amr_qc/bin/ \
        1> ${OUTPUT_PATH}"/log/kmerfinder/"${sample}.out \
        2> ${OUTPUT_PATH}"/log/kmerfinder/"${sample}.err
        
        
        
        ##########################################################
        #
        # 03 Assessment of the genomic sequence quality
        #
        ##########################################################
        
        # ASSEMBLY
        ${CONDAPATH}/dgsp_amr_assembly/bin/spades.py \
        -1 ${OUTPUT_PATH}"/out/0_fastq/"${sample}_1.fastq.gz \
        -2 ${OUTPUT_PATH}"/out/0_fastq/"${sample}_2.fastq.gz \
        -s ${OUTPUT_PATH}"/out/0_fastq/"${sample}_unpaired.fastq.gz \
        -k 21,33,55,77,99,127 \
        --careful --cov-cutoff auto \
        -t ${THREADS} -m ${SPADESMEM} \
        -o ${OUTPUT_PATH}"/out/1_assembly/"${sample} \
        1> ${OUTPUT_PATH}"/log/spades/"${sample}.out \
        2> ${OUTPUT_PATH}"/log/spades/"${sample}.err

        # RENAME FILE
        mv ${OUTPUT_PATH}"/out/1_assembly/"${sample}/contigs.fasta ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta  

        # QUAST
        ${CONDAPATH}/dgsp_amr_qc/bin/quast.py \
        -t ${THREADS} -o ${OUTPUT_PATH}"/qc/assembly_metrics/"${sample} \
        ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta  \
        1> ${OUTPUT_PATH}"/log/quast/"${sample}.out 2> ${OUTPUT_PATH}"/log/quast/"${sample}.err
        

        
        ##########################################################
        #
        # 04 AMR gene and point mutation prediction
        #
        ##########################################################
        
        ##################################################### AMR
        # RESFINDER
        ${CONDAPATH}/dgsp_amr_mlst_detection/bin/run_resfinder.py \
        -db_res ${CGE_RESFINDER_RESGENE_DB} \
        -o ${OUTPUT_PATH}"/out/2_amr/"${sample} \
        -l 0.6 -t 0.8 --acquired \
        -ifa ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta \
        1> ${OUTPUT_PATH}"/log/resfinder/"${sample}.out \
        2> ${OUTPUT_PATH}"/log/resfinder/"${sample}.err
        
        ############################################### ANNOTATION
        # PROKKA
        ${CONDAPATH}/dgsp_amr_mlst_detection/bin/prokka \
        -cpus ${THREADS} --prefix ${sample} \
        --strain ${sample} \
        ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta \
        --outdir ${OUTPUT_PATH}"/out/3_annotation" --force \
        1> ${OUTPUT_PATH}"/log/prokka/"${sample}.out \
        2> ${OUTPUT_PATH}"/log/prokka/"${sample}.err
 
        ##################################################### MLST
        source activate ${AMR_MLST_ENV}
   
        # PARA MLST-CGE necesitamos especificar especie, ayudamos con MLST Torsten Seeman  

        # MLST_Torsten
        mkdir -p ${OUTPUT_PATH}"/out/4_more_analysis/mlst/"${sample}
        ${AMR_MLST_ENV}/bin/perl ${AMR_MLST_ENV}/bin/mlst \
        --label ${sample} \
        ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta \
        1> ${OUTPUT_PATH}"/out/4_more_analysis/mlst/"${sample}"/"${sample}"_mlst.tsv" \
        2> ${OUTPUT_PATH}"/log/mlst_torsten/"${sample}".log"

        # MLST-CGE
        # Multi Locus Sequence Typing (MLST) from an assembled genome or from a set of reads.
        
        # Comprobamos si fichero está vacío
        # Check if the file is empty
        if [ ! -s ${OUTPUT_PATH}"/out/4_more_analysis/mlst/"${sample}"/"${sample}"_mlst.tsv" ]; then
          echo "File MLST Torsten is empty"
        else
          #echo "File is not empty"
          species=`cat ${OUTPUT_PATH}"/out/4_more_analysis/mlst/"${sample}"/"${sample}"_mlst.tsv" | awk '{print $2'} | awk -F "_" '{print $1}'`
          # #${CONDAPATH}/dgsp_amr_detection/bin/python ${CONDAPATH}/dgsp_amr_detection/bin/
          mlst.py \
          -i ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta \
          -o ${OUTPUT_PATH}"/out/4_more_analysis/mlst/"${sample} \
          -s ${species} \
          -d 10 \
          -p ${CGE_MLST_DB} \
          -x \
          1> ${OUTPUT_PATH}"/log/mlst_cge/"${sample}".out" \
          2> ${OUTPUT_PATH}"/log/mlst_cge/"${sample}".err"
          
          ##################################################### IF ECOLI SEROTYPE
          if echo "$species" | grep -q "ecoli"; then
            mkdir -p ${OUTPUT_PATH}"/out/4_more_analysis/serotype/"${sample}
            serotypefinder -i ${OUTPUT_PATH}"/out/1_assembly/"${sample}/${sample}.fasta \
            -o ${OUTPUT_PATH}"/out/4_more_analysis/serotype/"${sample} \
            -p ${CGE_SEROTYPEFINDER_DB} \
            -x \
            1> ${OUTPUT_PATH}"/log/serotypefinder/"${sample}".out" \
            2> ${OUTPUT_PATH}"/log/serotypefinder/"${sample}".err"
          fi 
        fi

        conda deactivate       
        
    done
} < ${SAMPLESHEET}



##########################################################
#
# 05 Multiqc & report
#
##########################################################
${CONDAPATH}/dgsp_amr_qc_multiqc/bin/multiqc ${OUTPUT_PATH}/out ${OUTPUT_PATH}/qc -o ${OUTPUT_PATH}"/report/multiqc" --title ${RUN} --force
