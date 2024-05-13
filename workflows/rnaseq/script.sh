#!/bin/bash
#import pandas as pd
#df = pd.read_csv('config/sampletable.tsv', sep='\t')
#groups = df.loc[:,'bio_group'].unique()
#subsamples = df[df['bio_group'].str.contains(groups[1])]
#
#if len(subsamples >= 2):
#    
#salmon: 'data/rnaseq_samples/{sample}/{sample}.salmon/quant.sf'

SAMPLETABLE='config/sampletable.tsv'
THREADS=12

INDEX_DIR='references_data/references_data/human/gencode-v28/transcriptome/salmon/human_gencode-v28/'
SAMPLES=`cut -d '	' -f 5 $SAMPLETABLE | uniq | awk 'NR>1'`
    
for group_id in $SAMPLES; do
    OUTPUT_DIR="data/rnaseq_samples/${group_id}.salmon/"
    
    R1files=$(awk -F'\t' -v s="${group_id}" '$5==s {print$7}' $SAMPLETABLE)
    R2files=$(awk -F'\t' -v s="${group_id}" '$5==s {print$8}' $SAMPLETABLE)
    
    LOG="data/rnaseq_samples/${group_id}.salmon/quant.sf.log"
    
    # Create the argument
    fastq_arg="-1 $R1files -2 $R2files"


    echo "salmon quant --index "${INDEX_DIR}" --output "${OUTPUT_DIR}" --threads "${THREADS}" --libType='A' --numGibbsSamples 30 --gcBias --seqBias --validateMappings "${fastq_arg}" &> "${LOG}""
done
