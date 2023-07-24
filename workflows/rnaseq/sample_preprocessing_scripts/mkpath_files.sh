#!/bin/bash

find /data/CARD_ARDIS/2012_G3boydenRNAseq/raw_data/fastq210201/  -maxdepth 1 -name "*.fastq.gz" > all_lanes_samples.tsv
find /data/CARD_ARDIS/TDP43kdRNAseq/NovaSeq221007notrebalanced_JR7961/ -maxdepth 1 -name "*.fastq.gz"  >> all_lanes_samples.tsv
find /data/CARD_ARDIS/TDP43kdRNAseq/NovaSeq221007notrebalanced_JR7961/sample1320/ -maxdepth 2 -name "*.fastq.gz"  >> all_lanes_samples.tsv
find /data/CARD_ARDIS/FUSkdRNAseq/NovaSeq_VR8153/ -maxdepth 1 -name "*.fastq.gz"  >> all_lanes_samples.tsv
find /data/CARD_ARDIS/HNRNPA1kdRNAseq/NovaSeq_VR8147/ -maxdepth 1 -name "*.fastq.gz"  >> all_lanes_samples.tsv
