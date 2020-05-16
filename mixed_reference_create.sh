#used to create human-mouse mixed reference for kallisto and salmon

cut -d '|' -f1 gencode.v33.transcripts.fa > gencode.v33.transcripts_headerTrimmed.fa
cut -d '|' -f1 gencode.vM24.transcripts.fa > gencode.vM24.transcripts_headerTrimmed.fa

cat gencode.v32.transcripts_headerTrimmed.fa gencode.vM23.transcripts_headerTrimmed.fa > mixed.fa

kallisto index -i human.idx gencode.v33.transcripts_headerTrimmed.fa
kallisto index -i human_mouse.idx mixed.fa

xenome index -M 15 -T 8 -P XenomeIdx -H GRCm38.primary_assembly.genome.fa -G GRCh38.primary_assembly.genome.fa

xenome classify -M 15 -T 8 -l log.txt -v -i "$SampleName"_R1.fastq -i "$SampleName"_R2.fastq --pairs --output-filename-prefix "$SampleName" -P "$xenomeRef" > "$SampleName"_stat.txt

cmd = 'kallisto quant -t 8 -i ' + KALLISTO_HUMAN_INDEX + ' -o ' + OUTPUT_FOLDER + ' ' + HUMAN_FASTQ_FILE_1 + ' ' + HUMAN_FASTQ_FILE_2

cmd = 'kallisto quant --pseudobam -t 8 -i ' + KALLISTO_MIXED_INDEX + ' -o ' + OUTPUT_FOLDER + ' ' + MIXED_FASTQ_FILE_1 + ' ' + MIXED_FASTQ_FILE_2

cmd = SALMON_PATH + ' quant -i /media/soheil/PDX_Disk2/RefGenomes/SalmonIndex_human -l A -1 ' + HUMAN_FASTQ_FILE_1 + ' -2 ' + HUMAN_FASTQ_FILE_2 + ' --validateMappings -o Output_human_salmon'

cmd = SALMON_PATH + ' quant -i /media/soheil/PDX_Disk2/RefGenomes/SalmonIndex_mixed -l A -1 ' + MIXED_FASTQ_FILE_1 + ' -2 ' + MIXED_FASTQ_FILE_2 + ' --validateMappings -o Output_mixed_salmon'

cmd = 'cat XenomeOutput__graft_1.fastq XenomeOutput__both_1.fastq XenomeOutput__ambiguous_1.fastq > No_Mouse_1.fastq'
cmd = 'cat XenomeOutput__graft_2.fastq XenomeOutput__both_2.fastq XenomeOutput__ambiguous_2.fastq > No_Mouse_2.fastq'


head -227913 Output_mixed/abundance.tsv > Output_mixed/abundance_human.tsv
