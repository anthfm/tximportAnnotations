library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)

options(stringsAsFactors = FALSE)

txdb <- makeTxDbFromGFF("/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33_lift/gencode.v33lift37.annotation.gtf")
saveDb(txdb, "/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33_lift/gencode.v33lift37.annotation.txdb")

txdb <- loadDb("/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33_lift/gencode.v33lift37.annotation.txdb")
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene) <- c("transcripts", "genes")
#re-adding these transcripts, as makeTxDbFromGFF removes them due to them being transcripts with "unmappable stop codons".
tx2gene <- rbind(tx2gene, c("ENST00000422803.2_2","ENSG00000119844.15_4"))
tx2gene <- rbind(tx2gene, c("ENST00000618549.1_2","ENSG00000197013.10_7"))
tx2gene <- rbind(tx2gene, c("ENST00000619291.4_2","ENSG00000048405.10_7"))
tx2gene <- rbind(tx2gene, c("ENST00000621077.1_2","ENSG00000122482.21_5"))
tx2gene <- rbind(tx2gene, c("ENST00000621229.1_2","ENSG00000131899.11_4"))
tx2gene <- rbind(tx2gene, c("ENST00000631326.2_2","ENSG00000105953.15_5"))

features <- rtracklayer::import('/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33_lift/gencode.v33lift37.annotation.gtf')

features_gene = as.data.frame(features)
features_gene = features_gene[which(features_gene$type == "gene"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","transcript_name","havana_transcript", "ccdsid","ont", "protein_id","exon_id", "exon_number", "transcript_status", "remap_original_location")
features_gene <- features_gene[,!names(features_gene) %in% remove]
rownames(features_gene) <- features_gene$gene_id


features_transcript = as.data.frame(features)
features_transcript = features_transcript[which(features_transcript$type == "transcript"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("exon_id", "exon_number","exon_version", "remap_original_location")
features_transcript <- features_transcript[,!names(features_transcript) %in% remove]
rownames(features_transcript) <- features_transcript$transcript_id
#remove_transcripts <- unique(features_transcript$transcript_id)[which(!unique(features_transcript$transcript_id) %in% unique(tx2gene$transcripts))]
#features_transcript <- features_transcript[-which(features_transcript$transcript_id %in% remove_transcripts),]


save(tx2gene, features_gene, features_transcript, file="/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33_lift/Gencode.v33lift37.annotation.RData")











