library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)

options(stringsAsFactors = FALSE)

#human tx2gene
txdb <- makeTxDbFromGFF("/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33/gencode.v33.annotation.gtf")
saveDb(txdb, "/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33/gencode.v33.annotation.txdb")

txdb <- loadDb("/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33/gencode.v33.annotation.txdb")
k <- keys(txdb, keytype="TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene) <- c("transcripts", "genes")

#mouse tx2gene

txdb <- makeTxDbFromGFF("/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/mouse/gencode.vM24.annotation.gtf")
saveDb(txdb, "/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/mouse/gencode.vM24.annotation.txdb")

txdb <- loadDb("/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/mouse/gencode.vM24.annotation.txdb")
k <- keys(txdb, keytype="TXNAME")
tx2gene_mouse <- select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene_mouse) <- c("transcripts", "genes")

#combine human and mouse tx2gene

tx2gene <- rbind(tx2gene, tx2gene_mouse)


#feature data human
features <- rtracklayer::import('/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/gen_v33/gencode.v33.annotation.gtf')

features_gene = as.data.frame(features)
features_gene = features_gene[which(features_gene$type == "gene"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","transcript_name","havana_transcript", "ccdsid","ont", "protein_id","exon_id", "exon_number", "transcript_status")
features_gene <- features_gene[,!names(features_gene) %in% remove]
rownames(features_gene) <- features_gene$gene_id
features_gene_human <- features_gene

features_transcript = as.data.frame(features)
features_transcript = features_transcript[which(features_transcript$type == "transcript"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("exon_id", "exon_number","exon_version", "hgnc_id")
features_transcript <- features_transcript[,!names(features_transcript) %in% remove]
rownames(features_transcript) <- features_transcript$transcript_id
features_transcript_human <- features_transcript

#feature data mouse


features <- rtracklayer::import('/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/mouse/gencode.vM24.annotation.gtf')

features_gene = as.data.frame(features)
features_gene = features_gene[which(features_gene$type == "gene"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","transcript_name","havana_transcript", "ccdsid","ont", "protein_id","exon_id", "exon_number", "transcript_status","mgi_id")
features_gene <- features_gene[,!names(features_gene) %in% remove]
rownames(features_gene) <- features_gene$gene_id
features_gene_mouse <- features_gene

features_transcript = as.data.frame(features)
features_transcript = features_transcript[which(features_transcript$type == "transcript"),]
#features_df <- features_df[order(match(features_df$transcript_id,tx2gene$TXNAME)),]
remove <- c("exon_id", "exon_number","exon_version", "mgi_id")
features_transcript <- features_transcript[,!names(features_transcript) %in% remove]
rownames(features_transcript) <- features_transcript$transcript_id
features_transcript_mouse <- features_transcript


#combine human and mouse
features_gene <- rbind(features_gene_human, features_gene_mouse)
features_transcript <- rbind(features_transcript_human, features_transcript_mouse)

save(tx2gene, features_gene, features_transcript, file="/Users/anthmam/Desktop/Projects/BHK/ORCESTRA/Transcriptome/mouse/Gencode.v33.m24.mouse.mixed.annotation.RData")
