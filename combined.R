TximportFeatureCreate <- function (
                             db=c("Gencode.v33","Gencode.v33lift37","Gencode.v99"),
                             mixed.mouse=c("yes","no"),
                             path.gtf.human,
                             path.gtf.mouse,
                             path.annotation.out) {

db <- match.arg(db)
mixed.mouse <- match.arg(mixed.mouse)

stopifnot(!missing(db), !missing(path.gtf.human), !missing(path.annotation.out), !missing(mixed.mouse))

stopifnot(file.exists(path.gtf.human), length(path.gtf.human) == 1)

library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)
 
options(stringsAsFactors = FALSE)

#make txdb
txdb <- makeTxDbFromGFF(path.gtf.human)
k <- keys(txdb, keytype="TXNAME")
tx2gene_human <- select(txdb, k, "GENEID", "TXNAME")
colnames(tx2gene_human) <- c("transcripts", "genes")

if (db == "Gencode.v33lift37"){
  #re-adding these genes/transcripts, as makeTxDbFromGFF removes them due to them being transcripts with "unmappable stop codons".
  tx2gene_human <- rbind(tx2gene_human, c("ENST00000422803.2_2","ENSG00000119844.15_4"))
  tx2gene_human <- rbind(tx2gene_human, c("ENST00000618549.1_2","ENSG00000197013.10_7"))
  tx2gene_human <- rbind(tx2gene_human, c("ENST00000619291.4_2","ENSG00000048405.10_7"))
  tx2gene_human <- rbind(tx2gene_human, c("ENST00000621077.1_2","ENSG00000122482.21_5"))
  tx2gene_human <- rbind(tx2gene_human, c("ENST00000621229.1_2","ENSG00000131899.11_4"))
  tx2gene_human <- rbind(tx2gene_human, c("ENST00000631326.2_2","ENSG00000105953.15_5"))
}

#isolate gene features
features <- rtracklayer::import(path.gtf.human)
features_gene = as.data.frame(features)
features_gene = features_gene[which(features_gene$type == "gene"),]

if (db == "Gencode.v33"){
  remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","transcript_name","havana_transcript", "ccdsid","ont", "protein_id","exon_id", "exon_number", "transcript_status")
} else if (db == "Gencode.v33lift37"){
  remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","transcript_name","havana_transcript", "ccdsid","ont", "protein_id","exon_id", "exon_number", "transcript_status", "remap_original_location", "remap_substituted_missing_target","transcript_status","remap_target_status","remap_status","remap_num_mappings")
} else {
  remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","tag","transcript_name","transcript_version","transcript_source","havana_transcript","transcript_biotype" ,"ccds_id","protein_version","ont", "protein_id","exon_id", "exon_number", "exon_version")
}
features_gene <- features_gene[,!names(features_gene) %in% remove]
rownames(features_gene) <- features_gene$gene_id
features_gene_human <- features_gene

#isolate transcript features
features_transcript = as.data.frame(features)
features_transcript = features_transcript[which(features_transcript$type == "transcript"),]
if (db == "Gencode.v33"){
  remove <- c("exon_id", "exon_number","exon_version")
} else if (db == "Gencode.v33lift37"){
  remove <- c("exon_id", "exon_number","exon_version", "remap_original_location","remap_target_status","remap_num_mappings","gene_status","remap_substituted_missing_target","transcript_status")
} else {
  remove <- c("exon_id", "exon_number", "exon_version")
}
features_transcript <- features_transcript[,!names(features_transcript) %in% remove]
rownames(features_transcript) <- features_transcript$transcript_id
features_transcript_human <- features_transcript

if(mixed.mouse == "yes"){
  
  txdb <- makeTxDbFromGFF(path.gtf.mouse)
  k <- keys(txdb, keytype="TXNAME")
  tx2gene_mouse <- select(txdb, k, "GENEID", "TXNAME")
  colnames(tx2gene_mouse) <- c("transcripts", "genes")
  
  features <- rtracklayer::import(path.gtf.mouse)
  features_gene = as.data.frame(features)
  
  #isolate gene features
  features_gene = features_gene[which(features_gene$type == "gene"),]
  remove <- c("hgnc_id","transcript_id","transcript_support_level","transcript_type","transcript_name","havana_transcript", "ccdsid","ont", "protein_id","exon_id", "exon_number", "transcript_status","mgi_id", "tag")
  features_gene <- features_gene[,!names(features_gene) %in% remove]
  rownames(features_gene) <- features_gene$gene_id
  features_gene <- features_gene
  
  #isolate transcript features
  features_transcript = as.data.frame(features)
  features_transcript = features_transcript[which(features_transcript$type == "transcript"),]
  remove <- c("exon_id", "exon_number","exon_version", "mgi_id","transcript_support_level")
  features_transcript <- features_transcript[,!names(features_transcript) %in% remove]
  rownames(features_transcript) <- features_transcript$transcript_id
  
  
  #combine human and mouse
  tx2gene <- rbind(tx2gene_human, tx2gene_mouse)
  features_gene <- rbind(
    data.frame(c(features_gene, sapply(setdiff(names(features_gene_human), names(features_gene)), function(x) NA))),
    data.frame(c(features_gene_human, sapply(setdiff(names(features_gene), names(features_gene_human)), function(x) NA)))
  )
  rownames(features_gene) <- features_gene$gene_id
  
  features_transcript <- rbind(
    data.frame(c(features_transcript, sapply(setdiff(names(features_transcript_human), names(features_transcript)), function(x) NA))),
    data.frame(c(features_transcript_human, sapply(setdiff(names(features_transcript), names(features_transcript_human)), function(x) NA)))
  )
  rownames(features_transcript) <- features_transcript$transcript_id
  
  save(tx2gene, features_gene, features_transcript, file=paste0(path.annotation.out,"/", db,"_","mixed.mouse","_annotation",".RData"))
  
} else {

tx2gene <- tx2gene_human
features_gene <- features_gene_human
features_transcript <- features_transcript_human
save(tx2gene, features_gene, features_transcript, file=paste0(path.annotation.out,"/", db,"_annotation",".RData"))
}

}


TximportFeatureCreate(db="Gencode.v33", mixed.mouse = "no", path.gtf.human = "/path/to/gencode.v33.annotation.gtf", path.annotation.out = "/path/to/outdir")
