#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: Rscript 1_sample_qc.R analysis_root")
root <- normalizePath(args[1],mustWork=TRUE)
suppressPackageStartupMessages(library(DESeq2))
options(stringsAsFactors=FALSE)
set.seed(20261112)
meta <- read.delim(file.path(root,"inputs/sample_metadata.tsv"),check.names=FALSE)
rownames(meta)<-meta$sample_id
meta$sex<-factor(meta$sex,levels=c("Female","Male"))
meta$breed<-factor(meta$breed,levels=c("CB","TS"))
meta$tissue<-factor(meta$tissue,levels=c("VAF","LIV","GM","BM"))
counts_all<-as.matrix(read.csv(gzfile(file.path(root,"inputs/all_tissues_counts.csv.gz")),row.names=1,check.names=FALSE))
storage.mode(counts_all)<-"integer"
stopifnot(identical(colnames(counts_all),rownames(meta)))

write_gz <- function(x,path,row.names=FALSE) {
    con<-gzfile(path,"wt");on.exit(close(con));write.csv(x,con,row.names=row.names)
}
qc_pca <- function(mat,md,prefix) {
    selected<-rowSums(mat>=10L)>=6L
    dds<-DESeqDataSetFromMatrix(mat[selected,,drop=FALSE],md,~1)
    dds<-estimateSizeFactors(dds)
    vsd<-vst(dds,blind=TRUE,nsub=1000L)
    v<-assay(vsd)
    vars<-matrixStats::rowVars(v)
    chosen<-order(vars,decreasing=TRUE)[seq_len(min(500L,nrow(v)))]
    pc<-prcomp(t(v[chosen,,drop=FALSE]),center=TRUE,scale.=FALSE)
    pct<-100*pc$sdev^2/sum(pc$sdev^2)
    scores<-data.frame(sample_id=rownames(pc$x),pc$x[,seq_len(min(5,ncol(pc$x))),drop=FALSE],check.names=FALSE)
    scores<-merge(md,scores,by="sample_id",sort=FALSE)
    scores<-scores[match(colnames(mat),scores$sample_id),]
    write.csv(scores,file.path(root,"qc",paste0(prefix,"_pca_scores.csv")),row.names=FALSE)
    write.csv(data.frame(component=paste0("PC",seq_along(pct)),variance_percent=pct),file.path(root,"qc",paste0(prefix,"_pca_variance.csv")),row.names=FALSE)
    write.csv(data.frame(feature_id=rownames(v)[chosen],variance=vars[chosen]),file.path(root,"qc",paste0(prefix,"_pca_top500_features.csv")),row.names=FALSE)
    write.csv(data.frame(sample_id=colnames(mat),size_factor=sizeFactors(dds)),file.path(root,"qc",paste0(prefix,"_size_factors.csv")),row.names=FALSE)
    distmat<-as.matrix(dist(t(v)))
    cormat<-cor(v,method="pearson")
    write.csv(distmat,file.path(root,"qc",paste0(prefix,"_sample_distances.csv")))
    write.csv(cormat,file.path(root,"qc",paste0(prefix,"_sample_correlations.csv")))
    write_gz(v,file.path(root,"qc",paste0(prefix,"_vst.csv.gz")),row.names=TRUE)
    saveRDS(list(vsd=vsd,pca=pc,top500=rownames(v)[chosen]),file.path(root,"qc",paste0(prefix,"_qc.rds")))
    message(prefix,": ",sum(selected)," genes for QC; PC1 ",round(pct[1],2),"%; PC2 ",round(pct[2],2),"%")
    list(v=v,vars=vars,dist=distmat,cor=cormat)
}
global<-qc_pca(counts_all,meta,"all_tissues")
top<-order(global$vars,decreasing=TRUE)[seq_len(min(5000L,length(global$vars)))]
v<-global$v[top,,drop=FALSE]
centroid_rows<-list()
for(s in colnames(v)) {
    distances<-sapply(levels(meta$tissue),function(tissue) {
        ids<-meta$sample_id[meta$tissue==tissue & meta$sample_id!=s]
        sqrt(sum((v[,s]-rowMeans(v[,ids,drop=FALSE]))^2))
    })
    corvec<-global$cor[s,];corvec[s]<-NA
    best<-names(which.max(corvec))
    centroid_rows[[s]]<-data.frame(sample_id=s,recorded_tissue=as.character(meta[s,"tissue"]),
        nearest_tissue_centroid=names(which.min(distances)),nearest_other_sample=best,
        nearest_other_sample_tissue=as.character(meta[best,"tissue"]),nearest_other_sample_correlation=max(corvec,na.rm=TRUE),
        distance_to_VAF=distances['VAF'],distance_to_LIV=distances['LIV'],distance_to_GM=distances['GM'],distance_to_BM=distances['BM'])
}
centroids<-do.call(rbind,centroid_rows);rownames(centroids)<-NULL
write.csv(centroids,file.path(root,"qc/tissue_centroid_diagnostics.csv"),row.names=FALSE)
print(centroids[centroids$sample_id=="LIV_TS_16",])
for(tissue in levels(meta$tissue)) {
    md<-droplevels(meta[meta$tissue==tissue,,drop=FALSE])
    within<-qc_pca(counts_all[,rownames(md),drop=FALSE],md,tissue)
    own<-rowMeans(within$dist)
    write.csv(data.frame(sample_id=names(own),mean_distance_to_tissue_samples=own),file.path(root,"qc",paste0(tissue,"_distance_summary.csv")),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(root,"repro/sessionInfo_QC.txt"))
message("QC COMPLETE")
