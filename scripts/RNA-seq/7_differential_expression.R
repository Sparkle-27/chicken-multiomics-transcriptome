#!/usr/bin/env Rscript
# True DESeq2 fits from frozen, unnormalized StringTie/prepDE gene counts.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript 2_differential_expression.R root VAF,LIV,GM,BM")
root<-normalizePath(args[1],mustWork=TRUE)
tissues<-strsplit(args[2],",",fixed=TRUE)[[1]]
options(stringsAsFactors=FALSE)
set.seed(20261112)
suppressPackageStartupMessages(library(DESeq2))
meta<-read.delim(file.path(root,"inputs/sample_metadata.tsv"),check.names=FALSE)
rownames(meta)<-meta$sample_id
meta$sex<-factor(meta$sex,levels=c("Female","Male"))
meta$breed<-factor(meta$breed,levels=c("CB","TS"))
ann<-read.delim(file.path(root,"inputs/gene_annotation.tsv"),check.names=FALSE)
config<-jsonlite::read_json(file.path(root,"inputs/configuration.json"),simplifyVector=TRUE)
stopifnot(!any(meta$sample_id %in% config$main_exclusions))

write_gz<-function(x,path) {
    con<-gzfile(path,"wt");on.exit(close(con));write.csv(x,con,row.names=TRUE)
}
fit_one<-function(ct,md,tissue,label,weight=FALSE) {
    started<-Sys.time()
    dest<-file.path(root,"results",tissue,label)
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    stopifnot(identical(colnames(ct),rownames(md)))
    md<-droplevels(md)
    md$weight_z<-as.numeric(scale(md$weight_g))
    formula <- if(weight) ~sex+weight_z+breed else ~sex+breed
    mm<-model.matrix(formula,md)
    stopifnot(qr(mm)$rank==ncol(mm),nrow(mm)>ncol(mm))
    keep<-rowSums(ct>=config$min_count)>=config$min_samples
    stopifnot(sum(keep)>0L,identical(rownames(ct),ann$feature_id))
    message(format(Sys.time())," START ",tissue," ",label," n=",nrow(md)," retained=",sum(keep)," ",deparse(formula))
    dds<-DESeqDataSetFromMatrix(ct[keep,,drop=FALSE],md,formula)
    dds<-DESeq(dds,test="Wald",fitType="parametric",sfType="ratio",betaPrior=FALSE,parallel=FALSE,quiet=FALSE)
    initial_nonconverged<-sum(!mcols(dds)$betaConv,na.rm=TRUE)
    if(initial_nonconverged>0L) {
        message("Refitting Wald coefficients with maxit=1000 for convergence")
        dds<-nbinomWaldTest(dds,betaPrior=FALSE,maxit=1000,quiet=FALSE)
    }
    res<-results(dds,contrast=c("breed","TS","CB"),alpha=config$alpha,
                 pAdjustMethod="BH",independentFiltering=TRUE)
    rdf<-as.data.frame(res)
    full<-ann
    full$prefilter_pass<-keep
    for(k in names(rdf))full[[k]]<-rdf[[k]][match(full$feature_id,rownames(rdf))]
    full$log2FC_CI95_low<-full$log2FoldChange-1.95996398454*full$lfcSE
    full$log2FC_CI95_high<-full$log2FoldChange+1.95996398454*full$lfcSE
    full$beta_converged<-mcols(dds)$betaConv[match(full$feature_id,rownames(dds))]
    full$fdr_significant<-!is.na(full$padj)&full$padj<config$alpha
    full$direction<-ifelse(full$fdr_significant,ifelse(full$log2FoldChange>0,"TS_higher","CB_higher"),"not_significant")
    full$test_status<-ifelse(!keep,"prefiltered",ifelse(is.na(full$pvalue),"no_pvalue_count_outlier_or_unavailable",
                                 ifelse(is.na(full$padj),"independent_filtered","tested_with_padj")))
    full$max_Cooks_distance<-if("maxCooks" %in% names(mcols(dds)))mcols(dds)$maxCooks[match(full$feature_id,rownames(dds))] else NA_real_
    full$counts_replaced_for_fit<-if("replace" %in% names(mcols(dds)))mcols(dds)$replace[match(full$feature_id,rownames(dds))] else NA
    full<-full[order(full$padj,full$pvalue,full$feature_id,na.last=TRUE),]
    write.csv(full,file.path(dest,"results_all.csv"),row.names=FALSE,na="")
    write.csv(full[full$fdr_significant,,drop=FALSE],file.path(dest,"results_FDR_0.05.csv"),row.names=FALSE,na="")
    for(direction in c("TS_higher","CB_higher"))write.csv(full[full$direction==direction,,drop=FALSE],file.path(dest,paste0(direction,"_FDR_0.05.csv")),row.names=FALSE,na="")
    write.csv(full[!is.na(full$pvalue)&full$pvalue<0.05&!full$fdr_significant,,drop=FALSE],file.path(dest,"nominal_only_exploratory.csv"),row.names=FALSE,na="")
    write_gz(counts(dds,normalized=TRUE),file.path(dest,"normalized_counts.csv.gz"))
    write.csv(data.frame(sample_id=colnames(dds),size_factor=sizeFactors(dds)),file.path(dest,"size_factors.csv"),row.names=FALSE)
    write.table(md,file.path(dest,"samples_used.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
    dg<-data.frame(feature_id=rownames(dds),baseMean=rowMeans(counts(dds,normalized=TRUE)),
        dispersion_gene=mcols(dds)$dispGeneEst,dispersion_fit=mcols(dds)$dispFit,
        dispersion_MAP=dispersions(dds),beta_converged=mcols(dds)$betaConv)
    write.csv(dg,file.path(dest,"gene_fit_diagnostics.csv"),row.names=FALSE)
    cooks<-assays(dds)[["cooks"]]
    csummary<-data.frame(sample_id=colnames(dds),median_Cooks=apply(cooks,2,median,na.rm=TRUE),
        q95_Cooks=apply(cooks,2,quantile,probs=.95,na.rm=TRUE),max_Cooks=apply(cooks,2,max,na.rm=TRUE),
        genes_Cooks_gt1=colSums(cooks>1,na.rm=TRUE))
    write.csv(csummary,file.path(dest,"sample_Cooks_diagnostics.csv"),row.names=FALSE)
    saveRDS(dds,file.path(dest,"dds_fitted.rds"))
    saveRDS(res,file.path(dest,"DESeqResults.rds"))
    info<-list(tissue=tissue,analysis=label,formula=deparse(formula),contrast="TS / CB",n_samples=nrow(md),
        CB_n=unname(sum(md$breed=="CB")),TS_n=unname(sum(md$breed=="TS")),
        original_features=nrow(ct),prefilter_retained=sum(keep),pvalue_available=sum(!is.na(full$pvalue)),
        padj_available=sum(!is.na(full$padj)),independent_filtered=sum(full$prefilter_pass&!is.na(full$pvalue)&is.na(full$padj)),
        no_pvalue=sum(full$prefilter_pass&is.na(full$pvalue)),significant=sum(full$fdr_significant),
        TS_higher=sum(full$direction=="TS_higher"),CB_higher=sum(full$direction=="CB_higher"),
        initial_nonconverged=initial_nonconverged,final_nonconverged=sum(!mcols(dds)$betaConv,na.rm=TRUE),
        replaced_features=if("replace" %in% names(mcols(dds)))sum(mcols(dds)$replace,na.rm=TRUE) else 0,
        model_matrix_rank=qr(mm)$rank,model_columns=ncol(mm),model_kappa=kappa(mm,exact=TRUE),
        independent_filter_baseMean_threshold=unname(metadata(res)$filterThreshold),
        weight_mean_g=mean(md$weight_g),weight_sd_g=sd(md$weight_g),
        weight_z_VIF=if(weight)1/(1-summary(lm(weight_z~sex+breed,data=md))$r.squared) else NA_real_,
        coefficient_names=resultsNames(dds),contrast_description=mcols(res)$description[2],
        DESeq2_version=as.character(packageVersion("DESeq2")),R_version=R.version.string,
        elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
        samples=md$sample_id)
    jsonlite::write_json(info,file.path(dest,"model_summary.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
    writeLines(capture.output(summary(res)),file.path(dest,"DESeq2_summary.txt"))
    message(format(Sys.time())," DONE ",tissue," ",label," significant=",sum(full$fdr_significant)," TS_high=",info$TS_higher," CB_high=",info$CB_higher)
    invisible(info)
}

for(tissue in tissues) {
    all_ct<-as.matrix(read.csv(gzfile(file.path(root,"inputs",paste0(tissue,"_counts.csv.gz"))),row.names=1,check.names=FALSE))
    storage.mode(all_ct)<-"integer"
    all_md<-meta[colnames(all_ct),,drop=FALSE]
    ids<-colnames(all_ct)
    ct<-all_ct[,ids,drop=FALSE];md<-all_md[ids,,drop=FALSE]
    fit_one(ct,md,tissue,"primary",FALSE)
    fit_one(ct,md,tissue,"weight_adjusted",TRUE)
}
writeLines(capture.output(sessionInfo()),file.path(root,"repro/sessionInfo_DESeq2.txt"))
message("REQUESTED FITS COMPLETE: ",paste(tissues,collapse=","))
