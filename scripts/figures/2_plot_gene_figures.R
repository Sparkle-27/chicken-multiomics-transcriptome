#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript 2_plot_gene_figures.R analysis_root")
root<-normalizePath(args[1],mustWork=TRUE)
suppressPackageStartupMessages({library(ggplot2);library(ggrepel);library(pheatmap)})
options(stringsAsFactors=FALSE)
set.seed(20261112)
out<-file.path(root,"figures/pdf")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
meta<-read.delim(file.path(root,"inputs/sample_metadata.tsv"),check.names=FALSE)
tissues<-c("VAF","LIV","GM","BM")
breed_colors<-c(CB="#5CA7BA",TS="#FE4365")
de_colors<-c(CB_higher="#5CA7BA",not_significant="#CECECE",TS_higher="#FE4365")
theme_set(theme_classic(base_size=11,base_family="Helvetica")+
    theme(plot.title=element_text(size=12,face="bold"),legend.title=element_blank(),
          legend.key=element_blank(),axis.title=element_text(color="black"),
          axis.text=element_text(color="black"),plot.margin=margin(9,12,9,9)))
files<-character()
saveplot<-function(p,name,w=5.6,h=4.7){
    path<-file.path(out,paste0(name,".pdf"))
    grDevices::pdf(path,width=w,height=h,family="Helvetica",useDingbats=FALSE)
    print(p);dev.off();files<<-c(files,path)
}
savebase<-function(name,expr,w=8.2,h=6.2){
    path<-file.path(out,paste0(name,".pdf"));grDevices::pdf(path,width=w,height=h,family="Helvetica",useDingbats=FALSE)
    tryCatch(force(expr),finally=dev.off());files<<-c(files,path)
}
draw_panels<-function(plots,nrow=1,ncol=length(plots)){
    grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow,ncol)))
    for(i in seq_along(plots))print(plots[[i]],vp=grid::viewport(layout.pos.row=ceiling(i/ncol),layout.pos.col=(i-1)%%ncol+1))
}
read_result<-function(t,m="primary")read.csv(file.path(root,"results",t,m,"results_all.csv"),check.names=FALSE)
summary<-read.csv(file.path(root,"tables/model_summary.csv"))
primary<-summary[summary$analysis=="primary",]
primary<-primary[match(tissues,primary$tissue),]

# Main figures: breed contrasts and global sample structure.
long<-rbind(data.frame(tissue=primary$tissue,direction="CB higher",genes=primary$CB_higher),
            data.frame(tissue=primary$tissue,direction="TS higher",genes=primary$TS_higher))
long$tissue<-factor(long$tissue,levels=tissues)
p<-ggplot(long,aes(tissue,genes,fill=direction))+geom_col(width=.64)+
    scale_fill_manual(values=c("CB higher"="#5CA7BA","TS higher"="#FE4365"))+
    labs(x=NULL,y="Genes with adjusted P < 0.05")+theme(legend.position="top")
totals<-data.frame(tissue=factor(primary$tissue,levels=tissues),genes=primary$significant)
p<-p+geom_text(data=totals,aes(tissue,genes,label=format(genes,big.mark=",")),inherit.aes=FALSE,vjust=-.6,size=3.5)+
    scale_y_continuous(expand=expansion(mult=c(0,.12)))
saveplot(p,"01_DEG_counts",6.0,4.7)

# Exact color/pch mapping from the supplied 32-uMAP.r, joined by sample ID.
# Original AA_F/AA_L/AA_G/AA_P = CB_VAF/CB_LIV/CB_GM/CB_BM; WJ = TS.
pca_style<-data.frame(original_group=c("AA_F","AA_L","AA_G","AA_P","WJ_F","WJ_L","WJ_G","WJ_P"),
    group=paste(rep(c("CB","TS"),each=4),rep(tissues,2),sep="_"),
    color=rep(c(rgb(246,112,136,max=255),rgb(150,163,49,max=255),
                rgb(53,172,164,max=255),rgb(163,140,244,max=255)),2),
    pch=c(3,5,10,12,15,16,23,25))
write.csv(pca_style,file.path(root,"figures/PCA_style_mapping.csv"),row.names=FALSE)
pca_plot<-function(prefix,name){
    df<-read.csv(file.path(root,"qc",paste0(prefix,"_pca_scores.csv")))
    stopifnot(!"LIV_TS_16" %in% df$sample_id)
    var<-read.csv(file.path(root,"qc",paste0(prefix,"_pca_variance.csv")))$variance_percent
    groups<-paste(df$breed,df$tissue,sep="_")
    ix<-match(groups,pca_style$group);stopifnot(!anyNA(ix))
    legends<-which(pca_style$group %in% groups)
    savebase(name,{
        layout(matrix(1:2,nrow=1),widths=c(1,0.35),height=2)
        par(mai=c(0.7,0.7,0.1,0.1),omi=c(0,0,0,0),mgp=c(2,0.5,0),tcl=-0.3)
        plot(df$PC1,df$PC2,xlab=sprintf("PC1 (%.1f%%)",var[1]),ylab=sprintf("PC2 (%.1f%%)",var[2]),
             col=pca_style$color[ix],pch=pca_style$pch[ix],cex=1.5,lwd=1.2,cex.lab=1.2,cex.axis=1.2)
        par(mai=c(0.7,0,0.1,0.1))
        plot(x=0,y=0,axes=FALSE,type="n",col="white",xlab="",ylab="")
        legend("topleft",legend=pca_style$group[legends],col=pca_style$color[legends],
            pch=pca_style$pch[legends],bty="n",cex=1,text.font=1)
    },8,6)
}
pca_plot("all_tissues","02_all_tissues_PCA")
for(tissue in tissues)pca_plot(tissue,paste0("03_",tissue,"_PCA"))

for(tissue in tissues){
    df<-read_result(tissue)
    pd<-df[is.finite(df$log2FoldChange)&is.finite(df$pvalue),]
    pd$Significant<-factor(ifelse(pd$fdr_significant,ifelse(pd$log2FoldChange>0,"Up","Down"),"Stable"),
                          levels=c("Down","Stable","Up"))
    pd$plot_logp<--log10(pmax(pd$pvalue,1e-300))
    # Match original figure dimensions, geometry, colors, labels, legend and typography.
    # Updated statistical classification: BH padj < .05, no post hoc fold-change cutoff.
    # Current analysis contrast is TS/CB; Up=TS higher, Down=CB higher.
    p<-ggplot(pd,aes(x=log2FoldChange,y=plot_logp))+geom_point(aes(color=Significant),size=2)+
        scale_color_manual(values=c(Down="#5CA7BA",Stable="grey",Up="#FE4365"),drop=FALSE)+
        labs(x="log2(fold change)",y="-log10 (p-value)")+theme_bw()+
        theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
            legend.position="bottom",legend.key.width=grid::unit(0.8,"cm"),legend.key.height=grid::unit(0.8,"cm"),
            axis.text.x=element_text(face="bold",colour="black",angle=0,vjust=0.5,hjust=1,size=18),
            axis.text.y=element_text(face="bold",colour="black",size=18),
            axis.title.x=element_text(face="bold",colour="black",size=24),
            axis.title.y=element_text(face="bold",colour="black",size=24),
            axis.line.x=element_line(colour="black",linewidth=0.5),axis.line.y=element_line(colour="black",linewidth=0.5),
            legend.text=element_text(face="bold",colour="black",size=18),
            legend.title=element_text(face="bold",colour="black",size=18),
            legend.background=element_rect(fill="transparent",colour=NA),
            plot.margin=margin(5.5,5.5,5.5,5.5))
    saveplot(p,paste0("04_",tissue,"_volcano"),8,8)
}

for(tissue in tissues){
    df<-read_result(tissue)
    df<-df[df$prefilter_pass&is.finite(df$log2FoldChange)&df$baseMean>0,]
    df$direction<-factor(df$direction,levels=names(de_colors))
    p<-ggplot(df,aes(baseMean,log2FoldChange))+geom_point(aes(color=direction),size=.85,alpha=.65)+
        scale_x_log10()+scale_color_manual(values=de_colors,labels=c("CB higher","Not significant","TS higher"),drop=FALSE)+
        geom_hline(yintercept=0,linewidth=.3,color="#777777")+
        labs(x="Mean normalized count",y=expression(log[2](TS/CB)),title=tissue)+
        theme(legend.position="top",legend.text=element_text(size=9))
    saveplot(p,paste0("05_",tissue,"_MA"),6.1,4.8)
}

ints<-read.csv(file.path(root,"tables/DEG_overlap_intersections.csv"))
ints<-ints[order(ints$genes,decreasing=TRUE),]
savebase("06_DEG_overlap",{
    layout(matrix(c(1,2),2,1),heights=c(3.5,1.6))
    par(mar=c(.3,4.4,2,1))
    mid<-barplot(ints$genes,col="#474D58",border=NA,axes=FALSE,ylim=c(0,max(ints$genes)*1.16))
    axis(2,las=1);mtext("Genes in each exact intersection",2,2.8)
    text(mid,ints$genes,labels=ints$genes,pos=3,cex=.7)
    par(mar=c(2,4.4,.2,1))
    plot(NA,xlim=c(min(mid)-.6,max(mid)+.6),ylim=c(.5,4.5),axes=FALSE,xlab="",ylab="")
    axis(2,at=4:1,labels=tissues,las=1,tick=FALSE)
    for(j in seq_len(nrow(ints))){
        vals<-vapply(ints[j,tissues],function(x)tolower(as.character(x))%in%c("true","1"),logical(1));points(rep(mid[j],4),4:1,pch=16,col="#DEDEDE",cex=1.05)
        yy<-(4:1)[vals]
        if(length(yy)>1)segments(mid[j],min(yy),mid[j],max(yy),lwd=1.7,col="#474D58")
        points(rep(mid[j],length(yy)),yy,pch=16,col="#474D58",cex=1.05)
    }
},8.3,5.4)

# Sample and model QC are kept separate from the publication contrasts.
qc<-read.delim(file.path(root,"qc/raw_count_sample_qc.tsv"))
qc$tissue<-factor(qc$tissue,levels=tissues)
p1<-ggplot(qc,aes(tissue,library_gene_count_sum/1e6,color=breed))+geom_jitter(width=.14,height=0,size=1.8,alpha=.8)+
    scale_color_manual(values=breed_colors)+labs(x=NULL,y="Gene count sum (millions)")+theme(legend.position="top")
p2<-ggplot(qc,aes(tissue,detected_genes_count_ge1,color=breed))+geom_jitter(width=.14,height=0,size=1.8,alpha=.8)+
    scale_color_manual(values=breed_colors)+labs(x=NULL,y="Genes with count >= 1")+theme(legend.position="top")
savebase("07_sample_count_QC",draw_panels(list(p1,p2)),10.4,4.7)

cent<-read.csv(file.path(root,"qc/tissue_centroid_diagnostics.csv"))
own<-sapply(seq_len(nrow(cent)),function(i)cent[i,paste0("distance_to_",cent$recorded_tissue[i])])
other<-sapply(seq_len(nrow(cent)),function(i)min(unlist(cent[i,paste0("distance_to_",setdiff(tissues,cent$recorded_tissue[i]))])))
cent$other_minus_own<-other-own
cent<-merge(cent,meta[,c("sample_id","breed")],by="sample_id",sort=FALSE)
cent$recorded_tissue<-factor(cent$recorded_tissue,levels=rev(tissues))
p1<-ggplot(cent,aes(other_minus_own,recorded_tissue,color=breed))+geom_jitter(height=.1,width=0,size=1.8,alpha=.75)+
    geom_vline(xintercept=0,linetype="dashed",linewidth=.35)+scale_color_manual(values=breed_colors)+
    labs(x="Nearest other-tissue distance - own-tissue distance",y=NULL)+theme(legend.position="top")
saveplot(p1,"08_tissue_identity_QC",8.2,4.6)

for(tissue in tissues){
    prefix<-tissue
    dm<-as.matrix(read.csv(file.path(root,"qc",paste0(prefix,"_sample_distances.csv")),row.names=1,check.names=FALSE))
    md<-meta[match(colnames(dm),meta$sample_id),];rownames(md)<-md$sample_id
    at<-md[,c("breed","sex"),drop=FALSE]
    ac<-list(breed=breed_colors,sex=c(Female="#8D86AF",Male="#C2B7A2"))
    hc<-hclust(as.dist(dm))
    ph<-pheatmap(dm,cluster_rows=hc,cluster_cols=hc,annotation_col=at,annotation_colors=ac,
        color=colorRampPalette(c("#FFFFFF","#B7CFDC","#263F59"))(100),border_color=NA,
        show_colnames=FALSE,fontsize_row=6.8,main=tissue,silent=TRUE)
    savebase(paste0("09_",tissue,"_sample_distance"),{grid::grid.newpage();grid::grid.draw(ph$gtable)},7.6,7.1)
}

for(tissue in tissues){
    d<-read.csv(file.path(root,"results",tissue,"primary/gene_fit_diagnostics.csv"))
    co<-read.csv(file.path(root,"results",tissue,"primary/sample_Cooks_diagnostics.csv"))
    r<-read_result(tissue)
    md<-meta[match(co$sample_id,meta$sample_id),]
    savebase(paste0("10_",tissue,"_fit_diagnostics"),{
        par(mfrow=c(2,2),mar=c(4.2,4.3,2,1),cex=.86)
        valid<-d$baseMean>0&d$dispersion_gene>0
        plot(d$baseMean[valid],d$dispersion_gene[valid],log="xy",pch=16,cex=.24,col=adjustcolor("#555555",.25),xlab="Mean normalized count",ylab="Dispersion",main=tissue)
        o<-order(d$baseMean);lines(d$baseMean[o],d$dispersion_fit[o],col="#C83752",lwd=1.2)
        plot(r$baseMean,r$lfcSE,log="xy",pch=16,cex=.25,col=adjustcolor("#444444",.3),xlab="Mean normalized count",ylab="LFC standard error")
        barplot(co$q95_Cooks,names.arg=co$sample_id,las=2,cex.names=.4,col=breed_colors[md$breed],border=NA,ylab="95th percentile Cook's distance")
        hist(r$pvalue,breaks=40,col="#8AA6B9",border="white",xlab="Nominal P value",main="",ylab="Genes")
    },9.4,7.4)
}

for(tissue in tissues){
    d<-read.csv(file.path(root,"tables",paste0(tissue,"_weight_sensitivity_gene_comparison.csv")))
    d<-d[is.finite(d$log2FoldChange_primary)&is.finite(d$log2FoldChange_alternative),]
    a<-!is.na(d$padj_primary)&d$padj_primary<.05;b<-!is.na(d$padj_alternative)&d$padj_alternative<.05
    d$class<-factor(ifelse(a&b,"Both",ifelse(a,"Primary only",ifelse(b,"Weight only","Neither"))),levels=c("Neither","Primary only","Weight only","Both"))
    p<-ggplot(d,aes(log2FoldChange_primary,log2FoldChange_alternative,color=class))+
        geom_point(size=.85,alpha=.55)+geom_hline(yintercept=0,color="#999999",linewidth=.25)+
        geom_vline(xintercept=0,color="#999999",linewidth=.25)+geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=.35)+
        scale_color_manual(values=c(Neither="#D3D3D3","Primary only"="#445F78","Weight only"="#C8A451",Both="#8B4478"),drop=FALSE)+
        labs(x="Primary-model log2(TS/CB)",y="Weight-adjusted log2(TS/CB)",title=tissue)+
        theme(legend.position="top",legend.text=element_text(size=8.4))
    saveplot(p,paste0("11_",tissue,"_weight_sensitivity"),6.2,5.0)
}
stopifnot(length(files)==29L)
writeLines(files,file.path(root,"figures/pdf_order.txt"))
message("PDF FIGURES COMPLETE: ",length(files))
if(!is.null(warnings()))writeLines(capture.output(warnings()),file.path(root,"logs/05_figure_warnings.txt"))
