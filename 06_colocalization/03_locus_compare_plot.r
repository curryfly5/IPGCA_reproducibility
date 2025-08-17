library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
### nohup Rscript coloc_func.r --gwas /disk201/chenzt/metaGWAS/M_BFT.txt.gz --eqtl /disk191_2/yupf/GUT_RNA/04_SNPCALLING/case/coloc/New_script/test_eQTL.txt --tis ileum_ent --trt M_BFT --outdir /disk191_2/yupf/GUT_RNA/04_SNPCALLING/case/coloc/New_script > coloc_func.log 2>&1 &
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  arg_list <- list()
  for (i in seq(1, length(args), by = 2)) {
    key <- args[i]
    value <- args[i + 1]
    arg_list[[key]] <- value
  }
  return(arg_list)
}
params <- parse_args(args)
gwasfile<-params[["--gwas"]]
eqtlfile<-params[["--eqtl"]]
gene<- params[["--gene"]]
gene_name <- if (!is.null(params[["--name"]])) {
  params[["--name"]]
} else {
  gene
}
target_snp <- params[["--snp"]]
trait <- params[["--trt"]]
tissue <- params[["--tiss"]]
outdir<-params[["--outdir"]]
dir.create(outdir,showWarnings = F)

if (grepl("_", target_snp))  {
  target_pos <- strsplit(as.character(target_snp), "_")[[1]][2]
}
if (grepl(":", target_snp))  {
  target_pos <- strsplit(as.character(target_snp), ":")[[1]][2]
  target_snp <- gsub(pattern = ":", replacement = "_", target_snp)
}

gene=gene_name
#format gwas
gwas<-read.table(gwasfile,header = TRUE)
gwas <- gwas[, -1]
colnames(gwas)<-c("chr","BP","logp","r2")
gwas <- gwas %>%
  mutate(across(c("BP","logp","r2"), as.numeric))
chr=unique(gwas$chr)

gwas$col[gwas$r2<=0.2]<-0.2
gwas$col[gwas$r2>0.2 & gwas$r2<=0.4]<-0.4
gwas$col[gwas$r2>0.4 & gwas$r2<=0.6]<-0.6
gwas$col[gwas$r2>0.6 & gwas$r2<=0.8]<-0.8
gwas$col[gwas$r2>0.8 & gwas$r2<=1]<-1
cols<-c("0.2"="#113161","0.4"="#418BBB", "0.6"= "#1B733A","0.8"="#F6B063","1"="#D61C1F")  
gwas$target_snp<-target_snp

#format eqtl
eqtl<-read.table(eqtlfile,header = TRUE)
eqtl <- eqtl[, -1]
colnames(eqtl)<-c("chr","BP","logp","r2")
eqtl <- eqtl %>%
  mutate(across(c("BP","logp","r2"), as.numeric))

eqtl$col[eqtl$r2<=0.2]<-0.2
eqtl$col[eqtl$r2>0.2 & eqtl$r2<=0.4]<-0.4
eqtl$col[eqtl$r2>0.4 & eqtl$r2<=0.6]<-0.6
eqtl$col[eqtl$r2>0.6 & eqtl$r2<=0.8]<-0.8
eqtl$col[eqtl$r2>0.8 & eqtl$r2<=1]<-1
cols<-c("0.2"="#113161","0.4"="#418BBB", "0.6"= "#1B733A","0.8"="#F6B063","1"="#D61C1F") 
eqtl$target_snp<-target_snp

##get the xlim 
start=min(min(eqtl$BP),min(gwas$BP))
end=max(max(eqtl$BP),max(gwas$BP))

##get the overlap between eqtl and gwas.
overlap<-merge(gwas,eqtl,by="BP")

#############################################################################
#plot combined plots.
subset=gwas[gwas$BP==target_pos,]
gwas_plot<-ggplot(overlap,aes(BP/1000000,logp.x,color=factor(col.x)))+geom_point(shape=16,size=2)+
  geom_point(data=subset,aes(BP/1000000,logp),shape=16,size=4,color="purple")+
  geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "grey", size = 0.8)+
  xlab("")+ylab(expression(GWAS-log(italic(P))))+
  theme_classic()+
  ggtitle(trait)+
  theme(        legend.position = "none",
                plot.title=element_text(size=16,color="black"),
                axis.title.y=element_text(size=16,color="black"),
                axis.title.x=element_text(size=16,color="black"),
                axis.text.y=element_text(size=16,color="black",hjust=0.5,vjust=0.5),
                axis.text.x=element_blank())+
  scale_color_manual(values=cols)+
#    geom_text_repel(
#    data = eqtl[eqtl$BP==target_pos,],
#    aes(label = target_snp),
#    size = 3,
#    box.padding = unit(0.35, "lines"),
#    point.padding = unit(0.3, "lines")
#   )+
#   scale_y_continuous(breaks=c(0,1,2,3,4,5))+
  theme(plot.margin=unit(c(2,1.4,0.4,0.4),"cm"))#top right botton left

subset=eqtl[eqtl$BP==target_pos,]
xlabel=paste(chr,"position (Mb)")
eqtl_plot<-ggplot(overlap,aes(BP/1000000,logp.y,color=factor(col.y)))+geom_point(shape=16,size=2)+
  geom_point(data=subset,aes(BP/1000000,logp),shape=16,size=4,color="purple")+
  geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "grey", size = 0.8)+
  xlab(xlabel)+ylab(expression(eQTL-log(italic(P))))+
  theme_classic()+
  ggtitle(paste0(tissue,": ",gene_name))+
  theme(        legend.position = "none",
                plot.title=element_text(size=16,color="black"),
                axis.title.y=element_text(size=16,color="black"),
                axis.title.x=element_text(size=16,color="black"),
                axis.text.y=element_text(size=16,color="black",hjust=0.5,vjust=0.5),
                axis.text.x=element_text(size=16,color="black",hjust=0.5,vjust=0.1))+
  scale_color_manual(values=cols)+
  #geom_text_repel(
  #  data = eqtl[eqtl$BP==target_pos,],
  #  aes(label = target_snp),
  #  size = 3,
   # box.padding = unit(0.35, "lines"),
  #  point.padding = unit(0.3, "lines")
  #)+
  #scale_y_continuous(breaks=c(0,1,2,3,4,5))+
  theme(plot.margin=unit(c(0.2,1.4,0.4,0.4),"cm"))#top right botton left

subset=overlap[overlap$BP==target_pos,]
res=cor.test(overlap$logp.y, overlap$logp.x)
p=round(res$p.value,3)
r=round(res$estimate,2)

overlap_plot<-ggplot(overlap,aes(logp.y,logp.x))+geom_point(aes(color=factor(col.x)),shape=16,size=2.8)+ #x eqtl; y gwas
  geom_point(data=subset,aes(logp.y,logp.x),shape=16,size=4,color="purple")+
  geom_text_repel(
    data = subset,
    aes(label = target_snp.x),
    color = "black",
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = "black",
    segment.size = 0.5,
    nudge_y = -0.3,
    force = 10
  ) +
  xlab(expression(eQTL-log(italic(P))))+ylab(expression(GWAS-log(italic(P))))+

  theme_classic()+
  theme(        plot.title=element_text(size=16,color="black"),
                axis.title.y=element_text(size=16,color="black"),
                axis.title.x=element_text(size=16,color="black"),
                axis.text.y=element_text(size=16,color="black",hjust=0.5,vjust=0.5),
                axis.text.x=element_text(size=16,color="black",hjust=0.5,vjust=0.1),
                legend.position="none")+
  scale_color_manual(values=cols)+
  # guides(fill=guide_legend(override.aes=list(shape=15,col=cols,size=4), title=bquote(~r^2)))+
   # geom_smooth(method=lm,se=FALSE)+
  #geom_text_repel(
  #  data = overlap[overlap$BP==target_pos,],
  #  aes(label = target_snp.x),
  #  size = 3,
   # box.padding = unit(0.35, "lines"),
   # point.padding = unit(0.3, "lines")
  #)+
  theme(plot.margin=unit(c(2,3,0.4,0.4),"cm"))
  #top right botton left
legend_box = data.frame(x = 0.8, y = seq(0.4, 0.28, -0.03))

overlap_plot=ggdraw(overlap_plot)+
             geom_rect(data = legend_box,
                      aes(xmin = x, xmax = x + 0.03, ymin = y, ymax = y + 0.03),
                      color = "black",
                      fill = rev(c("#113161", "#418BBB", "#1B733A", "#F6B063", "#D61C1F"))) +
            draw_label("0.8", x = legend_box$x[1] + 0.03, y = legend_box$y[1], hjust = -0.3, size = 14) +
            draw_label("0.6", x = legend_box$x[2] + 0.03, y = legend_box$y[2], hjust = -0.3, size = 14) +
            draw_label("0.4", x = legend_box$x[3] + 0.03, y = legend_box$y[3], hjust = -0.3, size = 14) +
            draw_label("0.2", x = legend_box$x[4] + 0.03, y = legend_box$y[4], hjust = -0.3, size = 14) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.03, y = legend_box$y[1], vjust = -2, size = 14)

#####################################
#combine multiple plot into one plot.

library(ggpubr)

pdf(paste0(outdir,"/",trait, "_",tissue, "_",gene,".pdf"),width=10,height=6)

ggarrange(ggarrange(gwas_plot, eqtl_plot, nrow = 2, labels = c("A", "B"),font.label=list(color="black",face = "bold", size=18)),
           ggarrange(overlap_plot, labels = "C",font.label=list(color="black",face = "bold", size=18)), # Second row with box and dot plots
          ncol = 2                                      # Labels of the scatter plot
          ) 

dev.off()