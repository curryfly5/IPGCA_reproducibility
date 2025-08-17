library(SingleCellExperiment)
suppressWarnings(library(BayesPrism))
library(biomaRt)
library(arrow)

setwd('/disk212/yupf/database/scRNA-seq/NewAtlas/figures')
# load("/disk212/yupf/database/scRNA-seq/NewAtlas/decon/JHtest/counts_test/sc_bk_counts.RData")

bk.data=read.table('/disk191_2/yupf/GUT_RNA/03_EXPRESSION/counts.txt',row.names = 1,header = TRUE)
bk.dat=t(bk.data)
sc.sce = readRDS('/disk222/yupf/RNA-seq/eQTL/deconvolution/atlas_sub_sce.rds')
cell.type.labels=as.vector(colData(sc.sce)$level3)
sc.dat=read_parquet('/disk222/yupf/RNA-seq/eQTL/deconvolution/altas_counts_sub.parquet')
sc.dat <- as.data.frame(sc.dat)
rownames(sc.dat)<-sc.dat$index
sc.dat <- sc.dat[, -1]
sc.dat=t(sc.dat)

dim(sc.dat)
#> [1] 117821  19123
sc.dat <- sc.dat[!rownames(sc.dat) %in% "__index_level_0__", ]
dim(sc.dat)
#> [1] 117820  19117

par(mar = c(2, 2, 2, 2))
pdf("cor_ct.pdf", width=18, height=18)
#细胞类型和状态标记的质量控制
plot.cor.phi (input=sc.dat, 
              input.labels=cell.type.labels, 
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=2, cexCol=2
)
dev.off()
# 过滤异常基因
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)


bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)

dim(sc.dat.filtered)
#> 117820  17876

# 在反卷积前，全部转化为EMSEMBLE IDs
my_mart<- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(my_mart)#查看有什么数据库
pig_dataset<- useDataset("sscrofa_gene_ensembl",mart = my_mart) 
pig_dataset@attributes$name[1:20]
attr1<- c("ensembl_gene_id","chromosome_name","external_gene_name","gene_biotype","description")
sc_gene=colnames(sc.dat)
count_ID<- getBM(attributes = attr1,
                 filters = "external_gene_name",
                 values = sc_gene, 
                 mart = pig_dataset)

sc_gene_id <- lapply(sc_gene, function(x) {
  if (x %in% count_ID$external_gene_name) {
    return(count_ID$ensembl_gene_id[count_ID$external_gene_name == x])
  } else {
    return(x)
  }
})

colnames(sc.dat) <- sc_gene_id


#检查不同基因表达的一致性 基因的类型——选择protein coding较合适
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")
# Select marker genes (Optional)
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)

sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)

# 构造一个棱镜对象
myPrism <- new.prism(
  reference=sc.dat.filtered, 
  mixture=bk.dat,
  input.type="count.matrix", 
  cell.type.labels = cell.type.labels, 
  cell.state.labels = cell.type.labels,
  key=NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
bp.res <- run.prism(prism = myPrism, n.cores=50)

### 获取细胞比例
theta <- get.fraction (bp=bp.res,
                       which.theta="final", #first
                       state.or.type="type")


### 获取某个细胞的表达量
Z.Enterocytes <- get.exp (bp=bp.res,
                          state.or.type="type",
                          cell.name="Enterocytes")
ct=levels(colData(sc.sce)$level3)
for (i in ct){
  Z.ct <- get.exp (bp=bp.res,
                            state.or.type="type",
                            cell.name=i)
  write.table(Z.ct,file=paste0('/disk222/yupf/RNA-seq/eQTL/deconvolution/celltype_exp/','exp_',i,'_decon.txt'),sep='\t')
  print(i)
}

### 几种变换方式

scale_theta <- as.data.frame(scale(theta))
write.table(scale_theta, file = "/disk222/yupf/RNA-seq/eQTL/deconvolution/scale_percent.txt", sep = "\t")
log_theta <- log(theta)
write.table(log_theta, file = "/disk222/yupf/RNA-seq/eQTL/deconvolution/log_percent.txt", sep = "\t")

inverse_normal_transform = function(x) {
  qnorm(rank(x) / (length(x)+1))
}

theta_inv= t(apply(theta, MARGIN = 1, FUN = inverse_normal_transform))
write.table(theta_inv, file = "/disk222/yupf/RNA-seq/eQTL/deconvolution/inv_percent.txt", sep = "\t")
