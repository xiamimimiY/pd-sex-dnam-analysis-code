rm(genes_uniquePD_)
cmPDann <- read.csv("cmPDann.csv")  
cmPDann1 <- cmPDann[!is.na(cmPDann$F_to_M.gene) & nchar(cmPDann$F_to_M.gene) > 0, ]
uniquePDann1 <- uniquePDann[!is.na(uniquePDann$F_to_M.gene) & nchar(uniquePDann$F_to_M.gene) > 0, ]

# Load required packages
library(clusterProfiler)
?clusterProfiler
?enrichGO
?enrichKEGG


# Read gene list from a file (assuming one gene per line)
genes_cmPD <- cmPDann1$F_to_M.gene
genes_uniquePD <- uniquePDann1$F_to_M.gene
class(cmPD_genes)
write.table(genes_cmPD, file = "genes_cmPD.csv",sep = ",", row.names=FALSE,col.names=TRUE)
write.table(genes_uniquePD, file = "genes_uniquePD.csv",sep = ",", row.names=FALSE,col.names=TRUE)
write.table(genes_uniquePD_down, file = "genes_uniquePD_down.csv",sep = ",", row.names=FALSE,col.names=TRUE)
write.table(genes_uniquePD_up, file = "genes_uniquePD_up.csv",sep = ",", row.names=FALSE,col.names=TRUE)
uniquePD_Promoter <- uniquePDann1[uniquePDann1$F_to_M.feature %in% c("TSS200", "TSS1500", "5'UTR", "1stExon"), ] #978
AA <- uniquePDann1[uniquePDann1$F_to_M.feature %in% c("Body", "3'UTR"), ]  #615
rm(AA)
genes_uniquePD_Promoter <- uniquePD_Promoter$F_to_M.gene
write.table(genes_uniquePD_Promoter, file = "genes_uniquePD_Promoter.csv",sep = ",", row.names=FALSE,col.names=TRUE)


# Perform Gene Ontology analysis using enrichGO() function,
Go_cmPD <- enrichGO(gene         = genes_cmPD,       # input gene list
                    OrgDb    = "org.Hs.eg.db", # genome background
                    ont         = "ALL",        # ontology to test (Biological Process)
                    pAdjustMethod = "BH",      # multiple testing adjustment method
                    qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                    keyType = "SYMBOL")
print(Go_cmPD)
Go_cmPD = as.data.frame(Go_cmPD)

Go_uniquePD <- enrichGO(gene         = genes_uniquePD,       # input gene list
                    OrgDb    = "org.Hs.eg.db", # genome background
                    ont         = "ALL",        # ontology to test (Biological Process)
                    pAdjustMethod = "BH",      # multiple testing adjustment method
                    qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                    keyType = "SYMBOL")  
print(Go_uniquePD)
write.table(Go_cmPD, 'Go_cmPD.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(Go_uniquePD, 'Go_uniquePD.txt', sep = '\t', row.names = FALSE, quote = FALSE)
Go_uniquePD = as.data.frame(Go_uniquePD)

#————————————————————————————————————————————————————————————————————————————————————————————
class(uniquePDann1)
uniquePDann1_up <- subset(uniquePDann1, Direction == "++")
uniquePDann1_down <- subset(uniquePDann1, Direction == "--")
genes_uniquePD_up <- uniquePDann1_up$F_to_M.gene
genes_uniquePD_down <- uniquePDann1_down$F_to_M.gene
# Perform Gene Ontology analysis using enrichGO() function,
Go_uniquePD_up <- enrichGO(gene         = genes_uniquePD_up,       # input gene list
                    OrgDb    = "org.Hs.eg.db", # genome background
                    ont         = "ALL",        # ontology to test (Biological Process)
                    pAdjustMethod = "BH",      # multiple testing adjustment method
                    qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                    keyType = "SYMBOL")
print(Go_uniquePD_up)

Go_uniquePD_down <- enrichGO(gene         = genes_uniquePD_down,       # input gene list
                        OrgDb    = "org.Hs.eg.db", # genome background
                        ont         = "ALL",        # ontology to test (Biological Process)
                        pAdjustMethod = "BH",      # multiple testing adjustment method
                        qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                        keyType = "SYMBOL")  
print(Go_uniquePD_down)
write.table(Go_uniquePD_up, 'Go_uniquePD_up.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(Go_uniquePD_down, 'Go_uniquePD_down.txt', sep = '\t', row.names = FALSE, quote = FALSE)
Go_uniquePD_up = as.data.frame(Go_uniquePD_up)
Go_uniquePD_down = as.data.frame(Go_uniquePD_down)

Go_uniquePD_Promoter <- enrichGO(gene         = genes_uniquePD_Promoter,       # input gene list
                           OrgDb    = "org.Hs.eg.db", # genome background
                           ont         = "ALL",        # ontology to test (Biological Process)
                           pAdjustMethod = "BH",      # multiple testing adjustment method
                           qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                           keyType = "SYMBOL")
print(Go_uniquePD_Promoter)
write.table(Go_uniquePD_Promoter, 'Go_uniquePD_Promoter.txt', sep = '\t', row.names = FALSE, quote = FALSE)
Go_uniquePD_Promoter = as.data.frame(Go_uniquePD_Promoter)

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————
uniquePD2.4 = read.csv("uniquePD2.4.csv")
rownames(uniquePD2.4)=uniquePD2.4[,1]
colnames(uniquePD2.4)[1] = 'CpGs'
uniquePDann2.4 = merge(uniquePD2.4,PDDMP,by="CpGs")
uniquePDann2.4 <- uniquePDann2.4[!is.na(uniquePDann2.4$F_to_M.gene) & nchar(uniquePDann2.4$F_to_M.gene) > 0, ]
genes_uniquePD2.4 <- uniquePDann2.4$F_to_M.gene
# Perform Gene Ontology analysis using enrichGO() function,
Go_uniquePD2.4 <- enrichGO(gene         = genes_uniquePD2.4,       # input gene list
                        OrgDb    = "org.Hs.eg.db", # genome background
                        ont         = "ALL",        # ontology to test (Biological Process)
                        pAdjustMethod = "BH",      # multiple testing adjustment method
                        qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                        keyType = "SYMBOL")  
print(Go_uniquePD2.4)
write.table(Go_uniquePD2.4, 'Go_uniquePD2.4.txt', sep = '\t', row.names = FALSE, quote = FALSE)
Go_uniquePD2.4 = as.data.frame(Go_uniquePD2.4)

uniquePDann2.4_up <- subset(uniquePDann2.4, Direction == "++")
uniquePDann2.4_down <- subset(uniquePDann2.4, Direction == "--")
genes_uniquePD2.4_up <- uniquePDann2.4_up$F_to_M.gene
genes_uniquePD2.4_down <- uniquePDann2.4_down$F_to_M.gene
# Perform Gene Ontology analysis using enrichGO() function,
Go_uniquePD2.4_up <- enrichGO(gene         = genes_uniquePD2.4_up,       # input gene list
                           OrgDb    = "org.Hs.eg.db", # genome background
                           ont         = "ALL",        # ontology to test (Biological Process)
                           pAdjustMethod = "BH",      # multiple testing adjustment method
                           qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                           keyType = "SYMBOL")
print(Go_uniquePD2.4_up)

Go_uniquePD2.4_down <- enrichGO(gene         = genes_uniquePD2.4_down,       # input gene list
                             OrgDb    = "org.Hs.eg.db", # genome background
                             ont         = "ALL",        # ontology to test (Biological Process)
                             pAdjustMethod = "BH",      # multiple testing adjustment method
                             qvalueCutoff = 0.05,      # q-value cutoff for significant terms
                             keyType = "SYMBOL")  
print(Go_uniquePD2.4_down)
write.table(Go_uniquePD2.4_up, 'Go_uniquePD2.4_up.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(Go_uniquePD2.4_down, 'Go_uniquePD2.4_down.txt', sep = '\t', row.names = FALSE, quote = FALSE)
Go_uniquePD2.4_up = as.data.frame(Go_uniquePD2.4_up)
Go_uniquePD2.4_down = as.data.frame(Go_uniquePD2.4_down)














enrichGO(
  gene,
  OrgDb,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = FALSE
)

# If GitHub installation requires credentials, configure them locally and never commit tokens.
BiocManager::install("KEGG.db")
devtools::install_local("KEGG.db_3.2.4.tar.gz")
library(clusterProfiler)
library(KEGG.db)
packageVersion("clusterProfiler")


?enrichKEGG
library(org.Hs.eg.db)
entrezid_cmPD <- mapIds(org.Hs.eg.db, keys = genes_cmPD, column = "ENTREZID", keytype = "SYMBOL")
entrezid_uniquePD <- mapIds(org.Hs.eg.db, keys = genes_uniquePD, column = "ENTREZID", keytype = "SYMBOL")
class(entrezid_cmPD)
head(entrezid_cmPD)
entrezid_cmPD <- na.omit(entrezid_cmPD, subset = c("entrezid_cmPD"))
entrezid_cmPD <- as.data.frame(entrezid_cmPD)
entrezid_cmPD <- as.numeric(entrezid_cmPD)
entrezid_cmPD <- as.vector(as.numeric(entrezid_cmPD)) #509

entrezid_uniquePD <- na.omit(entrezid_uniquePD, subset = c("entrezid_uniquePD"))
entrezid_uniquePD <- as.data.frame(entrezid_uniquePD)
entrezid_uniquePD <- as.numeric(entrezid_uniquePD)
entrezid_uniquePD <- as.vector(as.numeric(entrezid_uniquePD))  #1372

# Perform KEGG analysis using the reference genome for human (org.Hs.eg.db)
KEGG_cmPD <- enrichKEGG(gene = entrezid_cmPD$entrezid_cmPD,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5,
                        pAdjustMethod = "BH",
                        use_internal_data = T)
dim(KEGG_cmPD) #[1] 20  9
head(KEGG_cmPD@result,10)

KEGG_uniquePD <- enrichKEGG(gene = entrezid_uniquePD$entrezid_uniquePD,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.5,
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.5,
                        use_internal_data = T)
dim(KEGG_uniquePD) #[1] 23  9
# View the top 10 enriched pathways
head(KEGG_uniquePD@result, 10)

write.table(KEGG_cmPD, 'KEGG_cmPD.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(KEGG_uniquePD, 'KEGG_uniquePD.txt', sep = '\t', row.names = FALSE, quote = FALSE)
















#——————————————————————————————————————————————————————————————————————————————
cmPD_genelist <- read_excel("cmPD_genelist.xlsx")
uniquePD_genelist <- read_excel("uniquePD_genelist.xlsx")
write.table(cmPD_genes, file = "cmPD_genes.txt", sep = "", row.names = FALSE)
cmPD_genes <- read.delim('cmPD_genes.txt', header = TRUE, stringsAsFactors = FALSE)[[1]]
org <- "org.Hs.eg.db"
enrich_result <- enrichGO(gene = cmPD_genelist$Gene, universe = org, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, keyType = "SYMBOL", OrgDb = org)


###__________________________________________________________________________________________________
## KEGG pathway analysis
run_kegg <- function(gene_up,gene_down,geneList=F,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.up.csv'))
  
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.down.csv'))
  
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  kk=kk.diff
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.diff.csv'))
  
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = paste0(pro,'_kegg_up_down.png') )
  
if(geneList){
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
  kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  tmp=kk@result
  write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
  
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.png'))
  # 
}
  
}

### GO database analysis 
run_go <- function(gene_up,gene_down,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  g_list=list(gene_up=gene_up,
              gene_down=gene_down,
              gene_diff=gene_diff)
  if(T){
    go_enrich_results <- lapply( g_list , function(gene) {
      lapply( c('BP','MF','CC') , function(ont) {
        cat(paste('Now process ',ont ))
        ego <- enrichGO(gene          = gene,
                        #universe      = gene_all,
                        OrgDb         = org.Hs.eg.db,
                        ont           = ont ,
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.99,
                        qvalueCutoff  = 0.99,
                        readable      = TRUE)
        
        print( head(ego) )
        return(ego)
      })
    })
    save(go_enrich_results,file =paste0(pro, '_go_enrich_results.Rdata'))
    
  }
  load(file=paste0(pro, '_go_enrich_results.Rdata'))
  
  n1= c('gene_up','gene_down','gene_diff')
  n2= c('BP','MF','CC') 
  for (i in 1:3){
    for (j in 1:3){
      fn=paste0(pro, '_dotplot_',n1[i],'_',n2[j],'.png')
      cat(paste0(fn,'\n'))
      png(fn,res=150,width = 1080)
      print( dotplot(go_enrich_results[[i]][[j]] ))
      dev.off()
    }
  }
  
  
}

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}


