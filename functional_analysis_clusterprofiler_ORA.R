##ref: https://rpubs.com/jrgonzalezISGlobal/enrichment
##ORA using clusterprofiler
library(readxl)
library(org.Hs.eg.db)
library(clusterProfiler)

genes <- read_excel("H:/genes.xlsx")
down_genes<-genes[-c(582,583)]
up_genes <- genes[c(582,583)]
gene_id <- mapIds(org.Hs.eg.db, down_genes, 'ENTREZID', 'SYMBOL')
gene_id_up <- mapIds(org.Hs.eg.db, up_genes, 'ENTREZID', 'SYMBOL')
##GO-BP_down
go_bp <- enrichGO(gene = gene_id, 
                   ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   keyType = "ENTREZID",
                   readable=TRUE,
                   pvalueCutoff = 0.05)
tab_go_bp <- as.data.frame(go_bp)
write.csv(tab_go_bp,"H:/go_bp.csv")
##GO-BP_up
go_bp_up <- enrichGO(gene = gene_id_up, 
                  ont = "BP",
                  OrgDb ="org.Hs.eg.db",
                  keyType = "ENTREZID",
                  readable=TRUE,
                  pvalueCutoff = 0.05)
tab_go_bp_up <- as.data.frame(go_bp_up)
write.csv(tab_go_bp,"H:/go_bp_up.csv")

##Dotplot
#ref:https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html, https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
png("H:/dotplot_GO_BP.png",height = 1500, width = 1500)
dotplot_down<-dotplot(go_bp, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("GO:Biological Process Downregulated") +
  theme(text = element_text(size = 17))
dotplot_up<-dotplot(go_bp_up, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("GO:Biological Process Upregulated")+
  theme(text = element_text(size = 17))
cowplot::plot_grid(dotplot_down, dotplot_up, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, .8))
dev.off()

##GO_MF_down
go_mf <- enrichGO(gene = gene_id, 
                  ont = "MF",
                  OrgDb ="org.Hs.eg.db",
                  keyType = "ENTREZID",
                  readable=TRUE,
                  pvalueCutoff = 0.05)
tab_go_mf <- as.data.frame(go_mf)
write.csv(tab_go_mf,"H:/go_mf.csv")
#GO_MF_up
go_mf_up <- enrichGO(gene = gene_id_up, 
                  ont = "MF",
                  OrgDb ="org.Hs.eg.db",
                  keyType = "ENTREZID",
                  readable=TRUE,
                  pvalueCutoff = 0.05)
tab_go_mf_up <- as.data.frame(go_mf_up)
write.csv(tab_go_mf_up,"H:/go_mf_up.csv")

##Dotplot
#ref:https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html, https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
png("H:/dotplot_GO_MF.png",height = 1500, width = 1500)
dotplot_down<-dotplot(go_mf, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("GO:Molecular Function Downregulated") +
  theme(text = element_text(size = 17))
dotplot_up<-dotplot(go_mf_up, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("GO:Molecular Function Upregulated")+
  theme(text = element_text(size = 17))
cowplot::plot_grid(dotplot_down, dotplot_up, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, .8))
dev.off()

##Reactome Pathway
#ref: http://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
library(ReactomePA)
reactome_pathway <- enrichPathway(gene=gene_id, 
                                  pvalueCutoff = 0.05, 
                                  readable=TRUE,
                                  organism = "human")
tab_reactome_pathway <- as.data.frame(reactome_pathway)
write.csv(tab_reactome_pathway,"H:/reactome_pathway.csv")

##Reactome Pathway Upregulated
reactome_pathway_up <- enrichPathway(gene=gene_id_up, 
                                  pvalueCutoff = 0.05, 
                                  readable=TRUE,
                                  organism = "human")
tab_reactome_pathway_up <- as.data.frame(reactome_pathway_up)
write.csv(tab_reactome_pathway,"H:/reactome_pathway.csv")
##Dotplot
#ref:https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html, https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
png("H:/dotplot_reactome_pathway.png",height = 2200, width = 2200)
dotplot_down<-dotplot(reactome_pathway, showCategory=30,orderBy = "p.adjust",font.size=22) + ggtitle("Reactome Pathway - Downregulated") +
  theme(text = element_text(size = 22))
dotplot_up<-dotplot(reactome_pathway_up, showCategory=30,orderBy = "p.adjust",font.size=22) + ggtitle("Reactome Pathway - Upregulated")+
  theme(text = element_text(size = 22))
cowplot::plot_grid(dotplot_down, dotplot_up, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, .8))
dev.off()


##disease from DOSE package
#ref: http://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html
library(DOSE)
disease_DO <- enrichDO(gene_id,
                       ont = "DO",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

disease_DOSE <- as.data.frame(disease_DO)
write.csv(disease_DOSE,"H:/disease_DOSE.csv")
##Dotplot
#ref:https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html, https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
png("H:/dotplot_disease_DO.png",height = 1281, width = 961)
dotplot(disease_DO, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("Disease Ontology - Downregulated") +
  theme(text = element_text(size = 17))+
  theme(plot.title = element_text(hjust = 0.5))

dev.off()



##disease using enrich DisGeNET
#ref:http://yulab-smu.top/biomedical-knowledge-mining-book/dose-enrichment.html
disease_DGN <-enrichDGN(gene_id,
                        pvalueCutoff = 0.05,
                        readable=TRUE)
                      
disease_DGN_db <- as.data.frame(disease_DGN)
write.csv(disease_DGN_db,"H:/disease_DGN_db.csv")

##disease using enrich DisGeNET for upregulated
disease_DGN_up <-enrichDGN(gene_id_up,
                        pvalueCutoff = 0.05,
                        readable=TRUE)

disease_DGN_db_up <- as.data.frame(disease_DGN_up)
write.csv(disease_DGN_db_up,"H:/disease_DGN_db_up.csv")
##Dotplot
#ref:https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html, https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
png("H:/dotplot_disease.png",height = 1500, width = 1500)
dotplot_down<-dotplot(disease_DGN, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("Disease Association : DisGeNET - Downregulated") +
  theme(text = element_text(size = 17))
dotplot_up<-dotplot(disease_DGN_up, showCategory=30,orderBy = "p.adjust",font.size=17) + ggtitle("Disease Association : DisGeNET - Upregulated")+
  theme(text = element_text(size = 17))
cowplot::plot_grid(dotplot_down, dotplot_up, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, .8))
dev.off()


