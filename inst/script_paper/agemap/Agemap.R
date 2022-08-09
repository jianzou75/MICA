## run in server
library(parallel)
agemap_mcc_test = mscc.permute.test(agemap_data_list, agemap_label_list)


### after server
library(ggpubr)
library(org.Mm.eg.db)
gene_map <- select(org.Mm.eg.db, keys=colnames(agemap_data_list[[1]]),
                   columns=c("GENENAME","SYMBOL"), keytype="UNIGENE")

## 01. two genes compare
gene1 = scale(sapply(agemap_data_list, function(x) x[,"Mm.277792"])) %>% data.frame()  %>% mutate(time =  agemap_label_list[[1]])
gene2 = scale(sapply(agemap_data_list, function(x) x[,"Mm.27897"])) %>% data.frame()  %>% mutate(time =  agemap_label_list[[1]])
colnames(gene1)[1:5] = colnames(gene2)[1:5] = c("Ms", "Ad", "Hi", "Sp", "S")

gene1.long = gene1 %>% tibble::rownames_to_column("ID") %>% tidyr::gather(., tissue, expression, -c(ID, time), factor_key=TRUE) %>% mutate(gene = "Mm.277792")
gene2.long = gene2 %>% tibble::rownames_to_column("ID") %>% tidyr::gather(., tissue, expression, -c(ID, time), factor_key=TRUE) %>% mutate(gene = "Mm.27897")
gene.cmb = rbind(gene1.long, gene2.long)
gene.cmb$tissue = as.character(gene.cmb$tissue)
gene.cmb$tissue[gene.cmb$tissue == "Ms"] = "Muscle"
gene.cmb$tissue[gene.cmb$tissue == "Ad"] = "Adrenal Glands"
gene.cmb$tissue[gene.cmb$tissue == "Hi"] = "Hippocampus"
gene.cmb$tissue[gene.cmb$tissue == "Sp"] = "Spinal Cord"
gene.cmb$tissue[gene.cmb$tissue == "S"]  = "Spleen"
gene.cmb$gene[gene.cmb$gene == "Mm.277792"] = "Col1a2"
gene.cmb$gene[gene.cmb$gene == "Mm.27897"]  = "Dnaja1"

fig1 = ggboxplot(gene.cmb, x = "time", y = "expression",
          add = "jitter", add.params = list(alpha = 1.2, size = 1),
          facet.by = c("gene", "tissue"), color = "tissue",
          font.label = list(face = "bold")) +
  theme_bw() +
  theme(legend.position="none", strip.text.x = element_text(face = "bold", size = 12),
        strip.text.y = element_text(face = "bold", size = 12),
        axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(size=12)) +
  ylim(c(-2,2)) +
  stat_summary(fun = median, geom = "line", aes(group=1, color=tissue), size = 1) +
  ylab("Expression") + xlab("Age")
ggsave(fig1, width = 15, height = 6, unit = "in", file = "figure1_Col1a2_Dnaja1_pattern.svg")

