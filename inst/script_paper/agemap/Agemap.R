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

figs4 = ggboxplot(gene.cmb, x = "time", y = "expression",
          add = "jitter", add.params = list(alpha = 1.2, size = 1),
          facet.by = c("gene", "tissue"), color = "gene",
          font.label = list(face = "bold", size = 20)) +
  theme_bw() +
  theme(legend.position="none", strip.text.x = element_text(face = "bold", size = 20),
        strip.text.y = element_text(face = "bold", size = 20),
        axis.title = element_text(size=20,face="bold"),
        axis.text = element_text(size=20)) +
  ylim(c(-2,2)) +
  stat_summary(fun = median, geom = "line", aes(group=1, color=gene), size = 1) +
  ylab("Expression") + xlab("Age") +
  scale_color_manual(values = c("#e15f41", "#546de5", "#3dc1d3"))
ggsave(figs4, width = 15, height = 6, unit = "in", file = "../../202203_multi_class_correlation/manuscript_v4/figures_raw/figures4_Col1a2_Dnaja1_pattern.pdf")


## 02. MC-TC identified genes
tissue_gene_boxplot <- function(tissue_label, gene, tissue_name){

  tissue1 = scale(agemap_data_list[[tissue_label[1]]][,gene]) %>% data.frame()  %>% mutate(time =  agemap_label_list[[tissue_label[1]]])
  tissue2 = scale(agemap_data_list[[tissue_label[2]]][,gene]) %>% data.frame()  %>% mutate(time =  agemap_label_list[[tissue_label[2]]])

  tissue1.long = tissue1 %>% tibble::rownames_to_column("ID") %>% tidyr::gather(., gene, expression, -c(ID, time), factor_key=TRUE) %>% mutate(tissue = tissue_name[1])
  tissue2.long = tissue2 %>% tibble::rownames_to_column("ID") %>% tidyr::gather(., gene, expression, -c(ID, time), factor_key=TRUE) %>% mutate(tissue = tissue_name[2])
  tissue = rbind(tissue1.long, tissue2.long)
 # tissue$gene = gene_map$SYMBOL[match(tissue$gene,gene_map$UNIGENE)]

  figure <- ggboxplot(tissue, x = "time", y = "expression",
                      add = "jitter", add.params = list(alpha = 0.6, size = 1),
                      facet.by = c("tissue", "gene"), color = "gene",
                      font.label = list(face = "bold", size = 25)) +
    theme_bw() +
    theme(legend.position="none", strip.text.x = element_text(face = "bold", size = 15),
          strip.text.y = element_text(face = "bold", size = 15),
          axis.title = element_text(size=15,face="bold"),
          axis.text = element_text(size=15)) +
    ylim(c(-2,2)) +
    stat_summary(fun = median, geom = "line", aes(group=1, color=gene), size = 1) +
    ylab("") + xlab("")

  return(figure)
}

grp_name = combn(c("Ms", "Ad", "Hi", "Sp", "S"), 2)
grp_num  = combn(1:5, 2)
bp_list <- c()

for (i in 1:10){
  tissue_num  = grp_num[,i]
  tissue_name = grp_name[,i]
  gene_identify = which(agemap_result$MCTC_tbl$qval < 0.05 & agemap_result$pairs_MI_tbl_list[[i]]$qval < 0.05)
  bp_list <-c(bp_list, list(tissue_gene_boxplot(tissue_num, gene_identify, tissue_name)))
}

bp1 = plot_grid(bp_list[[1]], NULL, nrow = 1, rel_widths = c(13, 0))
bp2 = plot_grid(bp_list[[2]], NULL, nrow = 1, rel_widths = c(5, 8))
bp4 = plot_grid(bp_list[[4]], NULL, nrow = 1, rel_widths = c(5, 8))
bp8 = plot_grid(bp_list[[8]], NULL, nrow = 1, rel_widths = c(9, 4))
bp10 = plot_grid(bp_list[[10]], NULL, nrow = 1, rel_widths = c(5, 8))

supfig1A = plot_grid(bp1, bp2, bp4, bp8, bp10,
                    ncol = 1, align = "v", axis = "tblr")
supfig1B = plot_grid(bp5, bp6, bp7, bp8, bp10,
                     ncol = 1, align = "v", axis = "tblr")


ggsave(supfig1A, width = 10, height = 20, unit = "in", file = "fig3_mcmi.pdf")
ggsave(supfig1B, width = 20, height = 25, unit = "in", file = "supfig1B_mcmi.pdf")
