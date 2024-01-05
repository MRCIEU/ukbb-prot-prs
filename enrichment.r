library(here)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(parallel)

# load mr results
load(here("data", "all.rdata"))
options(mc.cores=30)

# PRS universe
a <- read.csv(here("data", "Bristol_and_Biogen_PRS_Associations", "PRScs_Olink_Associaton_Biogen_AD.csv"))
a <- a %>% tidyr::separate(UKBPPP_ProteinID, sep=":", into=c("symbol", "uniprot", "v1", "v2"))
prs_universe <- unique(a$symbol)

enr_cis <- mclapply(unique(res_cis$id.outcome), \(id) {
    message(id)
    x <- subset(res_cis, id.outcome==id)
    f <- subset(x, pval < 0.05)$id.exposure
    if(length(f) == 0) return(NULL)
    ego <- enrichGO(gene          = f,
                    keyType       = "SYMBOL",
                    universe      = unique(x$id.exposure),
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "none",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
            readable = TRUE)
    ego <- as.data.frame(ego)
    ego$trait <- id
    return(ego)
}) %>% bind_rows()

enr_trans <- mclapply(unique(res_trans$id.outcome), \(id) {
    message(id)
    x <- subset(res_trans, id.outcome==id)
    f <- subset(x, pval < 0.05)$id.exposure
    if(length(f) == 0) return(NULL)
    ego <- enrichGO(gene          = f,
                    keyType       = "SYMBOL",
                    universe      = unique(x$id.exposure),
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "none",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
            readable = TRUE)
    ego <- as.data.frame(ego)
    ego$trait <- id
    return(ego)
}) %>% bind_rows()

enr_prs <- mclapply(unique(prs_pairs$opengwasid), \(id) {
    message(id)
    f <- subset(prs_pairs, opengwasid==id)$prot
    if(length(f) == 0) return(NULL)
    ego <- enrichGO(gene          = f,
                    keyType       = "SYMBOL",
                    universe      = prs_universe,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "ALL",
                    pAdjustMethod = "none",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
            readable = TRUE)
    ego <- as.data.frame(ego)
    ego$trait <- id
    return(ego)
}) %>% bind_rows()

enr_prs <- bind_rows(enr_prs)
enr_cis <- bind_rows(enr_cis)
enr_trans <- bind_rows(enr_trans)

save(enr_prs, enr_cis, enr_trans, file=here("data", "gsea_enrichments.rdata"))
