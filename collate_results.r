library(dplyr)
library(here)
library(tidyr)
library(stringr)

load(here("data", "all.rdata"))
coloc_cis <- readRDS(here("data", "cis_coloc_res.rds"))
res_cis

# Primary causal results - cis MR + coloc
dim(res_cis)



# all MR results
head(res_cis)

res <- bind_rows(
    subset(res_cis %>% mutate(cistrans="cis"), select=c(exposure, outcome, cistrans, method, nsnp, b, se, pval)),
    subset(res_trans %>% mutate(cistrans="trans"), select=c(exposure, outcome, cistrans, method, nsnp, b, se, pval))
) %>% 
    arrange(pval)

head(res$outcome)
temp <- do.call(rbind, str_split(res$outcome, " \\|\\| "))
res$outcome <- temp[,1]
head(temp)
head(res)

head(res)
write.csv(res, file=here("data", "all_mr_results.csv"))

head(coloc_cis)
table(coloc_cis$outcome)

coloc_res <- coloc_cis %>% select(exposure, outcome, target_variant=snp, coloc_h4=PP.H4.abf, coloc_nsnp=nsnps) %>% left_join(
    ., select(res, exposure, outcome, mr_beta = b, mr_se = se, mr_pval = pval)
) %>% 
    arrange(mr_pval) %>% 
    filter(!duplicated(paste(exposure, outcome, target_variant))) %>%
    arrange(desc(coloc_h4))
coloc_res

write.csv(coloc_res, file=here("data", "cis_coloc_results.csv"))
