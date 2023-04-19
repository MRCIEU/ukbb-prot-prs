library(readxl)
library(dplyr)
library(ieugwasr)
library(tidyr)
library(TwoSampleMR)
library(ggplot2)
library(here)

download.file("https://www.biorxiv.org/content/biorxiv/early/2022/06/18/2022.06.17.496443/DC2/embed/media-2.xlsx?download=true", here("data", "media-2.xlsx"))
a <- read_xlsx(here("data", "media-2.xlsx"), sheet="ST6", skip=2)

b <- excel_sheets(here("data", "Supplementary_Table1.xlsx")) %>% lapply(., function(s) read_xlsx(here("data", "Supplementary_Table1.xlsx"), sheet=s) %>% mutate(trait=s, P_value=as.numeric(P_value))) %>% bind_rows()

b1 <- read_xlsx(here("data", "supplementary_tables.xlsx"), sheet=1, skip=1)


inst <- a %>% 
  filter(`Assay Target` %in% b$HGNC.symbol, rsID != "-") %>%
  dplyr::select(
    exposure="Assay Target",
    id.exposure="Assay Target",
    SNP="rsID", 
    beta.exposure="BETA (discovery, wrt. A1)", 
    se.exposure="SE (discovery)", 
    pval.exposure="log10(p) (discovery)", 
    eaf.exposure="A1FREQ (discovery)", 
    snpid="Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)",
    cistrans.exposure=`cis/trans`) %>% 
  mutate(pval.exposure = 10^-pval.exposure) %>%
  tidyr::separate(snpid, sep=":", into=c("chr.exposure","pos.exposure","other_allele.exposure", "effect_allele.exposure", "imp", "v"))

traits <- read_xlsx(here("data", "Supplementary_Table5.xlsx"))
b <- left_join(b, traits, by=c("trait"="code"))
b1 <- left_join(b1, traits, by=c("cancer"="code"))
prs_pairs <- bind_rows(
    b %>%
        dplyr::select(prot=HGNC.symbol, opengwasid, beta=Estimate, pval=P_value) %>%
        filter(!duplicated(paste(prot, opengwasid))),
    b1 %>%
        dplyr::select(prot=Assay.Target, opengwasid, beta=Estimate, se=Std_Error, pval=P) %>%
        filter(!duplicated(paste(prot, opengwasid)))
)
prs_pairs


lookups <- inner_join(inst, prs_pairs, by=c("exposure" = "prot"))

newids <- c("ieu-a-1121", "ieu-a-1127", "ieu-a-1128", "ieu-a-1228", "ieu-a-822", "ukb-a-60", "ukb-b-1316", "ukb-b-8837")

outcome_dat_all <- lapply(unique(traits$opengwasid), function(id)
{
    extract_outcome_data(unique(inst$SNP), id)
})

o <- bind_rows(outcome_dat_all)
dim(o)

dat <- harmonise_data(inst, o, action=1)
str(dat)
length(unique(dat$exposure))


eo_pairs <- unique(paste(prs_pairs$prot, prs_pairs$opengwasid))
length(eo_pairs)
dat_prs <- subset(dat, paste(id.exposure, id.outcome) %in% eo_pairs)
dim(dat_prs)
length(unique(paste(dat$id.exposure, dat$id.outcome)))
dat <- subset(dat, !duplicated(paste(dat$SNP, dat$id.exposure, dat$id.outcome)))


# Steiger filtering
sd_ids <- c("ieu-b-38", "ieu-b-39", "ieu-a-299", "ieu-a-300", "ieu-a-301", "ieu-a-302")
logor_ids <- traits$opengwasid[! traits$opengwasid %in% sd_ids]
dat$unit.exposure <- "SD"
dat$samplesize.exposure <- 35571
dat$unit.outcome <- "SD"
dat$unit.outcome[dat$id %in% logor_ids] <- "logOR"
dat <- add_rsq(dat)
dat$rsq.exposure <- get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)^2
dat$effective_n.outcome <- dat$samplesize.outcome
dat$effective_n.outcome[dat$unit == "logOR"] <- effective_n(dat$ncase.outcome[dat$unit == "logOR"], dat$ncontrol.outcome[dat$unit == "logOR"])
dat$effective_n.exposure <- dat$samplesize.exposure

st <- psych::r.test(n = dat$effective_n.exposure, n2 = dat$effective_n.outcome, 
        r12 = sqrt(dat$rsq.exposure), r34 = sqrt(dat$rsq.outcome))
dat$steiger_dir <- dat$rsq.exposure > dat$rsq.outcome
dat$steiger_pval <- pnorm(-abs(st$z)) * 2
table(is.na(dat$rsq.exposure))
table(is.na(dat$rsq.outcome))
table(dat$rsq.exposure > dat$rsq.outcome, dat$cistrans.exposure)
table(dat$steiger_pval < 0.05, dat$cistrans.exposure)


res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))
res <- subset(res, !is.na(pval))
res$fdr <- p.adjust(res$pval, "fdr")
table(res$fdr < 0.05)

res_prs <- subset(res, paste(id.exposure, id.outcome) %in% eo_pairs)
res_prs$fdr <- p.adjust(res_prs$pval, "fdr")
table(res_prs$pval < 0.05/nrow(res_prs))

res_cis <- subset(dat, cistrans.exposure == "cis") %>% mr()
res_trans <- subset(dat, cistrans.exposure == "trans") %>% mr(method_list=c("mr_wald_ratio", "mr_ivw"))
res_trans_sf <- subset(dat, cistrans.exposure == "trans" & steiger_dir & steiger_pval < 0.05) %>% mr(method_list=c("mr_wald_ratio", "mr_ivw"))

save(dat, traits, res, res_cis, res_trans, res_trans_sf, res_prs, prs_pairs, file=here("data", "all.rdata"))



