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
prs_pairs <- b %>%
  dplyr::select(HGNC.symbol, opengwasid, beta=Estimate, pval=P_value) %>%
  filter(!duplicated(paste(HGNC.symbol, opengwasid)))

lookups <- inner_join(inst, prs_pairs, by=c("exposure" = "HGNC.symbol"))

outcome_dat <- lapply(unique(lookups$opengwasid), function(id)
{
  x <- subset(lookups, opengwasid == id)
  extract_outcome_data(x$SNP, x$opengwasid)
})

o <- bind_rows(outcome_dat)
dim(o)

dat <- harmonise_data(inst, o, action=1)
str(dat)
length(unique(dat$exposure))
eo_pairs <- unique(paste(prs_pairs$HGNC.symbol, prs_pairs$opengwasid))
length(eo_pairs)
dat <- subset(dat, paste(id.exposure, id.outcome) %in% eo_pairs)
dim(dat)
length(unique(paste(dat$id.exposure, dat$id.outcome)))
dat <- subset(dat, !duplicated(paste(dat$SNP, dat$id.exposure, dat$id.outcome)))
dim(dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))
res <- subset(res, !is.na(pval))

save(dat, res, prs_pairs, file=here("data", "all.rdata"))
