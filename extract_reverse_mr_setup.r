library(ieugwasr)
library(dplyr)
library(here)
library(readxl)
library(data.table)
library(tidyr)

load(here("data", "all.rdata"))
traits <- read_xlsx(here("data", "Supplementary_Table5.xlsx"))
prs <- fread(here("data", "Full_PRS_Results_All_Disease_All_Pvalue.tsv"))


# Get instruments for all opengwasid in prs_pairs

inst <- ieugwasr::tophits(unique(traits$opengwasid))
inst
saveRDS(inst, file=here("data", "inst.rds"))
inst <- readRDS(here("data", "inst.rds"))


stopifnot(all(traits$Disease_alt2 %in% prs$Disease_Name))
prs <- left_join(prs, select(traits, Disease_Name=Disease_alt2, opengwasid))
prs <- tidyr::separate(prs, UKBPPP_ProteinID, sep=":", into=c("prot", "v1", "v2", "v3"), remove=FALSE)

table(is.na(prs$opengwasid))
prs <- subset(prs, !is.na(opengwasid))
prs$fdr <- p.adjust(prs$P_value, "fdr")
prs$sig <- prs$fdr < 0.05
table(prs$sig)
head(prs)

lookups <- inner_join(
    prs %>% filter(sig) %>% select(prot, id=opengwasid),
    inst %>% select(chr, position, rsid, ea, nea, id),
    by="id",
    relationship="many-to-many"
)

dim(lookups)
head(lookups)
saveRDS(lookups, file=here("data", "reverse_mr_lookups.rds"))
