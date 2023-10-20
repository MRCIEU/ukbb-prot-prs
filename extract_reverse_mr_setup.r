library(ieugwasr)
library(dplyr)
library(here)
library(readxl)

load(here("data", "all.rdata"))

traits <- read_xlsx(here("data", "Supplementary_Table5.xlsx"))

# Get instruments for all opengwasid in prs_pairs

inst <- ieugwasr::tophits(unique(traits$opengwasid))
inst
saveRDS(inst, file=here("data", "inst.rds"))

lookups <- inner_join(
    prs_pairs %>% select(prot, id=opengwasid),
    inst %>% select(chr, position, rsid, ea, nea, id),
    by="id",
    relationship="many-to-many"
)

lookups
saveRDS(lookups, file=here("data", "reverse_mr_lookups.rds"))
