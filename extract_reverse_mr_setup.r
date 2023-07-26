library(ieugwasr)
library(dplyr)
library(here)

load(here("data", "all.rdata"))

head(prs_pairs)

# Get instruments for all opengwasid in prs_pairs

inst <- ieugwasr::tophits(unique(prs_pairs$opengwasid))
inst

lookups <- left_join(
    prs_pairs %>% select(prot, id=opengwasid),
    inst %>% select(chr, position, rsid, ea, nea, id),
    by="id",
    relationship="many-to-many"
)

lookups
saveRDS(lookups, file=here("data", "reverse_mr_lookups.rds"))
