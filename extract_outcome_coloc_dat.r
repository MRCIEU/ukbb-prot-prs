library(dplyr)
library(here)
library(glue)
library(data.table)
library(readxl)
library(tidyr)
library(ieugwasr)
select_api("private")


cl <- readRDS(here("data", "coloclist.rds"))
out <- readRDS(here("data", "regions_for_coloc_hg38.rds"))
out <- tidyr::separate(out, 
    `Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, 
    sep=":", 
    into=c("chr","pos","other_allele", "effect_allele", "imp", "v"))


cl$pos <- as.numeric(cl$pos)
a <- associations(paste0(cl$chr[1], ":", cl$pos[1]-500000, "-", cl$pos[1]+500000), cl$id.outcome[1])
o <- readRDS(here("data", "exposure_coloc_extract.rds"))
o1 <- subset(o, id.exposure==cl$exposure[1]) %>% mutate(pos.exposure=as.numeric(pos.exposure))
inner_join(o1, a, by=c("pos.exposure"="position"))

l <- list()
for(i in 1:nrow(cl)) {
    message(i)
    r <- paste0(cl$chr.exposure[i], ":", as.numeric(cl$pos.exposure)[i]-500000, "-", as.numeric(cl$pos.exposure)[i]+500000)
    l[[i]] <- associations(r, cl$id.outcome[i])
}
o <- bind_rows(l)
saveRDS(o, here("data", "outcome_coloc_extract.rds"))
