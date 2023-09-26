library(dplyr)
library(tidyr)
library(here)
library(data.table)

## Organise exposure
loc <- here("data", "exposure_ranges")
filelist <- list.files(loc)

fn <- function(f)
{
    a <- fread(file.path(loc, f)) %>%
        filter(A1FREQ > 0.01 & A1FREQ < 0.99) %>%
        tidyr::separate(ID, sep=":", into=c("chr.exposure","pos.exposure","other_allele.exposure", "effect_allele.exposure", "imp", "v")) %>%
        tidyr::separate(UKBPPP_ProteinID, sep=":", into=c("exposure", "code1", "code2", "v")) %>%
        mutate(pval=10^-LOG10P, SNP=paste0(chr.exposure, ":", pos.exposure)) %>%
        select(SNP, chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure, eaf.exposure=A1FREQ, beta.exposure=BETA, se.exposure=SE, pval.exposure=pval, exposure=exposure, id.exposure=exposure)
    a
}

e <- lapply(filelist, function(f) {message(f); fn(f)}) %>% bind_rows()

saveRDS(e, file=here("data", "exposure_coloc_dat.rds"))
