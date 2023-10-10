library(dplyr)
library(here)
library(glue)
library(data.table)
library(readxl)
library(tidyr)

load(here("data", "all.rdata"))
a <- read_xlsx(here("data", "media-2.xlsx"), sheet="ST6", skip=2)

ukbpppdir <- "/mnt/storage/private/mrcieu/data/UKB-PPP/"
protmap <- fread(file.path(ukbpppdir, "Metadata/Protein annotation/olink_protein_map_1.5k_v1.tsv"))
ssdir <- file.path(ukbpppdir, "UKB-PPP pGWAS summary statistics/European (discovery)/")
protmap$path <- file.path(ssdir, paste0(gsub(":", "_", protmap$UKBPPP_ProteinID), "_", protmap$Panel, ".tar"))
stopifnot(all(file.exists(protmap$path)))

tt

tt <- subset(res_cis, paste(id.exposure, id.outcome) %in% paste(prs_pairs$prot, prs_pairs$opengwasid))
tt <- subset(tt, p.adjust(pval, "fdr") < 0.05)
dim(tt)

tt2 <- subset(res_cis, p.adjust(pval, "fdr") < 0.05)
dim(tt2)

tt3 <- bind_rows(tt, tt2)
tt3 <- subset(tt3, !duplicated(paste(id.exposure, id.outcome)))
dim(tt3)

head(tt3)

tt4 <- subset(dat, paste(id.exposure, id.outcome) %in% paste(tt3$id.exposure, tt3$id.outcome) & cistrans.exposure == "cis")
saveRDS(tt4, file=here("data", "coloclist.rds"))

tt5 <- subset(a, paste(rsID, `Assay Target`) %in% paste(tt4$SNP, tt4$exposure))
out <- tt5 %>% select(rsID, `Assay Target`, CHROM, `GENPOS (hg38)`, `Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`) %>% mutate(start=`GENPOS (hg38)` - 500000, end = `GENPOS (hg38)` + 500000)
saveRDS(out, file=here("data", "regions_for_coloc_hg38.rds"))

table(out$`Assay Target` %in% protmap$Assay)


out
out <- left_join(out, protmap %>% select(Assay, path), by=c("Assay Target"="Assay"))

extract <- function(path, chr, start, end) {
    td <- tempdir()
    glue("tar xvf '{path}' -C {td}") %>% system()
    p2 <- gsub(".tar", "", path) %>% basename() %>% file.path(td, .)

    fn <- list.files(p2, full.names=TRUE) %>% grep(paste0("chr", chr, "_"), ., value=TRUE)
    d <- fread(fn) %>% 
        filter(GENPOS > start & GENPOS < end) %>%
        mutate(prot=basename(fn) %>% gsub(paste0("discovery_chr", chr, "_"), "", .)) %>%
        filter(A1FREQ > 0.01 & A1FREQ < 0.99) %>%
        tidyr::separate(ID, sep=":", into=c("chr.exposure","pos.exposure","other_allele.exposure", "effect_allele.exposure", "imp", "v")) %>%
        tidyr::separate(prot, sep=":", into=c("exposure", "code1", "code2", "v", "panel")) %>%
        mutate(pval=10^-LOG10P, SNP=paste0(chr.exposure, ":", pos.exposure)) %>%
        select(SNP, chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure, eaf.exposure=A1FREQ, beta.exposure=BETA, se.exposure=SE, pval.exposure=pval, exposure=exposure, id.exposure=exposure)
    unlink(p2)
    d
}

# dir.create(here("data", "coloc_exposure_extract"), recursive=TRUE)

l <- list()
for(i in 120:nrow(out)) {
    l[[i]] <- extract(out$path[i], out$CHROM[i], out$start[i], out$end[i])
}
o <- bind_rows(l)

saveRDS(o, file=here("data", "exposure_coloc_extract.rds"))
# saveRDS(e, file=here("data", "exposure_coloc_dat.rds"))
