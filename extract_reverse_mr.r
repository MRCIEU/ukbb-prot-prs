library(ieugwasr)
library(dplyr)
library(readxl)
library(here)
library(tidyr)
library(data.table)
library(R.utils)
library(furrr)
library(glue)

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

lookups <- readRDS(here("data", "reverse_mr_lookups.rds"))
prots <- unique(lookups$prot)
prot_dir <- "/projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl"

alldir <- list.files(prot_dir)
head(alldir)

paths <- tibble(pdirst=grep(".tar", alldir, value=TRUE)) %>%
    mutate(pdirs = gsub(".tar", "", pdirst)) %>%
    tidyr::separate(pdirs, sep="_", into=c("prot", "c1", "c2", "v1", "type"), remove=FALSE)


paths$prot[paths$prot == "MICB"] <- "MICB_MICA"
paths$prot[paths$prot == "IL12A"] <- "IL12A_IL12B"
paths$prot[paths$prot == "DEFA1"] <- "DEFA1_DEFA1B"
paths$prot[paths$prot == "LGALS7"] <- "LGALS7_LGALS7B"
paths$prot[paths$prot == "DEFB4A"] <- "DEFB4A_DEFB4B"
paths$prot[paths$prot == "FUT3"] <- "FUT3_FUT5"
paths$prot[paths$prot == "EBI3"] <- "EBI3_IL27"

table(prots %in% paths$prot)

dim(paths)
paths <- subset(paths, paths$prot %in% prots)
dim(paths)

paths$output <- here("data", "pqtl_extract", paste0(paths$prot, ".rds"))
dim(paths)
paths <- subset(paths, !file.exists(output))
dim(paths)

stopifnot(all(file.exists(file.path(prot_dir, paths$pdirst))))

a <- subset(lookups, !duplicated(rsid))
a <- paste0("chr", a$chr, ":", a$position)
write.table(a, file="lookups_b37.txt", row=F, qu=F)

hg38 <- read.table(here("data", "hglft_genome_30ad4_c311e0.bed"), he=F)
hg38 <- tidyr::separate(hg38, V4, sep=":", into=c("chr", "pos1"))
hg38 <- tidyr::separate(hg38, pos1, sep="-", into=c("pos1", "pos2"))
hg38$chr <- gsub("chr", "", hg38$chr)
hg38 <- hg38 %>% select(chr, pos38=V2, position=pos2)
hg38$position <- as.numeric(hg38$position)
hg38$pos38 <- as.numeric(hg38$pos38)
lookups <- inner_join(lookups, hg38, by=c("chr", "position"))
dim(lookups)
head(lookups)

lookup_txt <- function(fn, pos) {
    tf <- tempfile()
    tf2 <- tempfile()
    write.table(unique(pos), file=tf, row=F, col=F, qu=F)
    cmd <- glue("zgrep -wf {tf} {fn} > {tf2}")
    system(cmd)
    fread(tf2)
}


plan(multisession, workers = 32)
furrr::future_map(1:nrow(paths), \(i)
{
    p <- paths$prot[i]
    message("Iteration ", i)
    x <- subset(lookups, prot == p)
    fnt <- file.path(prot_dir, subset(paths, prot==p)$pdirst)
    if(!file.exists(fnt))
    {
        next
    }
    cmd <- paste0("tar xvf ", fnt)
    system(cmd)
    l <- list()
    for(ch in unique(x$chr))
    {
        message(p, " ", ch)
        fn <- file.path(subset(paths, prot==p)$pdirs,
            paste0("discovery_chr", ch, "_", gsub("_", ":", subset(paths, prot==p)$pdirs), ".gz")
        )
        if(!file.exists(fn))
        {
            next
        }
        # d <- fread(fn) %>% mutate(prot=p)
        # %>% 
        #     tidyr::separate(ID, sep=":", into=c("chr", "pos", "a1", "a2", "imp", "v1")) 
        # d <- subset(d, GENPOS %in% subset(x, x$chr == ch)$pos38) %>% mutate(prot=p)
        d <- lookup_txt(fn, subset(x, x$chr == ch)$pos38)
        l[[ch]] <- d
    }
    l <- bind_rows(l)
    con <- gzfile(fn)
    names(l) <- scan(con, nlines=1, what=character())
    close(con)
    l$prot <- p
    saveRDS(l, file=paths$output[i])
    system(paste0("rm -r ", subset(paths, prot==p)$pdirs))
})

fn <- list.files(here("data", "pqtl_extract"))

pqtl_extract <- lapply(fn, \(x){
    readRDS(here("data", "pqtl_extract", x))
}) %>% bind_rows()
saveRDS(pqtl_extract, file=here("data", "pqtl_extract.rds"))


# Harmonise exposure / outcome

library(dplyr)
library(here)
library(ieugwasr)
library(TwoSampleMR)
library(tidyr)

load(here("data", "all.rdata"))
pqtl_extract <- readRDS(here("data", "pqtl_extract.rds"))
lookups <- readRDS(here("data", "reverse_mr_lookups.rds"))

head(pqtl_extract)
head(lookups)

inst <- ieugwasr::tophits(unique(prs_pairs$opengwasid))
expdat <- format_data(
    inst, "exposure", snp_col="rsid", phenotype="id", effect_allele="ea", other_allele="nea", pval="p"
)

head(pqtl_extract)
pqtl_extract <- pqtl_extract %>% tidyr::separate(ID, sep=":", into=c("chr", "pos", "nea", "ea", "imp", "v1"))
pqtl_extract$pos <- as.numeric(pqtl_extract$pos)
rsid <- subset(inst, !duplicated(rsid), select=c(rsid, chr, position))

pqtl_extract <- left_join(pqtl_extract, rsid, by=c("chr", "pos"="position"))
dim(pqtl_extract)
table(is.na(pqtl_extract$rsid))
pqtl_extract <- subset(pqtl_extract, !is.na(rsid))

outdat <- format_data(
    pqtl_extract, "outcome", snp_col="rsid", effect_allele="ea", other_allele="nea", beta="BETA", se="SE", eaf="A1FREQ", phenotype="prot"
)

prs_pairs
dat <- lapply(1:nrow(prs_pairs), \(i)
{
    harmonise_data(
        subset(expdat, exposure == prs_pairs$opengwasid[i]),
        subset(outdat, outcome == prs_pairs$prot[i]),
        action=1
    )
})
dat <- bind_rows(dat)
dat$id.outcome <- dat$outcome
dat$id.exposure <- dat$exposure
dat$exposure <- traits$code[match(dat$id.exposure, traits$opengwasid)]
dim(dat)
saveRDS(dat, file=here("data", "reverse_mr_dat.rds"))