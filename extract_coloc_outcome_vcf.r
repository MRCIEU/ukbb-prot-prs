library(gwasvcf)
library(gwasglue)
library(TwoSampleMR)
library(dplyr)
library(here)
library(coloc)
library(GenomicRanges)
library(ggplot2)
set_bcftools()

# set dir
vcfdir <- "/mnt/storage/private/mrcieu/data/IGD/data/public"

# load mr data
load(here("data", "all.rdata"))

# load exposure data
exposure_dat <- readRDS(here("data", "exposure_coloc_dat.rds"))

# get the snp-exposure-outcome sets to colocalise
harmonise_list <- dat %>% 
    select(SNP, id.outcome, pval.outcome) %>%
    filter(!duplicated(paste(SNP, id.outcome))) %>%
    filter(p.adjust(pval.outcome, "bonferroni") < 0.05) %>%
    left_join(., select(dat, id.exposure, id.outcome, SNP, chr.exposure, pos.exposure))
head(harmonise_list)
dim(harmonise_list)

# extract outcome data
extract_coloc_data <- function(vcffile, chr, pos, radius)
{
    a <- GRanges(chr, IRanges(pos-radius, pos+radius), target=paste0(chr, ":", pos))
    b <- reduce(a)
    query_gwas(vcffile, as.data.frame(b) %>% rename(chrom=seqnames))
    #  %>%
    #     gwasglue::gwasvcf_to_TwoSampleMR(., type="outcome") %>%
    #     mutate(rsid.outcome=SNP, SNP=paste0(chr.outcome, ":", pos.outcome), target=paste0(chr, ":", pos))
}

update_n <- function(x)
{
    x <- x %>% mutate(outcome=gsub("EBI", "ebi", outcome), id.outcome = outcome)
    x$samplesize.outcome[x$id.outcome == "ieu-b-2"] <- 63926
    x$ncase.outcome[x$id.outcome == "ieu-b-2"] <- 21982
    x$ncontrol.outcome[x$id.outcome == "ieu-b-2"] <- 41944
    x$ncase.outcome[x$id.outcome == "ieu-a-1024"] <- 9722
    x$ncontrol.outcome[x$id.outcome == "ieu-a-1024"] <- 17376
    x$ncontrol.outcome
    return(x)
}

outcome_coloc_dat <- lapply(unique(harmonise_list$id.outcome), function(id) {
    message(id)
    h <- subset(harmonise_list, id.outcome == id)
    chr <- h$chr.exposure
    pos <- as.numeric(h$pos.exposure)
    extract_coloc_data(
        file.path(vcfdir, id, paste0(id, ".vcf.gz")),
        chr,
        pos,
        500000
    ) %>%
        gwasglue::gwasvcf_to_TwoSampleMR(., type="outcome") %>%
        mutate(
            rsid.outcome=SNP, 
            SNP=paste0(chr.outcome, ":", pos.outcome), 
            target=paste0(chr.outcome, ":", pos.outcome))
}) %>% bind_rows()

outcome_coloc_dat <- update_n(outcome_coloc_dat)
saveRDS(outcome_coloc_dat, file=here("data", "outcome_coloc_dat.rds"))

do_coloc <- function(harmonise_list, row, exposure_dat, outcome_dat, radius = 500000, p1 = 1, p2 = 1, p12 = 1){
    chr <- harmonise_list$chr.exposure[row]
    pos <- harmonise_list$pos.exposure[row]
    eid <- harmonise_list$id.exposure[row]
    oid <- harmonise_list$id.outcome[row]

    e <- subset(exposure_dat, 
        id.exposure == eid &
        chr.exposure == chr &
        pos.exposure > (pos - radius) &
        pos.exposure < (pos + radius)
    )
    if(nrow(e) == 0) return(NULL)
    o <- subset(outcome_dat, 
        id.outcome == oid &
        chr.outcome == chr &
        pos.outcome > (pos - radius) &
        pos.outcome < (pos + radius)
    )
    if(nrow(o) == 0) return(NULL)
    d <- harmonise_data(e, o, action=1)
    if(nrow(d) == 0) return(NULL)

    d1 <- list(beta=e$beta.exposure, varbeta=e$se.exposure^2, MAF=e$eaf.exposure, N=31000, type="quant")

    d2 <- list(beta=d$beta.outcome, varbeta=d$se.outcome^2, N=mean(d$samplesize.outcome), MAF=d$eaf.outcome)
    s <- mean(d$ncase.outcome / d$samplesize.outcome)
    print(s)
    if(is.na(s))
    {
        d2$type <- "quant"
    } else {
        d2$type <- "cc"
        d2$s <- s
    }
    res <- coloc.abf(d1, d2, p1, p2, p12)
    out <- list(
        colocres = as.list(res$summary) %>% as_tibble() %>% {cbind(harmonise_list[row,], .)},
        dat = d
    )
    return(out)
}

plot_res <- function(d, m="split")
{
    p1 <- d %>% 
        mutate(ze=beta.exposure / se.exposure, zo=beta.outcome/se.outcome) %>%
        dplyr::select(pos.exposure, ze, zo) %>% 
        tidyr::gather("key", "value", ze, zo) %>% 
        ggplot(., aes(x=pos.exposure, y=abs(value))) +
            geom_point(aes(colour=key))

    # p1 <- d %>% 
    #     dplyr::select(pval.exposure, pval.outcome, pos.exposure) %>%
    #     tidyr::gather(., "key", "value", pval.exposure, pval.outcome) %>%
    #     ggplot(., aes(x=pos.exposure, y=-log10(value))) +
    #     geom_point(aes(colour=key))
    if(m == "split")
    {
        p1 <- p1 + facet_grid(key ~ ., scale="free_y")
    }
    p1
}

# missing eaf values
harmonise_list$pos.exposure <- as.numeric(harmonise_list$pos.exposure)
table(outcome_coloc_dat$id.outcome, is.na(outcome_coloc_dat$eaf.outcome))
outcome_coloc_dat$eaf.outcome[is.na(outcome_coloc_dat$eaf.outcome)] <- 0.2

# Which pairs actually have data
harmonise_list_r <- subset(harmonise_list, paste0(id.exposure, " ", chr.exposure, ":", pos.exposure) %in% paste(exposure_dat$exposure, exposure_dat$SNP))
dim(harmonise_list_r)

# run coloc
colocres <- list()
for(i in 1:nrow(harmonise_list_r))
{
    message(i)
    colocres[[i]] <- do_coloc(harmonise_list_r, i, exposure_dat, outcome_coloc_dat)
}

for(i in 1:nrow(harmonise_list_r))
{
    message(i)
    p <- colocres[[i]]$dat %>% plot_res
    ggsave(p, file=here("images", paste0(harmonise_list_r$id.exposure[i], "_", harmonise_list_r$id.outcome[i], ".pdf")))
}

lapply(colocres, function(x) x$colocres) %>% bind_rows %>% saveRDS(., file=here("data", "colocres.rds"))

####

## deprecated

l <- list()
for(id in unique(harmonise_list$id.outcome))
{
    temp <- subset(harmonise_list, id.outcome == id)
    v <- subset(exposure_dat, id.exposure %in% temp$id.exposure)
    x <- query_gwas(file.path(vcfdir, id, paste0(id, ".vcf.gz")), chrompos=unique(v$SNP))    
    x1 <- gwasglue::gwasvcf_to_TwoSampleMR(x, type="outcome") %>%
        mutate(rsid.outcome=SNP, SNP=paste0(chr.outcome, ":", pos.outcome))
    d <- harmonise_data(v, x1, action=1)
    l[[id]] <- d
}

saveRDS(l, file=here("data", "gwasvcf_extract.rds"))


update_n <- function(x)
{
    x <- bind_rows(l) %>%
        mutate(outcome=gsub("EBI", "ebi", outcome), id.outcome = outcome)
    x$samplesize.outcome[x$id.outcome == "ieu-b-2"] <- 63926
    x$ncase.outcome[x$id.outcome == "ieu-b-2"] <- 21982
    x$ncontrol.outcome[x$id.outcome == "ieu-b-2"] <- 41944
    x$ncase.outcome[x$id.outcome == "ieu-a-1024"] <- 9722
    x$ncontrol.outcome[x$id.outcome == "ieu-a-1024"] <- 17376
    x$ncontrol.outcome
    return(x)
}


table(coloc_dat$id.outcome, is.na(coloc_dat$ncase.outcome))

group_by(coloc_dat, id.outcome) %>% summarise(s=mean(ncase.outcome/ncontrol.outcome))
table(coloc_dat$id.outcome, mean(coloc_dat$ncase.outcome/coloc_dat$ncontrol.outcome))
table(coloc_dat$id.outcome, is.na(coloc_dat$eaf.))

do_coloc <- function(coloc_dat, p1, p2, p12)
{
    d1 <- list(beta=coloc_dat$beta.exposure, varbeta=coloc_dat$se.exposure^2, MAF=coloc_dat$eaf.exposure, N=31000, type="quant")

    d2 <- list(beta=coloc_dat$beta.outcome, varbeta=coloc_dat$se.outcome^2, N=mean(coloc_dat$samplesize.outcome), MAF=coloc_dat$eaf.outcome)
    s <- mean(coloc_dat$ncase.outcome / coloc_dat$samplesize.outcome)
    print(s)
    if(is.na(s))
    {
        d2$type <- "quant"
    } else {
        d2$type <- "cc"
        d2$s <- s
    }
    coloc.abf(d1, d2, p1, p2, p12)
}

cd <- subset(coloc_dat, paste(id.exposure, id.outcome) %in% unique(paste(dat$id.exposure, dat$id.outcome)))

dim(coloc_dat)
dim(cd)

summary(dat$pval.outcome)
dats <- dat %>%
    select(id.outcome, pval.outcome, chr.exposure, pos.exposure) %>% 
    filter(!duplicated(paste(chr.exposure, pos.exposure, id.outcome))) %>%
    filter(pval.outcome < 0.05 / n()) %>%
    left_join(., select(dat, id.exposure, id.outcome, chr.exposure, pos.exposure))
dim(dats)
head(dats)
subset(dat %>% filter(!duplicated(paste(SNP, id.exposure))), paste(chr.exposure, pos.exposure) == "1 109817590")

lapply(1:nrow(dats), function(i)
{
    subset(coloc_dat, 
        id.exposure == dats$id.exposure[i] &
        id.outcome == dats$id.outcome[i] &
        chr.exposure == dats$chr.exposure[i] &
        pos.exposure > (dats$pos.exposure[i] - 500000) &
        pos.exposure < (dats$pos.exposure[i] + 500000)
    )
})


npair <- length(unique(paste(dat$id.exposure, dat$id.outcome)))
thresh <- 0.05/npair

group_by(coloc_dat, id.exposure, id.outcome) %>%


n <- 
coloc_dat <- group_by(id.exposure, id.outcome)

a <- subset(coloc_dat, id.exposure == id.exposure[1] & id.outcome == id.outcome[1])

library(ggplot2)
p1 <- a %>% select(pos.exposure, pval.exposure, pval.outcome) %>% tidyr:: gather(., key="key", value="value", pval.exposure, pval.outcome) %>% ggplot(., aes(x=pos.exposure, y=-log10(value))) +
geom_point(aes(colour=key))
ggsave(p1, file="coloc.pdf")



o <- do_coloc(a, 1, 1, 1)
names(o)
o$summary

summary(a$ncase.outcome)
summary(a$ncontrol.outcome)
