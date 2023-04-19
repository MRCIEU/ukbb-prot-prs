library(dplyr)
library(ieugwasr)
library(tidyr)
library(TwoSampleMR)
library(ggplot2)
library(here)
library(metafor)
library(GenomicRanges)
library(readxl)


load(here("data", "all.rdata"))

a <- read_xlsx(here("data", "media-2.xlsx"), sheet="ST6", skip=2)

## Coloc list
temp <- subset(dat, !duplicated(paste(SNP, id.exposure, id.outcome)))
dim(temp)
temp <- subset(temp, p.adjust(pval.outcome, "fdr") < 0.05)
dim(temp)
temp <- subset(temp, !duplicated(paste(SNP, id.exposure)))
dim(temp)
length(unique(paste(temp$SNP, temp$id.exposure)))
str(temp)
dim(temp)
temp2 <- subset(a, paste(rsID, `Assay Target`) %in% paste(temp$SNP, temp$exposure))
dim(temp2)
head(temp2)
out <- temp2 %>% select(rsID, `Assay Target`, CHROM, `GENPOS (hg38)`, `Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`) %>% mutate(start=`GENPOS (hg38)` - 500000, end = `GENPOS (hg38)` + 500000)
saveRDS(out, file=here("data", "regions_for_coloc_hg38.rds"))

# opengwas lookups

temp <- subset(dat, !duplicated(paste(SNP, id.exposure, id.outcome)))
dim(temp)
temp <- subset(temp, p.adjust(pval.outcome, "fdr") < 0.05)
dim(temp)
temp <- subset(temp, !duplicated(SNP, id.outcome))
dim(temp)
lookup <- select(temp, chr=chr.exposure, pos=pos.exposure, SNP, id=id.outcome, pval=pval.outcome) %>%
    as_tibble() %>%    
    mutate(pos=as.numeric(pos), start=pos-500000, end=pos+500000)
head(lookup)
lookupgr <- lookup %>% group_by(id, chr) %>%
    do({
        x <- reduce(IRanges(start=.$start, end=.$end))
        tibble(id=.$id[1], start=start(x), end=end(x), r=paste0(.$chr[1], ":", start, "-", end))
    })
lookupgr
l <- list()
for(i in 1:nrow(lookupgr))
{
    message(i, " of ", nrow(lookupgr))
    l[[i]] <- associations(lookupgr$r[i], lookupgr$id[i])
}
outcome_coloc_dat <- bind_rows(l)
outcome_coloc_dat <- outcome_coloc_dat %>% 
    mutate(SNP=paste0(chr, ":", position)) %>%
    select(SNP, pval.outcome=p, chr.outcome=chr, beta.outcome=beta, se.outcome=se, samplesize.outcome=n, pos.outcome=position, id.outcome=id, rsid.outcome=rsid, effect_allele.outcome=ea, other_allele.outcome=nea, eaf.outcome=eaf, outcome=trait)

dim(outcome_coloc_dat)
outcome_coloc_dat %>% group_by(outcome) %>% summarise(n=n()) %>% {summary(.$n)}


saveRDS(outcome_coloc_dat, file=here("data", "outcome_coloc_dat.rds"))

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
outcome_coloc_dat <- update_n(outcome_coloc_dat)


## Organise exposure
loc <- here("data", "exposure_ranges")
filelist <- list.files(loc)


fn <- function(f)
{
    a <- fread(file.path(loc, f)) %>%
        filter(A1FREQ > 0.01 & A1FREQ < 0.99) %>%
        tidyr::separate(ID, sep=":", into=c("chr.exposure","pos.exposure","other_allele.exposure", "effect_allele.exposure", "imp", "v")) %>%
        tidyr::separate(UKBPPP_ProteinID, sep=":", into=c("exposure", "code1", "code2", "v")) %>%
        mutate(pval=10^-LOG10P, SNP=paste0(chr, ":", pos)) %>%
        select(chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure, eaf.exposure=A1FREQ, beta.exposure=BETA, se.exposure=SE, pval.exposure=pval, exposure, id.exposure=exposure)
    a
}

e <- lapply(filelist, function(f) {message(f); fn(f)}) %>% bind_rows()

saveRDS(e, file=here("data", "exposure_coloc_dat.rds"))

out <- out %>% 
    tidyr::separate(`Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)`, sep=":", into=c("chr", "pos", "a0", "a1", "imp", "v")) %>%
    mutate(SNP=paste0(chr, ":", pos))

head(e)

code <- paste(out$`Assay Target`, out$SNP) %>% unique
e$code <- paste(e$id.exposure, e$SNP)
ind <- code %in% e$code
table(ind)
table(code %in% pas)


## harmonise

harmonise_list <- dat %>% 
    select(SNP, id.outcome, pval.outcome) %>%
    filter(!duplicated(paste(SNP, id.outcome))) %>%
    filter(p.adjust(pval.outcome, "bonferroni") < 0.05) %>%
    left_join(., select(dat, id.exposure, id.outcome, SNP, chr.exposure, pos.exposure))
head(harmonise_list)
dim(harmonise_list)

do_coloc <- function(expdat, outdat, ide, ido, chr, pos, radius, p1, p2, p12)
{
    ed <- subset(expdat, id.exposure == ide & chr.exposure == chr & pos.exposure > (pos - radius) & pos.exposure < (pos + radius))
    if(nrow(ed) == 0) return(NULL)
    od <- subset(outdat, id.outcome == ido & chr.outcome == chr & pos.outcome > (pos - radius) & pos.outcome < (pos + radius))
    if(nrow(od) == 0) return(NULL)
    d <- harmonise_data(ed, od)
    if(nrow(d) == 0) return(NULL)

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

for(i in 1:nrow(harmonise_list))
{

}





table(e$exposure %in% temp$id.exposure)
exposure_coloc_dat <- subset(e, exposure %in% temp$id.exposure)
temp <- subset(dat, !duplicated(paste(SNP, id.exposure, id.outcome)))
temp <- subset(temp, p.adjust(pval.outcome, "fdr") < 0.05)

harmonise_list <- subset(temp, !duplicated(paste(id.exposure, id.outcome)), select=c(id.exposure, id.outcome))
dim(harmonise_list)
head(harmonise_list)

l <- list()
for(id in unique(harmonise_list$id.outcome))
{
    temp <- subset(harmonise_list, id.outcome == id)
    v <- subset(e, id.exposure %in% temp$id.exposure)$SNP %>% unique()
    x <- TwoSampleMR::extract_outcome_data(v, id, proxies=FALSE)
}



lapply(1:nrow(harmonise_list), function(i)
{
    a <- subset(exposure_coloc_dat, id.exposure == harmonise_list$id.exposure[i])
    b <- associations(a$SNP, harmonise_list$id.outcome[i]) %>%
        mutate(SNP=paste0(chr, ":", position)) %>%
        select(SNP, pval.outcome=p, chr.outcome=chr, beta.outcome=beta, se.outcome=se, samplesize.outcome=n, pos.outcome=position, id.outcome=id, rsid.outcome=rsid, effect_allele.outcome=ea, other_allele.outcome=nea, eaf.outcome=eaf, outcome=trait)
    c <- harmonise_data(a, b, action=1)
    dim(c)
    harmonise_data(
        ,
        subset(outcome_coloc_dat, id.outcome == harmonise_list$id.outcome[i]),
        action=1
    ) %>% dim
    coloc_dat <- harmonise_data(e, outcome_coloc_dat)
})


