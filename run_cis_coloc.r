library(dplyr)
library(ggplot2)
library(here)
library(ieugwasr)
library(purrr)
library(TwoSampleMR)

expdat <- readRDS(here("data", "exposure_coloc_dat.rds"))
expdat$pos.exposure <- as.numeric(expdat$pos.exposure)

load(here("data", "all.rdata"))
dat$pos.exposure <- as.numeric(dat$pos.exposure)

dat_cis <- subset(dat, cistrans.exposure == "cis") %>%
    filter(p.adjust(pval.outcome, "fdr") < 0.05)
dim(dat_cis)
all(unique(dat_cis$exposure) %in% expdat$id.exposure)

dat_cis$snpid <- paste0(dat_cis$chr.exposure, ":", dat_cis$pos.exposure)
table(unique(dat_cis$snpid) %in% expdat$SNP)

length(unique(expdat$exposure))

convert_to_outcome <- function(d)
{
    d %>% mutate(SNP=paste0(chr, ":", position)) %>% 
    dplyr::select(SNP,
        id.outcome=id, 
        outcome=trait, 
        beta.outcome=beta, 
        se.outcome=se, rsid, samplesize.outcome=n, effect_allele.outcome=ea, other_allele.outcome=nea, eaf.outcome=eaf)
}

make_coloc_dat <- function(i, dat_cis, expdat)
{
    r <- paste0(dat_cis$chr.exposure[i], ":", dat_cis$pos.exposure[i]-500000, "-", dat_cis$pos.exposure[i]+500000)
    a <- associations(r, dat_cis$id.outcome[i]) %>% convert_to_outcome
    b <- subset(expdat, 
        id.exposure == dat_cis$id.exposure[i] & 
        pos.exposure > (dat_cis$pos.exposure[i] - 500000) &
        pos.exposure < (dat_cis$pos.exposure[i] + 500000) &
        chr.exposure == dat_cis$chr.exposure[i]
    )
    d <- suppressMessages(harmonise_data(b, a, action=1))
    if(nrow(d) == 0) return(NULL)
    d$target_snp <- paste0(dat_cis$chr.exposure[i], ":", dat_cis$pos.exposure[i])
    return(d)
}

update_n <- function(x) {
    x <- x %>% mutate(outcome=gsub("EBI", "ebi", outcome), id.outcome = outcome)
    if("samplesize.outcome" %in% names(x))
    {
        x$samplesize.outcome[is.na(x$samplesize.outcome)] <- mean(x$samplesize.outcome, na.rm=T)
    }
    x$eaf.outcome[is.na(x$eaf.outcome)] <- 0.2
    x$eaf.outcome[x$eaf.outcome == 0] <- 0.001
    x$eaf.outcome[x$eaf.outcome == 1] <- 0.999
    if("ncase.outcome" %in% names(x))
    {
        x$ncase.outcome[x$id.outcome == "ieu-b-2"] <- 21982
        x$ncontrol.outcome[x$id.outcome == "ieu-b-2"] <- 41944
        x$ncase.outcome[x$id.outcome == "ieu-a-1024"] <- 9722
        x$ncontrol.outcome[x$id.outcome == "ieu-a-1024"] <- 17376
        x$ncase.outcome[is.na(x$ncase.outcome)] <- mean(x$ncase.outcome, na.rm=T)
        x$ncontrol.outcome[is.na(x$ncontrol.outcome)] <- mean(x$ncontrol.outcome, na.rm=T)
    }
    return(x)
}

do_coloc <- function(d, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5){
    if(is.null(d)) return(NULL)
    d <- update_n(d)
    d1 <- list(beta=d$beta.exposure, varbeta=d$se.exposure^2, MAF=d$eaf.exposure, N=31000, type="quant")

    d2 <- list(beta=d$beta.outcome, varbeta=d$se.outcome^2, N=mean(d$samplesize.outcome), MAF=d$eaf.outcome)
    s <- mean(d$ncase.outcome / d$samplesize.outcome)

    if(is.na(s))
    {
        d2$type <- "quant"
    } else {
        d2$type <- "cc"
        d2$s <- s
    }
    res <- coloc.abf(d1, d2, p1, p2, p12)
    as.list(res$summary) %>% as_tibble() %>% mutate(snp=d$SNP[1], exposure=d$exposure[i], id.exposure=d$id.exposure[i], outcome=d$outcome[i], id.outcome=d$id.outcome[i]) %>% return()
}

plot_res <- function(d, m="split") {
    p1 <- d %>% 
        mutate(ze=beta.exposure / se.exposure, zo=beta.outcome/se.outcome) %>%
        dplyr::select(pos.exposure, ze, zo) %>% 
        tidyr::gather("key", "value", ze, zo) %>% 
        mutate(key = case_when(key == "ze" ~ d$exposure[1], key == "zo" ~ d$outcome[1], TRUE ~ "")) %>%
        ggplot(., aes(x=pos.exposure, y=abs(value))) +
            geom_point(aes(colour=key))

    # p1 <- d %>% 
    #     dplyr::select(pval.exposure, pval.outcome, pos.exposure) %>%
    #     tidyr::gather(., "key", "value", pval.exposure, pval.outcome) %>%
    #     ggplot(., aes(x=pos.exposure, y=-log10(value))) +
    #     geom_point(aes(colour=key))
    if(m == "split")
    {
        p1 <- p1 + facet_grid(key ~ ., scale="free_y") + theme(legend.position="none")
    }
    p1
}

o <- purrr::map(1:nrow(dat_cis), make_coloc_dat, dat_cis=dat_cis, expdat=expdat, .progress=TRUE)
length(o)
sapply(o, nrow)

res <- purrr::map(o, do_coloc, .progress=TRUE) %>% bind_rows()
res
table(res$PP.H4.abf > 0.8)

purrr::map(1:nrow(dat_cis), \(i){
    if(!is.null(o[[i]]))
    {
        plot_res(o[[i]]) %>% ggsave(., file=here("images", paste0(o[[i]]$id.exposure[1], "_", o[[i]]$id.outcome[1], ".pdf")))
    }
}, .progress=TRUE)


saveRDS(res, file=here("data", "cis_coloc_res.rds"))
saveRDS(o, file=here("data", "cis_coloc_dat.rds"))


