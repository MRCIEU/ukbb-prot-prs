library(ieugwasr)
library(dplyr)
library(readxl)
library(here)
library(tidyr)
library(data.table)

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

prots[!prots %in% paths$prot]


l <- list()
i <- 1
for(p in prots)
{
    message("Iteration ", i)
    x <- subset(lookups, prot == p)
    fnt <- file.path(prot_dir, subset(paths, prot==p)$pdirst)
    if(!file.exists(fnt))
    {
        next
    }
    cmd <- paste0("tar xvf ", fnt)
    system(cmd)
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
        d <- fread(fn) %>% 
            tidyr::separate(ID, sep=":", into=c("chr", "pos", "a1", "a2", "imp", "v1")) %>% mutate(prot=p)
        d <- subset(d, pos %in% subset(x, x$chr == ch)$position)
        l[[i]] <- d
        i <- i + 1
    }
    system(paste0("rm -r ", subset(paths, prot==p)$pdirs))
}

length(unique(paste(lookups$prot, lookups$chr)))
table(table(paste(lookups$prot, lookups$chr)))


for(i in 1:10)
{
    if(i > 5)
    {
        next
    }
    print(i)
}
