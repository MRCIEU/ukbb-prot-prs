# ukbb-prot-prs


1. MR + coloc results (Gib)
2. Pathway analysis (Danai)
3. Pharmaprojects lookup (Tom)
4. PRS and MR cross-disease heterogeneity (Phil)



## To run

`data/Supplementary_table5.xlsx` lists all the phenotypes to be included in the study.

Protein QTL summary stats are stored on bc4 here: `/mnt/storage/private/mrcieu/data/UKB-PPP/`.

Perform MR - cis, trans, steiger filtered etc. This uses pQTLs for 1500 proteins and then extracts outcome data from OpenGWAS. Does the MR and saves along with PRS results etc in a single Rdata file called `data/all.rdata`

```
Rscript all_cis_mr.r
```

Perform coloc for cis MR results

```
Rscript extract_exposure_coloc_dat.r
Rscript extract_outcome_coloc_dat.r
Rscript run_cis_coloc.r
Rscript collate_results.r
```

Analysis of forward MR

```
quarto render docs/all_cis_mr.qmd
```

Analysis of reverse MR. First extact data from pQTL summary stats

```
Rscript extract_reverse_mr_setup.r
Rscript extract_reverse_mr.r
quarto render docs/reverse_mr.rmd
```

Trial enrichment analysis

```
quarto render docs/drugs.qmd
```

Enrichments

```
Rscript enrichment.r
quarto render docs/enrichment.qmd
```

Note that I used singularity to run it:

```
singularity shell --bind /projects/MRC-IEU/research/projects/icep2/wp1/028/working/data/ukbb_pqtl ~/verse_latest.sif 
```



Document here: https://uob-my.sharepoint.com/:w:/r/personal/ph14916_bristol_ac_uk/_layouts/15/guestaccess.aspx?e=4%3ANpXgxH&at=9&share=EebxUnfMH-lLm2HGxiKId9EB5jzG0jr8uDUstUMACBPmIQ

