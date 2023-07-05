# ukbb-prot-prs


1. MR + coloc results (Gib)
2. Pathway analysis (Danai)
3. Pharmaprojects lookup (Tom)
4. PRS and MR cross-disease heterogeneity (Phil)



## To run

Perform MR - cis, trans, steiger filtered etc

```
Rscript all_cis_mr.r
```

Perform coloc for cis MR results

```
Rscript extract_exposure_coloc_dat.r
Rscript run_cis_coloc.r
Rscript collate_results.r
```

Analysis of forward MR

```
quarto render docs/all_cis_mr.qmd
```

Analysis of reverse MR

```
Rscript extract_reverse_mr.r
```



Document here: https://uob-my.sharepoint.com/:w:/r/personal/ph14916_bristol_ac_uk/_layouts/15/guestaccess.aspx?e=4%3ANpXgxH&at=9&share=EebxUnfMH-lLm2HGxiKId9EB5jzG0jr8uDUstUMACBPmIQ

