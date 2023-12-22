library(readxl)
library(here)
library(tidyr)
library(dplyr)
# cd ~/ukbb-prot-prs/data
setwd("~/ukbb-prot-prs/data")
# Library
library(hrbrthemes)
library(ggplot2)
# a2<-read.table("Full_PRS_Results_All_Disease_All_Pvalue.txt",head=TRUE,sep="\t",stringsAsFactors=FALSE)

b<-load_and_format_protein_info()
a<-load_and_format_prs_results(protein_info=b)
a2<-a[which(a$P_value<2.30e-6),]
length(unique(a$Disease))
length(unique(a2$Disease))

rmr<-load_and_format_rmr_results(protein_info=b)
fmr<-load_and_format_fmr_results(protein_info=b)

head(rmr)

join_and_format_all_results<-funcion(a=NULL){
  d<-merge(a,rmr[,c("z","pval","ppid","fdr")],by="ppid",all.x=TRUE)
  names(d)[names(d) == "z"]<-"z_rmr"
  names(d)[names(d) == "fdr"]<-"fdr_rmr"
  names(d)[names(d) == "pval"]<-"pval_rmr"
  d<-merge(d,fmr[,c("ppid","cistrans","method","z","pval","fdr")],by="ppid",all.x=TRUE)

  names(d)[names(d) == "method"] <-"method_fmr"
  names(d)[names(d) == "z"] <-"z_fmr"
  names(d)[names(d) == "pval"] <-"pval_fmr"
  names(d)[names(d) == "fdr"] <-"fdr_fmr"

  d$prs_p_sig<-FALSE
  Pos<-which(d$P_value<2.30e-6)
  d$prs_p_sig[Pos]<-TRUE

  d$rmr_p_sig<-NA
  Pos<-which(d$pval_rmr<0.05)
  d$rmr_p_sig[Pos]<-TRUE
  
  d$fmr_p_sig<-NA
  Pos<-which(d$fdr_fmr<0.05)
  d$fmr_p_sig[Pos]<-TRUE

  write.table(d,"~/ukbb-prot-prs/data/prs_rmr_fmr.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

  d$Fill<-NA
  d$Fill

  d$Fill[d$prs_p_sig]<-1
  d$Fill[!d$prs_p_sig]<-0
  d$Fill[which(d$prs_p_sig & d$rmr_p_sig)]<-2
  d$Fill[which(d$prs_p_sig & d$rmr_p_sig & d$fmr_p_sig)]<-3
  d$Fill[which(d$fmr_p_sig & d$prs_p_sig & is.na(d$rmr_p_sig))]<-4
  d$Fill[which(d$fmr_p_sig & !d$prs_p_sig )]<-5

  d$outcome<-d$Disease
  d$outcome[d$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"BC"
  d$outcome[d$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"ER-BC"
  d$outcome[d$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"ER+BC"
  d$outcome[d$outcome=="Endometrial cancer"]<-"END"
  d$outcome[d$outcome=="Epilepsy"]<-"EPI"
  d$outcome[d$outcome=="High grade and low grade serous ovarian cancer"]<-"SOC"
  d$outcome[d$outcome=="High grade serous ovarian cancer"]<-"HgSOC"
  d$outcome[d$outcome=="Lung cancer"]<-"LC"
  d$outcome[d$outcome=="Lung adenocarcinoma"]<-"LAD"
  d$outcome[d$outcome=="Ovarian cancer"]<-"OC"
  d$outcome[d$outcome=="Oral cavity and pharyngeal cancer"]<-"OCP"
  d$outcome[d$outcome=="Prostate cancer (advanced)"]<-"AdPC"
  d$outcome[d$outcome=="Prostate cancer (early-onset)"]<-"EoPC"
  d$outcome[d$outcome=="Squamous cell lung cancer"]<-"SqLC"
  d$outcome[d$outcome=="Systemic lupus erythematosus"]<-"SLE"
  d$outcome[d$outcome=="Pancreatic cancer"]<-"PanC"
  d$outcome[d$outcome=="Prostate cancer"]<-"ProC"

  d<-d[which(d$outcome !="AdPC"),]
  d<-d[which(d$outcome !="EoPC"),]


unique(d$outcome)
table(d$UniProt)
head(d)
which(is.na(d$Fill))

# Heatmap
head(d) 
ggplot(d, aes(Disease, UniProt, fill= Fill)) + 
  geom_tile()


ggplot(d, aes(outcome, UniProt, fill= Fill)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="blue") +
  theme_ipsum()


# d2[d2$outcome == "SLE",c("prs_p_sig","fmr_p_sig")]
# exclude outcomes with 0 sig prs results
Table<-table(d2$outcome,d2$prs_p_sig)
Excl<-names(which(Table[,2]==0))
d2<-d2[!d2$outcome %in% Excl,]

colors <- c("white","black","blue","green", "yellow", "red")
ggplot(d, 
    aes(x = UniProt, y =outcome , fill = factor(Fill))) + 
    geom_tile() + 
    scale_fill_manual(values=colors)

d2<-d[d$prs_p_sig, ]
colors <- c("black","blue","green","red")
ggplot(d2, 
    aes(x = UniProt, y =outcome , fill = factor(Fill))) + 
    geom_tile() + 
    scale_fill_manual(values=colors)

table(d2$Fill)
d2[,c("UniProt","outcome","Fill")]

head(d2)
coloc<-read.csv("cis_coloc_results.csv",head=TRUE,stringsAsFactors=FALSE)


# Dummy data
x <- LETTERS[1:20]
y <- paste0("var", seq(1,20))
data <- expand.grid(X=x, Y=y)
data$Z <- runif(400, 0, 5)
 
# Heatmap 
ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()



load_and_format_prs_results<-function(protein_info=NULL){
  a<-read.table("Full_PRS_Results_All_Disease_All_Pvalue.tsv",head=TRUE,sep="\t",stringsAsFactors=FALSE,quote="")
  # a4<-a[grep("ieu-a-985",a$Disease_Name),]
  # a3<-a[a$P_value<2.30e-6,]
  # head(a4)
  # unique(a3$Disease_Name)
  # min(a4$P_value)
  head(a)
  a<-a[,c(1:6,8:9)]
  # table(a$Disease)
  # a3<-a[a$P_value<2.30e-6,]
  # unique(a3$Disease)
  # a$id<-NA
  # a2<-a[grep("merged",a$Disease),]
  # table(a$Disease_Name)
  # a2$id<-unlist(strsplit(a2$Disease,split="_merged"))
  # a1<-a[grep("merged",a$Disease,invert=TRUE),]
  # ieugwasr::get_access_token()
  # # ieugwasr::check_access_token()
  # ao<-ieugwasr::gwasinfo()
  # IDS<-unique(a2$id) 
  # ao2<-ao[ao$id %in% IDS, c("trait","id","note")]
  # ao2$trait[which(ao2$note=="Lung cancer in ever smokers" ) ]<-"Lung cancer in ever smokers"
  # unique(a2$Disease)
  # length(unique(a$Disease))
  # a2<-merge(a2,ao2,by="id")
  # a2$Disease<-a2$trait 

  # a<-plyr::rbind.fill(a1,a2)
  # a<-a[!a$Disease %in% c("diastolic blood pressure","HDL cholesterol","LDL cholesterol","Total cholesterol","Triglycerides","systolic blood pressure"),]
  # a3<-a[which(a$P_value<2.30e-6),]
  
  proteins<-unlist(strsplit(a$UKBPPP_ProteinID,":"))
  proteins<-proteins[seq(2,length(proteins),by=4)]
  a$UniProt<-proteins
  a<-a[!a$Disease_Name %in%  c("Early-onset prostate cancer",  "Advanced prostate cancer"),    ]
  a.m<-merge(a,b,by.x="UKBPPP_ProteinID",by.y="UKBPPP.ProteinID")
  a.m$Disease<-a.m$Disease_Name
  a.m$Disease[a.m$Disease=="Coronary artery disease"] <-"CAD"
  a.m$Disease[a.m$Disease=="Type 2 Diabetes"] <-"T2D"
  a.m$Disease[a.m$Disease=="Chronic kidney disease"] <-"CKD"
  a.m$Disease[a.m$Disease=="Lung cancer (in ever smokers)"] <-"LCsmk"
  a.m$Disease[a.m$Disease == "Schizophrenia"]<-"SCZ"
  a.m$Disease[a.m$Disease == "Alzheimer’s disease"]<-"AD"
  a.m$Disease[a.m$Disease == "Parkinson’s disease"]<-"PD"
  a.m$Disease[a.m$Disease == "Amyotrophic lateral sclerosis" ]<-"ALS"
  a.m$Disease[a.m$Disease == "Multiple sclerosis"   ]<-"MS"
  a.m$Disease[a.m$Disease == "Major depressive disorder"  ]<-"MDD"
  a.m$ppid<-paste(a.m$UKBPPP_ProteinID,a.m$Disease)
  

  return(a.m)
}

# a<-a.m


load_and_format_protein_info<-function(){
  xls <- read_xlsx(here("~/ukbb-prot-prs/data", "media-2.xlsx"), sheet="ST3", skip=2)

  b<-data.frame(xls)
  b<-unique(b[,c("UniProt" ,"UKBPPP.ProteinID","Assay.Target","Olink.ID")])
  # b<-unique(b[,c("UniProt","Assay.Target")])
  Dups<-b$UniProt[duplicated(b$UniProt)]
  b[b$UniProt == Dups[1],]
  b$UniProt<-gsub(";","_",b$UniProt)
  return(b)
}

load_and_format_rmr_results<-function(protein_info=NULL){
  
  # readRDS("reverse_mr.rds")
  
  rmr<-readRDS("reverse_mr_outliers.rds")
  # rmr<-read.csv("~/ukbb-prot-prs/results/resnohet.csv",head=TRUE,stringsAsFactors=FALSE)
  rmr2<-merge(rmr,b,by.x="outcome",by.y="Assay.Target")
  rmr2$id2<-paste(rmr2$outcome,rmr2$exposure)
  # duplicates introduced because one some proteins (as indexed by Assay.Target or UniProt) have multiple olink IDs. The results (not surprising) and protein name (reassuring?) are same across the olink IDs. 
  # Dups<-unique(rmr2$id2[duplicated(rmr2$id2)])
  # rmr2[rmr2$id2 %in% Dups[2],]
  rmr2<-rmr2[!duplicated(rmr2$id2),]
  rmr2<-rmr2[!rmr2$exposure %in% c("Diastolic Blood Pressure","HDLc","LDLc","Systolic Blood Pressure","Total Cholesterol","Triglycerides"),]
  rmr2$exposure[rmr2$exposure=="Lung cancer (ever smokers)"] <-"LCsmk"
  rmr2$exposure[rmr2$exposure== "Schizophrenia"]<-"SCZ"
  rmr2$ppid<-paste(rmr2$UKBPPP.ProteinID,rmr2$exposure)
  ppid1<-unique(rmr2$ppid)
  rmr2$z<-rmr2$b/rmr2$se
  rmr2$fdr<-rmr2$pval<-p.adjust(rmr2$pval, "fdr")

  # ppid2<-unique(a$ppid)
  # ppid1[!ppid1 %in% ppid2]  
  return(rmr2)
}

load_and_format_fmr_results<-function(protein_info=b){
  fmr<-read.csv("~/ukbb-prot-prs/data/all_mr_results.csv",head=TRUE,stringsAsFactors=FALSE)
  # all fmr proteins are nested within b 
  fmr2<-merge(fmr,b,by.x="exposure",by.y="Assay.Target")
  prot1<-unique(fmr$exposure)
  prot2<-unique(b$Assay.Target)
  fmr2$id2<-paste(fmr2$exposure,fmr2$outcome,fmr2$cistrans)
  Dups<-unique(fmr2$id2[duplicated(fmr2$id2)])
  fmr2<-fmr2[!duplicated(fmr2$id2),]
  # dis1<-unique(fmr2$outcome)
  # dis2<-unique(d$Disease)
  # dis1[!dis1 %in% dis2]
  # table(a$Disease)
  fmr2$outcome[fmr2$outcome == "Parkinson's disease"]<-"PD"
  fmr2$outcome[fmr2$outcome == "multiple sclerosis"]<-"MS"
  fmr2$outcome[fmr2$outcome == "Major depressive disorder"]<-"MDD"
  fmr2$outcome[fmr2$outcome == "schizophrenia"]<-"SCZ"
  fmr2$outcome[fmr2$outcome == "Amyotrophic lateral sclerosis"]<-"ALS"
  fmr2$outcome[fmr2$outcome == "Alzheimer's disease"]<-"AD"
  fmr2$outcome[fmr2$outcome == "Type 2 diabetes" ]<-"T2D"
  fmr2$outcome[fmr2$outcome == "Chronic kidney disease" ]<-"CKD"
  fmr2$outcome[fmr2$outcome == "Coronary artery disease" ]<-"CAD"
  fmr2$outcome[fmr2$outcome == "ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"ER+ breast cancer"
  fmr2$outcome[fmr2$outcome ==    "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis)"]<-"Breast cancer"
  fmr2$outcome[fmr2$outcome ==   "High grade serous ovarian cancer"  ]<-"High-grade serous ovarian cancer"  
  fmr2$outcome[fmr2$outcome ==  "High grade and low grade serous ovarian cancer"  ]<-"Serous ovarian cancer"
  # 
  fmr2$ppid<-paste(fmr2$UKBPPP.ProteinID,fmr2$outcome)
  fmr2$fdr<-p.adjust(fmr2$pval, "fdr")
  fmr2$z<-fmr2$b/fmr2$se
  fmr2<-fmr2[fmr2$cistrans =="cis",]
  # ppid1<-unique(fmr2$ppid)
  # ppid2<-unique(d$ppid)
  return(fmr2)
}