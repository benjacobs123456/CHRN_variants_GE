# Code used to replicate findings from Briggs 2020 (MSJ) in UKB

### The original paper is [here](https://journals.sagepub.com/doi/full/10.1177/1352458520958361#bibr14-1352458520958361)

#### First, we prepared the cohort, covariate files, applied individual QC, and defined smoking status in R
````r
library(readr)
library(dplyr)
setwd("/data/Wolfson-UKBB-Dobson/smoking_gwis")
df = read_tsv("/data/Wolfson-UKBB-Dobson/ukb_pheno_911/ukb_pheno_final_MS_1301")

# remove any dups
df = df %>% distinct(EID,.keep_all=TRUE)
#remove ms cases with no age at dx / diagnosed before 20

df = df %>% filter(!(MS_status==1 & (age_at_ms_diagnosis<30 | is.na(age_at_ms_diagnosis))))
# select only unrelated white individuals
kin = read_table2("/data/Wolfson-UKBB-Dobson/helper_progs_and_key/ukb43101_rel_s488282.dat")
exclusion = kin %>% filter(Kinship>0.0884) %>% select(ID1) %>% rename("EID"="ID1") 
df = df %>% 
  mutate("FID"=EID) %>% 
  rename("IID"=EID) %>% 
  mutate("FID"=IID)  
df = df %>% distinct(IID,.keep_all = TRUE)
# relatedness and ethnic groups
df = df %>% filter(!IID %in% exclusion$EID) %>% filter(`Genetic ethnic grouping.0.0`=="Caucasian")

df = df %>% select(FID,IID,Sex.0.0,contains("enetic"),MS_status,"Smoking status.0.0","Age at recruitment.0.0",contains("Age started smoking in"))

#make young starters
young= df %>% filter(`Smoking status.0.0`=="Previous" &`Age started smoking in former smokers.0.0` < 20) %>% mutate("young_smoker"="2")
young2= df %>% filter(`Smoking status.0.0`=="Current" &`Age started smoking in current smokers.0.0` < 20) %>% mutate("young_smoker"="2")
nyoung= df %>% filter(!IID %in% young$IID) %>% filter(!IID %in% young2$IID) %>% mutate("young_smoker"="1")
df = bind_rows(young,young2,nyoung)


# whole unmatched cohort
# write ms pheno file
ms = df %>% select(FID,IID,MS_status) %>% mutate(MS_status = ifelse(MS_status==1,2,1))
write_tsv(ms,"whole_cohort_MS_pheno.txt")


# write individuals to keep
indivs = df %>% select(FID,IID) 
write_tsv(indivs,"whole_cohort_indivs_to_include",col_names=FALSE)

#write covars
covar=df %>% select(FID, IID,`Age at recruitment.0.0`,Sex.0.0,young_smoker,contains("Genetic principal"),-`Used in genetic principal components.0.0`)

colnames(covar)=c("FID","IID","Age","Sex","young_smoker",paste0("pc",c(1:40)))
covar$Sex = recode(covar$Sex,"Male"="1","Female"="2")

write_tsv(covar,"whole_cohort_covars.txt")
````

#### Next, we ran GxE interaction models in PLINKv2 for CHRN7A
        plink2 \
          --chr 15 \
          --covar /data/Wolfson-UKBB-Dobson/smoking_gwis/whole_cohort_covars.txt \
          --covar-name Age Sex young_smoker pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 \
          --from-bp 32312691 \
          --geno 0.1 \
          --glm interaction \
          --hwe 1e-5 \
          --maf 0.01 \
          --mind 0.1 \
          --out chnra7 \
          --parameters 1, 3-14, 17 \
          --pfile /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files/chr_15 \
          --pheno /data/Wolfson-UKBB-Dobson/smoking_gwis/whole_cohort_MS_pheno.txt \
          --to-bp 32474722

#### And for CHRN9A
        plink2 \
          --chr 4 \
          --covar /data/Wolfson-UKBB-Dobson/smoking_gwis/whole_cohort_covars.txt \
          --covar-name Age Sex young_smoker pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 \
          --from-bp 40327346 \
          --geno 0.1 \
          --glm interaction \
          --hwe 1e-5 \
          --maf 0.01 \
          --mind 0.1 \
          --out chnra9 \
          --parameters 1, 3-14, 17 \
          --pfile /data/Wolfson-UKBB-Dobson/imputed_ukb_genotypes/plink2_files/chr_4 \
          --pheno /data/Wolfson-UKBB-Dobson/smoking_gwis/whole_cohort_MS_pheno.txt \
          --to-bp 40367234

#### Now we have a good look in R

        # read in libs
        library(tidyverse)
        library(qqman)

        # read in GxE results
        chrn7 = read_table2("/data/Wolfson-UKBB-Dobson/smoking_gwis/chnra7.MS_status.glm.logistic.hybrid")
        chrn9 = read_table2("/data/Wolfson-UKBB-Dobson/smoking_gwis/chnra9.MS_status.glm.logistic.hybrid")

        # filter to GxE results and count how man ps < 0.05
        chrn7 = chrn7 %>% filter(TEST=="ADDxyoung_smoker")
        chrn7 %>% filter(P<0.05) %>% count()
        chrn9 = chrn9 %>% filter(TEST=="ADDxyoung_smoker")
        chrn9 %>% filter(P<0.05) %>% count()

        # save snps
        chrn7_snps = chrn7 %>% select(ID,P) %>% rename(SNP=ID)
        chrn9_snps = chrn9 %>% select(ID,P) %>% rename(SNP=ID)
        write_tsv(chrn7_snps,"chrn7_snps")
        write_tsv(chrn9_snps,"chrn9_snps")

        # clump
        system("module load plink;plink --bfile /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr4 --clump chrn9_snps --clump-p1 1 --clump-p2 1 --out chrn9" )
        system("module load plink;plink --bfile /data/Wolfson-UKBB-Dobson/1kg_reference/filtered_chr15 --clump chrn7_snps --clump-p1 1 --clump-p2 1 --out chrn7" )


        # define p threshold
        p_threshold = 0.05/(63+27)
        chrn9=chrn9 %>% mutate(below_threshold = ifelse(P<p_threshold,TRUE,FALSE))
        chrn7=chrn7 %>% mutate(below_threshold = ifelse(P<p_threshold,TRUE,FALSE))

        # read in briggs snps
        briggs_snps = read_csv("briggs_st2.csv")
        briggs_snps = briggs_snps %>% select(rsID, `Bootstrapped p-value`, `Interaction OR (Bootstrapped bias-corrected 95% CI)`,A1)

        # filter to overlapping snps
        briggs_snps =  briggs_snps %>% filter(rsID %in% chrn7_snps$SNP | rsID %in% chrn9_snps$SNP)
        combo =  briggs_snps %>% rename(ID = rsID) %>% left_join(bind_rows(chrn7,chrn9),by="ID")

        # find p proportion
        combo = combo %>% mutate("p_prop" = `Bootstrapped p-value`/P)

        # get ORs from briggs
        combo = combo %>% separate(`Interaction OR (Bootstrapped bias-corrected 95% CI)`,sep="\\(",into=c("OR_briggs"," "))

        # ensure effect alleles are aligned
        combo = combo %>% filter(A1.x==A1.y)


        combo$OR_briggs =  as.numeric(combo$OR_briggs)

        # get betas and compare
        combo$beta_briggs = log(combo$OR_briggs)
        combo$beta_ukb = log(combo$OR)
        combo = combo %>% mutate("beta_prop" = beta_briggs/beta_ukb)




