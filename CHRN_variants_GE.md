# Code used to replicate findings from Briggs 2020 (MSJ) in UKB

### The original paper is [here](https://journals.sagepub.com/doi/full/10.1177/1352458520958361#bibr14-1352458520958361)

#### First, we ran GxE interaction models in PLINKv2 for CHRN7A
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
