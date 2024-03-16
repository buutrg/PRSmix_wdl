# Introduction
The PRSmix package can be used to benchmark PRSs and compute the linear combination of trait-specific scores (PRSmix) and cross-trait scores (PRSmix+) to improve prediction accuracy.


NOTE: 
- You can run all steps or run each of them individually (set a_runStep1 or a_runStep2 to true or false).
- Please make sure ALL optional FILES for all steps are provided even if you do not run that step (you can use a dummy file)

# Manual
The pipeline contains 3 steps: 
- Harmonizing SNP effects: s1_harmonize_SNPeffects
- Compute PRS: s2_computePRS
- Combine PRS: s3_combine_PRS

### Step 1: Harmonizing SNP effects
- This step makes sure all set of SNP weights aligned to the same risk allele. (**Highly important**)
- For more detail, please read more at [PRSmix github part 1](https://github.com/buutrg/PRSmix?tab=readme-ov-file#1-harmonize-per-allele-effect-sizes-to-the-effects-of-alternative-allele-in-the-target-cohort).
- Output: A file with Harmonized SNP effect.

### Step 2: Compute PRS
- This step computes PRS with the SNP effects harmonized in step 1.
- For more detail, please read more at [PRSmix github part 2](https://github.com/buutrg/PRSmix?tab=readme-ov-file#2-compute-prss-for-all-scores).
- Output: A file with PRS computed

### Step 3: Combine PRS
- This step performs PRS combination.
- For more detail, please read more at [PRSmix github part 3](https://github.com/buutrg/PRSmix?tab=readme-ov-file#3-perform-linear-combination-trait-specific-prsmix-and-cross-trait-prsmix).
- Output: A tarball with:
	- The dataframe of training and testing sample split from the main dataframe: `_train_df.txt` and `_test_df.txt` 
	- The prediction accuracy in the training set for all PRSs: `_train_allPRS.txt` 
	- The prediction accuracy in the testing set for trait-specific PRSs: `_test_allPRS.txt` 
	- The prediction accuracy of PRSmix and PRSmix+: `_test_summary_traitPRS_withPRSmix.txt` and `_test_summary_traitPRS_withPRSmixPlus.txt` 

# References

Truong, B., Hull, L. E., Ruan, Y., Huang, Q. Q., Hornsby, W., Martin, H. C., van Heel, D. A., Wang, Y., Martin, A. R., Lee, H. and Natarajan, P. 2023. "Integrative Polygenic Risk Score Improves The Prediction Accuracy Of Complex Traits And Diseases". *doi.org/10.1016/j.xgen.2024.100523*.

# Contact information
This workspace is intended to be used for development of Consortium methods for PRSmix.

Please contact Buu Truong (btruong@broadinstitute.org) or Pradeep Natarajan (pradeep@broadinstitute.org) if you have any queries.