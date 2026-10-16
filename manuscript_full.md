# Admission AST/ALT ratio and 28-day mortality in sepsis without chronic liver disease

Xiaoye Xu1#, Tielong Chen2#, Pengqin Shu2#, Zhao Liu2#, Fei Ying2, Cheng Bao2, Chao Chen2, Jing Sun2*

# These authors contributed equally to this work and share first authorship.

1 Affiliated Mental Health Center & Hangzhou Seventh People's Hospital, Zhejiang University School of Medicine, Hangzhou, Zhejiang, China.

2 Hangzhou TCM Hospital Affiliated with Zhejiang Chinese Medical University, Hangzhou, Zhejiang, China.

*Corresponding author: Jing Sun, Hangzhou Traditional Chinese Medicine Hospital Affiliated with Zhejiang Chinese Medical University, No. 453 Tiyuchang Road, Xihu District, Hangzhou City, Zhejiang Province, China. Email: sj524258755@163.com.

Short title: AST/ALT ratio and 28-day mortality in sepsis

Article type: Original Research

## Abstract

Background: The AST/ALT ratio tracks mortality in sepsis. Whether the excess risk is concentrated in the first week, and whether the ratio rises because ALT falls, is not established. Methods: Adults with sepsis and no chronic liver disease in MIMIC-IV and eICU were grouped by the first AST/ALT pair within 24 hours, using MIMIC-IV cutpoints of 1.04, 1.44, and 2.03. Follow-up for 28-day death started on the enzyme day. Models included age, sex, organ-failure score, lactate, bilirubin, creatinine, care unit, and alcohol use disorder, with multiple imputation for missing covariates. The eICU outcome was hospital death. Results: Among 16,621 MIMIC-IV patients, 28-day mortality was 19.9%, 22.5%, 27.5%, and 32.5%. Imputed hazard ratios versus the lowest quartile were 1.05 (95% CI 0.96–1.15), 1.19 (1.09–1.30), and 1.45 (1.33–1.59). The highest-quartile hazard ratio was 1.83 (1.62–2.07) in the first week and 1.15 (1.01–1.30) thereafter (ratio 0.63, 0.53–0.75). Higher AST and lower ALT were each associated with death. After comorbidity, weight, albumin, vasopressors, ventilation, and dialysis, the highest-quartile hazard ratio was 1.36 (1.25–1.49). In 12,374 eICU patients, hospital mortality rose from 10.9% to 31.4% (adjusted odds ratio 1.64, 1.43–1.88). Conclusions: A higher admission AST/ALT ratio marks higher early mortality in sepsis without chronic liver disease. The excess hazard is larger in the first week than afterward, and the enzyme pattern is a lower ALT together with a higher AST.

Keywords: sepsis; AST/ALT ratio; De Ritis ratio; mortality; intensive care; MIMIC-IV; eICU

## Introduction

Sepsis remains a leading cause of death in intensive care. Sepsis-3 defines it as life-threatening organ dysfunction caused by a dysregulated host response to infection [1]. The 2021 Surviving Sepsis Campaign guidelines concentrate early care on infection control and restoration of perfusion [2]. Liver involvement sits beside that pathway rather than inside a specific treatment recommendation. The liver can amplify inflammation, alter drug handling, and fail as a target organ [3,4]. Hypoxic, or ischemic, hepatitis is the extreme circulatory form: aminotransferases rise sharply in cardiac, circulatory, or respiratory failure once other causes of hepatic necrosis are set aside [5–7]. In a medical intensive care series, hypoxic hepatitis was frequently accompanied by acute kidney injury and high short-term mortality [8].

AST and ALT do not describe the same injury. AST is present in hepatocytes, skeletal muscle, and myocardium, whereas ALT is more concentrated in the liver [5]. Their ratio, often called the De Ritis ratio, has therefore been used as a severity marker outside ordinary hepatology clinics, including in COVID-19 and in traumatic brain injury [9,10]. In sepsis, the recent evidence is concentrated in narrower syndromes. Chen and colleagues reported that the ratio was associated with mortality in sepsis-associated acute kidney injury [11]. Luo and colleagues then examined trajectories of the ratio in that kidney-injury population [12]. Yu and colleagues studied the ratio together with heparin exposure in sepsis-associated liver injury [13]. Palmowski and colleagues used the ratio to stratify mortality inside sepsis-associated acute liver injury [14]. Separate work has linked lactate trajectories, rather than the aminotransferase ratio, to outcome in the same liver-injury group [15].

The association between the ratio and mortality has also been examined in broader sepsis populations. Choi and colleagues, in a prospective multicenter registry of septic shock, found that an emergency-department AST/ALT ratio of at least 1.84 was associated with 28-day death and added information to the SOFA score [24]. Yang and colleagues, in older patients from MIMIC-IV, found a nonlinear association with hospital death that flattened above a ratio of about 1.8 [25]. Neither analysis asked whether the association changes across the first month, or whether a high ratio reflects a rising AST, a falling ALT, or both. A ratio above 2 is also the classic pattern of alcoholic liver injury, so alcohol use has to be in the model before that pattern is attributed to sepsis.

We therefore aimed to examine whether the association between the admission AST/ALT ratio and mortality varies over the first 28 days, and whether a high ratio reflects elevated AST, reduced ALT, or both. Quartile cutpoints fixed in MIMIC-IV were then applied to hospital mortality in eICU.

## Methods

### Study design and reporting

This was a retrospective two-database cohort study. Reporting followed the STROBE statement and the RECORD statement for routinely collected health data [16,23,26]. The development analysis was specified in MIMIC-IV. No sample-size calculation was performed.

### Data sources

MIMIC-IV version 3.1 contains electronic health-record data from intensive care stays at one United States academic medical center [17,18]. The eICU Collaborative Research Database version 2.0 contains stays from 208 United States hospitals in 2014–2015 [19,20]. Both files were used under the PhysioNet credentialed-data agreement.

### Participants

MIMIC-IV eligibility required a first intensive care stay, age of at least 18 years, and Sepsis-3 status in the derived sepsis table [1]. Chronic liver disease was identified from ICD-9 codes 571.2, 571.5, 571.6, 571.8, and 572.2, and from ICD-10 codes beginning with K70 or K74, or K71.7. Cardiac surgical intensive care stays were retained. In eICU, the first unit stay of a hospitalization was eligible when the diagnosis list contained “sepsis” or “septic.” We excluded cirrhosis, chronic liver disease, and ICD-9 571.2, 571.5, and 571.6. Cardiac, cardiothoracic, and cardiac surgical units were retained.

### Exposure

The exposure was the AST/ALT ratio from the first pair of results within 24 hours after ICU admission in MIMIC-IV or unit admission in eICU. This value is referred to as the admission ratio. Pairs were kept when AST lay between 5 and 5,000 U/L and ALT lay between 5 and 2,000 U/L. MIMIC-IV quartiles were formed at the 25th, 50th, and 75th percentiles and then fixed at 1.0357, 1.4444, and 2.0294. The same cutpoints were applied to the Kaplan–Meier curves and in eICU. Follow-up started on the calendar day of that enzyme pair. A death recorded on an earlier calendar day did not enter the exposed cohort.

### Covariates

The adjusted MIMIC-IV model included age, sex, the first-day SOFA score [21], the first lactate within 24 hours, the logarithm of the first bilirubin, the first creatinine, care unit, and alcohol use disorder. Alcohol use disorder was defined by ICD-10 F10 and ICD-9 303 or 305.0. Care unit was coded as medical, coronary, cardiac surgical, surgical or trauma, or other. A second model entered standardized log AST and standardized log ALT together. A further model added the Charlson comorbidity index, first-day weight, albumin, norepinephrine in the first 24 hours, invasive ventilation in the first 24 hours, and dialysis in the first 24 hours. Creatine kinase was available for a minority of stays and was examined only in that subset. eICU models used the APACHE score in place of SOFA and included age, sex, lactate, bilirubin, creatinine, care-unit group, and an alcohol-related diagnosis. That diagnosis was a string containing “alcohol” or ICD-9 291, 303, or 305.0.

### Outcomes

In MIMIC-IV, the outcome was death within 28 days of the enzyme measurement, using the recorded date of death. Patients alive at 28 days, or dying later, were censored at 28 days. Deaths on the day of measurement were assigned a time of half a day. eICU does not provide the same post-discharge vital status, so the external outcome was death before hospital discharge.

### Statistical analysis

Kaplan–Meier curves and crude rates use every patient with a usable ratio, and the same cutpoints in both databases. Lactate was missing in 31.9% of the 16,621 MIMIC-IV patients and was the main reason a complete-case model kept 10,894 stays. The adjusted hazard ratios come from 30 imputed datasets, because the largest missing fraction was about 32%. Imputation used chained equations with posterior draws. The imputation model contained the analysis covariates, the 28-day death indicator, and the Nelson–Aalen cumulative hazard. Alcohol use disorder came from diagnosis codes and had no missing values, so it was not imputed. The same imputed datasets were used for the 28-day model and for one time-varying Cox model. In that model the quartile coefficients change after day 7 and the other covariates stay constant. Complete-case estimates are shown beside the imputed ones. The creatine-kinase subset and the analysis that removed alcohol use disorder, viral hepatitis, and fatty liver are complete-case analyses. In the creatine-kinase subset, log creatine kinase was a covariate.

Hazard ratios for the second, third, and highest quartiles used the lowest quartile as the reference. A rank Schoenfeld test on the complete-case Cox model rejected proportional hazards for the highest quartile (p = 0.006), and also for SOFA, lactate, age, and bilirubin. The time-varying model lets the quartile coefficients change after day 7 for that reason. A sensitivity analysis let age, SOFA, lactate, and bilirubin change after day 7 as well. Restricted cubic splines of the log ratio described the shape of the association. The plotted curve shows the hazard ratio relative to the median ratio.

Subgroup models contrasted the highest and lowest quartiles and included a product term so that an interaction p value could be reported for sex, age of at least 65 years, lactate of at least 4 mmol/L, and cardiac surgical intensive care. These comparisons were not registered in advance. A sensitivity analysis excluded alcohol use disorder, viral hepatitis (ICD-10 B15–B19, ICD-9 070), and fatty liver (ICD-10 K76.0 or K75.8).

eICU hospital death was analyzed with the same chained-equation imputation used in MIMIC-IV, on all 12,374 patients with a ratio. The complete-case logistic model is reported beside it. Analyses were run in Python 3.14 with pandas 3.0.4, patsy 1.0.2, SciPy 1.18.0, statsmodels 0.14.6, scikit-learn 1.9.0, and lifelines 0.30.0. Figures were drawn in Matplotlib.

### Ethics statement

The collection of MIMIC-IV was approved by the Institutional Review Boards of Beth Israel Deaconess Medical Center and the Massachusetts Institute of Technology, and informed consent was waived because the records are de-identified [17,18]. The eICU Collaborative Research Database is a de-identified dataset released under the HIPAA safe-harbor provision; Privacert certified the re-identification risk (HIPAA Certification no. 1031219-2) [19,20]. This analysis used only those files and involved no direct patient contact.

## Results

### Cohort

After exclusion of chronic liver disease, MIMIC-IV contained 33,401 eligible sepsis stays, including cardiac surgical intensive care. A usable AST and ALT pair was available in 16,621 stays (49.8%). None of those stays had a death recorded on a calendar day before the enzyme measurement. Table 1 gives their characteristics. Median age was 68.9 years (interquartile range 57.4–79.6), 9,403 (56.6%) were men, and the median SOFA score was 6 (4–9). The median ratio was 1.4 (1.0–2.0). From the lowest to the highest quartile, median AST rose from 35 to 60 U/L and median ALT fell from 48 to 21 U/L. Median bilirubin stayed near 0.6–0.7 mg/dL. Crude mortality within 28 days of the enzyme measurement was 826/4,158 (19.9%), 939/4,171 (22.5%), 1,138/4,139 (27.5%), and 1,351/4,153 (32.5%).

Table 2 compares these stays with the 16,780 eligible stays that had no usable pair. Using 28-day death from ICU admission for both groups, mortality was 4,244/16,621 (25.5%) when both enzymes were measured and 2,285/16,780 (13.6%) when they were not. Measured stays had a higher SOFA score (median 6 versus 4) and were less often in the cardiac surgical intensive care unit (7.8% versus 28.4%). Alcohol use disorder was recorded in 1,250 of the 16,621 measured stays (7.5%), viral hepatitis in 677, and fatty liver in 237.

### Unadjusted survival

Kaplan–Meier curves for all 16,621 patients, timed from the enzyme day, separated in quartile order (Figure 2A). The highest quartile left the other curves during the first week.

### Adjusted association

In the imputed model, hazard ratios versus the lowest quartile were 1.05 (95% confidence interval 0.96–1.15), 1.19 (1.09–1.30), and 1.45 (1.33–1.59) (Table 3). The complete-case model gave a highest-quartile hazard ratio of 1.36 (1.23–1.51).

In the time-varying model the highest-quartile hazard ratio was 1.83 (1.62–2.07) during days 0–7 and 1.15 (1.01–1.30) during days 7–28. The later hazard was 0.63 times the first-week hazard (95% confidence interval 0.53–0.75; p < 0.001). When age, SOFA, lactate, and bilirubin were also allowed to change after day 7, the corresponding hazard ratios were 1.70 (1.50–1.93) and 1.24 (1.09–1.41), and their ratio was 0.73 (0.61–0.87).

On the 10,894 complete cases, a joint model entered log AST and log ALT separately. Per standard deviation, the hazard ratio was 1.38 (1.28–1.49) for AST and 0.81 (0.75–0.87) for ALT. The ratio is the contrast between those two directions. Adding the Charlson index, weight, albumin, norepinephrine, invasive ventilation, and dialysis, again with the 30 imputed datasets, gave a highest-quartile hazard ratio of 1.36 (1.25–1.49). Creatine kinase was available in 38% of stays. In that complete-case subset, with log creatine kinase in the model, the highest-quartile hazard ratio was 1.27 (1.03–1.57). Excluding alcohol use disorder, viral hepatitis, and fatty liver in a complete-case analysis yielded a highest-quartile hazard ratio of 1.40 (1.25–1.56; 9,529 patients).

### Subgroups

The association did not differ by sex (interaction hazard ratio 1.03, 95% confidence interval 0.84–1.27; p = 0.75) or by lactate (1.09, 0.85–1.41; p = 0.50). It tended to be weaker at age 65 years or older (0.80, 0.64–1.00; p = 0.054) and in the cardiac surgical intensive care unit (0.61, 0.36–1.03; p = 0.067) (Table 4).

### Replication in eICU

eICU contained 22,405 eligible first stays after the liver-disease exclusion, and 12,374 had a usable ratio. Hospital mortality in these 12,374 patients was 424/3,881 (10.9%), 463/3,033 (15.3%), 538/2,715 (19.8%), and 862/2,745 (31.4%) (Figure 2B). Imputed odds ratios versus the lowest quartile were 0.87 (0.76–1.01), 1.05 (0.91–1.21), and 1.64 (1.43–1.88). In the complete-case model (7,092 patients, 1,565 hospital deaths) the highest-quartile odds ratio was 1.96 (1.64–2.36).

## Discussion

### Principal findings

In this two-database cohort of patients with sepsis but without chronic liver disease, a higher admission AST/ALT ratio was associated with higher 28-day mortality. The findings extend earlier work in two ways [24,25]. First, the excess risk was larger early: the hazard ratio for the highest quartile was 1.83 during the first week and 1.15 thereafter. Second, the ratio carries information from both enzymes. Across quartiles, median ALT fell from 48 to 21 U/L while median AST rose from 35 to 60 U/L, and a higher AST and a lower ALT were each associated with death.

Crude hospital mortality in eICU rose from 10.9% to 31.4% across the same cutpoints. After adjustment, the association was confined to the highest quartile (odds ratio 1.64). In MIMIC-IV the highest-quartile hazard ratio remained 1.36 after comorbidity, weight, albumin, vasopressors, ventilation, and dialysis were added, and 1.40 after alcohol use disorder, viral hepatitis, and fatty liver were removed.

The early concentration of risk is consistent with the clinical course of the first ICU days, when circulatory failure, hepatic hypoperfusion, and AST release from muscle or myocardium may contribute most directly to death [5,8].

### What the enzyme pattern shows

The enzyme profile in the highest quartile was modest: median AST was 60 U/L and bilirubin was unchanged, far below the levels typical of hypoxic hepatitis [5–7,22]. The ratio was therefore carried by a low ALT as well as by a small rise in AST.

A low ALT can mark reduced physiological reserve. In older adults it tracks frailty, disability, and sarcopenia, and low ALT predicts higher long-term mortality [27,28]. The association here was not stronger at age 65 years or older (interaction hazard ratio 0.80, 0.64–1.00), so old-age frailty alone does not account for it. Albumin and weight, together with the Charlson index and markers of shock and respiratory failure, attenuated the highest-quartile hazard ratio from 1.45 to 1.36. ALT activity also depends on pyridoxal 5′-phosphate, and deficiency of that cofactor can lower the measured value [29]. We did not have vitamin B6 levels, so this remains a possible contributor rather than a demonstrated pathway. Where creatine kinase was available, adding it yielded a highest-quartile hazard ratio of 1.27. Excluding recorded alcohol use disorder, viral hepatitis, and fatty liver yielded a hazard ratio of 1.40.

### Relation to previous studies

Choi and colleagues reported an adjusted hazard ratio of 1.46 (95% CI 1.30–1.64) for an emergency-department ratio of at least 1.84 in septic shock, and a gain in discrimination when the ratio was added to SOFA [24]. Yang and colleagues reported a rise in hospital mortality up to a ratio of about 1.8 and little further rise above it, in MIMIC-IV patients aged 65 years or older [25]. Chen, Luo, Yu, and Palmowski studied the ratio inside sepsis-associated kidney or liver injury [11–14]. The present analysis uses an unselected sepsis cohort, fixes the cutpoints before looking at eICU, and shows that the mortality difference is larger in the first week and is carried by both a higher AST and a lower ALT.

### Implications

A ratio above 2.03 identified the quartile with 32.5% mortality at 28 days in MIMIC-IV and 31.4% crude hospital mortality in eICU. The ratio may help identify patients at increased risk during the first week of intensive care, when closer monitoring may be warranted.

### Limitations

Enzymes were measured in half of eligible stays, and those patients were more severely ill than patients without a pair. The estimates apply to patients in whom clinicians ordered both enzymes. ICD codes can miss alcohol use disorder and chronic liver disease. eICU identified chronic liver disease in 1.9% of first sepsis stays, against 11.0% in MIMIC-IV, so some chronic liver disease likely remains in the eICU cohort and could inflate the apparent association. First-day SOFA, norepinephrine, ventilation, and dialysis can be recorded after the enzyme measurement. If these variables lie on the pathway from the enzyme pattern to death, adjusting for them may attenuate the estimate; if they act as colliders, the direction of bias is uncertain. MIMIC-IV is a single center. Both databases are from the United States, the eICU stays date from 2014–2015, and the eICU outcome is hospital death rather than 28-day death. The cardiac surgical comparison was too imprecise to support a separate conclusion.

## Conclusion

In sepsis without chronic liver disease, a higher admission AST/ALT ratio was associated with higher 28-day mortality. The excess hazard was larger in the first week and was still present afterward. The ratio reflected a lower ALT as well as a higher AST, and it remained after alcohol use disorder, comorbidity, and shock were taken into account. In eICU, the highest quartile was likewise associated with higher hospital mortality (adjusted odds ratio 1.64), although the intermediate quartiles did not differ after adjustment.

## Acknowledgments

The authors thank the teams responsible for developing and maintaining MIMIC-IV, eICU-CRD, and PhysioNet.

## Data availability statement

MIMIC-IV v3.1 and eICU-CRD v2.0 are available through PhysioNet to credentialed users who complete the required training and data-use procedures. The analysis code is available at https://github.com/shupengqin/ast-alt-sepsis-mortality. Patient-level records are not included.

## Funding

The authors received no specific funding for this study.

## Conflict of interest

The authors declare that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.

## Author contributions

Xiaoye Xu, Tielong Chen, Pengqin Shu, and Zhao Liu contributed equally. Xiaoye Xu and Tielong Chen built the cohort and ran the analysis. Pengqin Shu interpreted the clinical findings and revised the manuscript. Zhao Liu prepared the figures. Fei Ying and Cheng Bao checked the analysis. Chao Chen contributed to data curation and clinical interpretation. Jing Sun supervised the study and led the writing. All authors read and approved the final manuscript.

## Figure legends

Figure 1. Cohort construction. Cardiac surgical and cardiothoracic units are retained. Quartile membership uses the fixed MIMIC-IV cutpoints 1.0357, 1.4444, and 2.0294 in both databases. The analysis counts are the patients with a usable ratio: 16,621 in MIMIC-IV and 12,374 in eICU-CRD.

Figure 2. Association of the admission AST/ALT ratio with short-term death, timed from the enzyme measurement in MIMIC-IV. (A) Unadjusted 28-day survival in all 16,621 MIMIC-IV patients with a ratio. (B) Crude death by quartile in those 16,621 patients and in all 12,374 eICU patients with a ratio. MIMIC-IV is 28-day death. eICU is hospital death. (C) Adjusted associations. MIMIC-IV estimates are multiply imputed hazard ratios. eICU estimates are multiply imputed hospital-death odds ratios. (D) Complete-case spline in MIMIC-IV, adjusted for age, sex, SOFA, lactate, bilirubin, creatinine, and care unit, relative to the median ratio.

Figure 3. (A) Crude 28-day mortality by lactate. (B) Interaction hazard ratios for the highest-versus-lowest contrast. A value below 1 means a smaller association in the named group. The cardiac surgical interval includes 1.

## References

1. Singer M, Deutschman CS, Seymour CW, Shankar-Hari M, Annane D, Bauer M, et al. The Third International Consensus Definitions for Sepsis and Septic Shock (Sepsis-3). JAMA. 2016;315:801-810. doi:10.1001/jama.2016.0287.

2. Evans L, Rhodes A, Alhazzani W, Antonelli M, Coopersmith CM, French C, et al. Surviving sepsis campaign: international guidelines for management of sepsis and septic shock 2021. Intensive Care Med. 2021;47:1181-1247. doi:10.1007/s00134-021-06506-y.

3. Strnad P, Tacke F, Koch A, Trautwein C. Liver - guardian, modifier and target of sepsis. Nat Rev Gastroenterol Hepatol. 2017;14:55-66. doi:10.1038/nrgastro.2016.168.

4. Nesseler N, Launey Y, Aninat C, Morel F, Mallédant Y, Seguin P. Clinical review: the liver in sepsis. Crit Care. 2012;16:235. doi:10.1186/cc11381.

5. Henrion J. Hypoxic hepatitis. Liver Int. 2012;32:1039-1052. doi:10.1111/j.1478-3231.2011.02655.x.

6. Henrion J, Schapira M, Luwaert R, Colin L, Delannoy A, Heller FR. Hypoxic hepatitis: clinical and hemodynamic study in 142 consecutive cases. Medicine (Baltimore). 2003;82:392-406. doi:10.1097/01.md.0000101573.54295.bd.

7. Smith JE, Rockey DC. Update on ischemic hepatitis. Curr Opin Gastroenterol. 2024;40:143-147. doi:10.1097/MOG.0000000000001017.

8. Drolz A, Horvatits T, Roedl K, Rutter K, Staufer K, Haider DG, et al. Outcome and features of acute kidney injury complicating hypoxic hepatitis at the medical intensive care unit. Ann Intensive Care. 2016;6:61. doi:10.1186/s13613-016-0162-4.

9. Mangoni AA, Zinellu A. An updated systematic review and meta-analysis of the association between the De Ritis ratio and disease severity and mortality in patients with COVID-19. Life (Basel). 2023;13:1324. doi:10.3390/life13061324.

10. Lin YC, Tsai CH, Su WT, Hsu SY, Hsieh CH, Lin CH. De Ritis ratio as a prognostic marker for mortality in moderate-to-severe traumatic brain injury: a propensity score-matched analysis. Diagnostics (Basel). 2025;15:2416. doi:10.3390/diagnostics15192416.

11. Chen S, Ye X, Guo L, Ma X, Li Y, Wu L, et al. Association of De Ritis ratio with mortality in ICU patients with sepsis-associated acute kidney injury: a retrospective multicenter cohort study. BMC Infect Dis. 2026 Feb 19. doi:10.1186/s12879-026-12569-3.

12. Luo X, Huang G, Zeng W, Lai N, Zhang M, Zhu D, et al. De Ritis ratio trajectories and mortality risk in sepsis-associated acute kidney injury. J Surg Res. 2026;321:96-106. doi:10.1016/j.jss.2026.02.009.

13. Yu B, Li D, Luan Y, Rajagopalan R, Guan Y, Xiang Y, et al. De Ritis ratio and heparin therapy in sepsis-associated liver injury: a multicenter cohort study. Thromb Haemost. 2026 Apr 23. doi:10.1055/a-2852-9288. Epub ahead of print.

14. Palmowski L, Westhus B, Witowski A, Nowak H, Traut I, Canbay A, et al. Subphenotypes and the De Ritis ratio for mortality risk stratification in sepsis-associated acute liver injury: a retrospective cohort study. EClinicalMedicine. 2025;82:103173. doi:10.1016/j.eclinm.2025.103173.

15. Wei H, Zhang J, Yu S, Zhang C, Weng Y, Shao M, et al. Association between lactate trajectories and clinical outcomes in critically ill patients with sepsis-associated liver injury. J Clin Med. 2026;15:6837. doi:10.3390/jcm15176837.

16. von Elm E, Altman DG, Egger M, Pocock SJ, Gøtzsche PC, Vandenbroucke JP, et al. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement: guidelines for reporting observational studies. Lancet. 2007;370:1453-1457. doi:10.1016/S0140-6736(07)61602-X.

17. Johnson AEW, Bulgarelli L, Shen L, Gayles A, Shammout A, Horng S, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10:1. doi:10.1038/s41597-022-01899-x.

18. Johnson A, Bulgarelli L, Pollard T, Gow B, Moody B, Horng S, et al. MIMIC-IV (version 3.1) [dataset]. PhysioNet; 2024. doi:10.13026/kpb9-mt58.

19. Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU Collaborative Research Database, a freely available multi-center database for critical care research. Sci Data. 2018;5:180178. doi:10.1038/sdata.2018.178.

20. Pollard T, Johnson A, Raffa J, Celi LA, Badawi O, Mark R. eICU Collaborative Research Database (version 2.0) [dataset]. PhysioNet; 2019. doi:10.13026/C2WM1R.

21. Vincent JL, Moreno R, Takala J, Willatts S, De Mendonça A, Bruining H, et al. The SOFA (Sepsis-related Organ Failure Assessment) score to describe organ dysfunction/failure. Intensive Care Med. 1996;22:707-710. doi:10.1007/BF01709751.

22. Drolz A, Horvatits T, Michl B, Roedl K, Schellongowski P, Holzinger U, et al. Statin therapy is associated with reduced incidence of hypoxic hepatitis in critically ill patients. J Hepatol. 2014;60:1187-1193. doi:10.1016/j.jhep.2014.01.019.

23. Vandenbroucke JP, von Elm E, Altman DG, Gøtzsche PC, Mulrow CD, Pocock SJ, et al. Strengthening the Reporting of Observational Studies in Epidemiology (STROBE): explanation and elaboration. Epidemiology. 2007;18:805-835. doi:10.1097/EDE.0b013e3181577511.

24. Choi S, Nah S, Suh GJ, Choi SH, Chung SP, Kim WY, et al. Prognostic value of the AST/ALT ratio in patients with septic shock: a prospective, multicenter, registry-based observational study. Diagnostics (Basel). 2025;15:1773. doi:10.3390/diagnostics15141773.

25. Yang H, Meng L, Xin S, Chang C, Zhao X, Guo B, et al. Non-linear association between aspartate aminotransferase to alanine aminotransferase ratio and mortality in critically ill older patients: a retrospective cohort study. PLoS One. 2023;18:e0293749. doi:10.1371/journal.pone.0293749.

26. Benchimol EI, Smeeth L, Guttmann A, Harron K, Moher D, Petersen I, et al. The REporting of studies Conducted using Observational Routinely-collected health Data (RECORD) statement. PLoS Med. 2015;12:e1001885. doi:10.1371/journal.pmed.1001885.

27. Vespasiani-Gentilucci U, De Vincentis A, Ferrucci L, Bandinelli S, Antonelli Incalzi R, Picardi A. Low alanine aminotransferase levels in the elderly population: frailty, disability, sarcopenia, and reduced survival. J Gerontol A Biol Sci Med Sci. 2018;73:925-930. doi:10.1093/gerona/glx126.

28. Ramaty E, Maor E, Peltz-Sinvani N, Brom A, Grinfeld A, Kivity S, et al. Low ALT blood levels predict long-term all-cause mortality among adults. A historical prospective cohort study. Eur J Intern Med. 2014;25:919-921. doi:10.1016/j.ejim.2014.10.019.

29. Youssef EM, Wu GY. Subnormal serum liver enzyme levels: a review of pathophysiology and clinical significance. J Clin Transl Hepatol. 2024;12:428-435. doi:10.14218/JCTH.2023.00446.
