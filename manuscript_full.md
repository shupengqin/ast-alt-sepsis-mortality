# Admission AST/ALT ratio and 28-day mortality in sepsis: a two-database cohort study

Xiaoye Xu1#, Tielong Chen2#, Pengqin Shu2#, Zhao Liu2#, Fei Ying2, Cheng Bao2, Chao Chen2, Jing Sun2*

# These authors contributed equally to this work and share first authorship.

1 Affiliated Mental Health Center & Hangzhou Seventh People's Hospital, Zhejiang University School of Medicine, Hangzhou, Zhejiang, China.

2 Hangzhou TCM Hospital Affiliated to Zhejiang Chinese Medical University, Hangzhou, Zhejiang, China.

*Corresponding author: Jing Sun, Hangzhou TCM Hospital Affiliated to Zhejiang Chinese Medical University, No. 453 Tiyuchang Road, Xihu District, Hangzhou City, Zhejiang Province, China. Email: sj524258755@163.com.

Short title: AST/ALT ratio and 28-day mortality in sepsis

Article type: Original Research

## Abstract

Background: The aspartate aminotransferase (AST) to alanine aminotransferase (ALT) ratio is associated with mortality in sepsis. However, whether this association is stronger early, and whether the ratio fully captures the prognostic information contained in AST and ALT, remains unclear. Methods: Adults with sepsis in the Medical Information Mart for Intensive Care IV (MIMIC-IV) and the eICU Collaborative Research Database (eICU) were included after exclusion of a prespecified set of chronic liver-disease codes. Patients were grouped by the first AST/ALT pair within 24 hours, using MIMIC-IV cutpoints of 1.04, 1.44, and 2.03. Follow-up for 28-day death started on the day of enzyme measurement. Models included age, sex, the Sequential Organ Failure Assessment (SOFA) score, lactate, bilirubin, creatinine, care unit, and alcohol use disorder, with multiple imputation for missing covariates. The eICU outcome was hospital death. Results: Among 16,616 MIMIC-IV patients, 28-day mortality was 19.8%, 22.5%, 27.4%, and 32.5%. Imputed hazard ratios versus the lowest quartile were 1.04 (95% CI 0.95–1.14), 1.17 (1.07–1.29), and 1.43 (1.31–1.57). The highest-quartile hazard ratio was 1.81 (1.60–2.05) in the first week and 1.13 (0.99–1.28) from day 7 to day 28 (ratio 0.62, 0.52–0.74). In the imputed model, each standard deviation of the log ratio had a hazard ratio of 1.16 (1.13–1.20) and each standard deviation of their log mean had a hazard ratio of 1.09 (1.06–1.13). In 12,374 eICU patients, hospital mortality rose from 10.9% to 31.4% (adjusted odds ratio for the highest quartile 1.95, 1.69–2.26). Conclusions: Among patients with sepsis who had both enzymes measured, a higher admission AST/ALT ratio was associated with higher short-term mortality, most strongly in the first week. The ratio did not capture all of the prognostic information in the two enzymes. The association was replicated for hospital mortality in eICU, but the ratio added little discrimination beyond severity of illness.

Keywords: sepsis; AST/ALT ratio; De Ritis ratio; mortality; intensive care; MIMIC-IV; eICU

## Introduction

Sepsis remains a leading cause of death in intensive care. Sepsis-3 defines it as life-threatening organ dysfunction caused by a dysregulated host response to infection [1]. The 2021 Surviving Sepsis Campaign guidelines concentrate early care on infection control and restoration of perfusion [2]. Liver dysfunction, however, is not the target of any specific treatment recommendation. The liver can amplify inflammation, alter drug handling, and fail as a target organ [3,4]. Hypoxic, or ischemic, hepatitis is the extreme circulatory form: aminotransferases rise sharply in cardiac, circulatory, or respiratory failure once other causes of hepatic necrosis are set aside [5–7]. In a medical intensive care series, hypoxic hepatitis was frequently accompanied by acute kidney injury and high short-term mortality [8].

Aspartate aminotransferase (AST) and alanine aminotransferase (ALT) reflect distinct tissue distributions and pathophysiological processes. AST is present in hepatocytes, skeletal muscle, and myocardium, whereas ALT is more concentrated in the liver [5]. Their ratio, often called the De Ritis ratio, has therefore been used as a severity marker beyond hepatology, including in COVID-19 and in traumatic brain injury [9,10]. In sepsis, most studies have examined the ratio within narrower syndromes, such as sepsis-associated acute kidney injury [11,12] and sepsis-associated liver injury [13,14].

The association between the ratio and mortality has also been examined beyond a single syndrome. Wang and Xu studied adult intensive care stays in MIMIC-IV, with eICU as a sensitivity cohort, and found a nonlinear association with 28-day death [15]. That cohort was not limited to sepsis. Choi and colleagues, in a prospective multicenter registry of septic shock, found that an emergency-department AST/ALT ratio of at least 1.84 was associated with 28-day death and added information to the SOFA score [16]. Yang and colleagues, in older patients from MIMIC-IV, found a nonlinear association with hospital death that flattened above a ratio of about 1.8 [17]. Those studies describe the overall association and its shape. They do not test whether the association in sepsis is stronger early than later, and they do not test whether the ratio contains the prognostic information in both enzymes. A ratio above 2 is also the classic pattern of alcoholic liver injury, so alcohol use must be accounted for before the pattern is attributed to sepsis.

We therefore aimed to examine whether the association between the admission AST/ALT ratio and 28-day mortality is stronger early than later, and whether the ratio fully captures the prognostic information contained in AST and ALT. Quartile cutpoints fixed in MIMIC-IV were applied to hospital mortality in eICU.

## Methods

### Study design and reporting

This was a retrospective cohort study in two databases. Reporting followed the STROBE statement and the RECORD statement for routinely collected health data [18–20]. The time-varying analysis was specified in MIMIC-IV. eICU was used to repeat the association with hospital death. No sample-size calculation was performed.

### Data sources

The Medical Information Mart for Intensive Care IV (MIMIC-IV) version 3.1 contains electronic health-record data from intensive care stays at one United States academic medical center [21,22]. The eICU Collaborative Research Database (eICU) version 2.0 contains stays from 208 United States hospitals in 2014–2015 [23,24]. Both databases were accessed under the PhysioNet credentialed-data agreement.

### Participants

MIMIC-IV eligibility required a first intensive care stay, age of at least 18 years, and Sepsis-3 status in the derived sepsis table [1]. The liver exclusion used ICD-9 codes 571.2, 571.5, 571.6, 571.8, and 572.2, and ICD-10 codes beginning with K70 or K74, or K71.7. Viral hepatitis and fatty liver were not in this list. Cardiac surgical intensive care stays were retained. In eICU, the first unit stay of a hospitalization was eligible when the diagnosis list contained “sepsis” or “septic.” We excluded records coded as cirrhosis or chronic liver disease, and ICD-9 571.2, 571.5, and 571.6. Cardiac, cardiothoracic, and cardiac surgical units were retained.

### Exposure

The exposure was the AST/ALT ratio from the first pair of results within 24 hours after ICU admission in MIMIC-IV or unit admission in eICU. This value is referred to as the admission ratio. In MIMIC-IV the two enzymes come from the same laboratory record. In eICU the first AST and the first ALT had the same time offset in 14,773 of 14,803 stays with both results (99.8%). Pairs were kept when AST lay between 5 and 5,000 U/L and ALT lay between 5 and 2,000 U/L. MIMIC-IV quartiles were formed at the 25th, 50th, and 75th percentiles and then fixed at 1.0357, 1.4444, and 2.0294. The same cutpoints were applied to the Kaplan–Meier curves and in eICU. Follow-up started on the calendar day of that enzyme pair. A death recorded on an earlier calendar day did not enter the exposed cohort. The time from ICU admission to the enzyme chart time was summarized. For deaths in hospital, that chart time was compared with the recorded death time. Five stays in which the chart time followed the death time were removed before the models were fit. Deaths after discharge have a date and no clock time.

### Covariates

The adjusted MIMIC-IV model included age, sex, the first-day Sequential Organ Failure Assessment (SOFA) score [25], the first lactate within 24 hours, the logarithm of the first bilirubin, the first creatinine, care unit, and alcohol use disorder. Alcohol use disorder was defined by ICD-10 F10 and ICD-9 303 or 305.0. Care unit was coded as medical, coronary, cardiac surgical, surgical or trauma, or other. A second model entered standardized log AST and standardized log ALT together. The same information was also entered as the standardized log ratio and the standardized mean of log AST and log ALT. These two forms are the same model with different parameters. The coefficient on the mean level is zero only when log AST and log ALT have coefficients of equal size and opposite sign, which is the case in which the ratio alone carries the association. A further model added the Charlson comorbidity index, first-day weight, albumin, norepinephrine in the first 24 hours, invasive ventilation in the first 24 hours, and dialysis in the first 24 hours. Creatine kinase was available for a minority of stays and was examined only in that subset. eICU models used the Acute Physiology and Chronic Health Evaluation (APACHE) score in place of SOFA and included age, sex, lactate, bilirubin, creatinine, care-unit group, and an alcohol-related diagnosis. That diagnosis was a string containing “alcohol” or ICD-9 291, 303, or 305.0.

### Outcomes

In MIMIC-IV, the outcome was death within 28 days of the enzyme measurement, using the recorded date of death. That date includes deaths matched from state records after hospital discharge. Patients alive at 28 days, or dying later, were censored at 28 days. Deaths on the day of measurement were assigned a time of half a day. eICU does not provide the same post-discharge vital status, so the replication outcome was death before hospital discharge. The same hospital-death endpoint was also fit in MIMIC-IV, so the two databases could be compared on one definition.

### Statistical analysis

Kaplan–Meier curves and crude rates use every patient with a valid ratio, and the same cutpoints in both databases. Lactate was missing in 31.9% of the 16,616 MIMIC-IV patients and was the main reason a complete-case model kept 10,890 stays. The adjusted hazard ratios come from 30 imputed datasets, because the largest missing fraction was about 32%. Imputation used chained equations with posterior draws. The imputation model contained the quartile indicators, the analysis covariates, the 28-day death indicator, and the Nelson–Aalen cumulative hazard. The quartiles were complete, so they entered the imputation as predictors and were not themselves imputed. Alcohol use disorder came from diagnosis codes and had no missing values, so it was not imputed. The same imputed datasets were used for the 28-day model and for one time-varying Cox model. In that model the quartile coefficients change after day 7 and the other covariates stay constant. Complete-case estimates are shown beside the imputed ones. The creatine-kinase subset and the analysis that removed alcohol use disorder, viral hepatitis, and fatty liver are complete-case analyses. In the creatine-kinase subset, log creatine kinase was a covariate.

Hazard ratios for the second, third, and highest quartiles used the lowest quartile as the reference. A rank Schoenfeld test on the complete-case Cox model rejected proportional hazards for the highest quartile (p = 0.006), and also for SOFA, lactate, age, and bilirubin. The time-varying model lets the quartile coefficients change after day 7 for that reason. A sensitivity analysis let age, SOFA, lactate, and bilirubin change after day 7 as well. Restricted cubic splines of the log ratio described the shape of the association. The plotted curve shows the hazard ratio relative to the median ratio.

Subgroup models contrasted the highest and lowest quartiles and included a product term so that an interaction p value could be reported for sex, age of at least 65 years, lactate of at least 4 mmol/L, and cardiac surgical intensive care. These comparisons were not registered in advance. A sensitivity analysis excluded alcohol use disorder, viral hepatitis (ICD-10 B15–B19, ICD-9 070), and fatty liver (ICD-10 K76.0 or K75.8).

eICU hospital death was analyzed with the same chained-equation imputation, including the quartile indicators, on all 12,374 patients with a ratio. Logistic models include an intercept. The reported standard errors are robust and clustered by hospital. The complete-case logistic model is reported beside the imputed one. A sensitivity analysis excluded patients who had already died within 24 hours after ICU admission and refit the day-7 split in that group. Follow-up was still counted from the day of enzyme measurement. A further sensitivity analysis removed patients who died on the calendar day of measurement and refit the same day-7 split, including the later-to-early ratio. The time-varying model was repeated with the split placed at day 3 and at day 14. Inverse-probability weights came from a logistic model for having a valid enzyme pair, using age, sex, SOFA, a missing-SOFA indicator, alcohol use disorder, and care unit. Weights were stabilized and truncated at the 1st and 99th percentiles, then applied to a Cox model of every measured stay. Enzyme bounds of AST 5–5,000 U/L and ALT 5–2,000 U/L were compared with wider bounds (1–10,000) and tighter bounds (10–1,000), keeping the same quartile cutpoints. Discrimination was Harrell’s C-index for the covariate model with and without the log ratio, averaged across the 30 imputed datasets, at 28 days and within the first 7 days. Analyses were run in Python 3.14.5 with pandas 3.0.4, patsy 1.0.2, SciPy 1.18.0, statsmodels 0.14.6, scikit-learn 1.9.0, and lifelines 0.30.0. Figures were drawn in Matplotlib.

### Ethics statement

The collection of MIMIC-IV was approved by the Institutional Review Boards of Beth Israel Deaconess Medical Center and the Massachusetts Institute of Technology, and informed consent was waived because the records are de-identified [21,22]. The eICU Collaborative Research Database is a de-identified dataset released under the HIPAA safe-harbor provision; Privacert certified the re-identification risk (HIPAA Certification no. 1031219-2) [23,24]. This analysis used only those files and involved no direct patient contact.

## Results

### Cohort

After exclusion of the prespecified chronic liver-disease codes, MIMIC-IV contained 33,401 eligible sepsis stays, including cardiac surgical intensive care (Figure 1). A valid AST and ALT pair was available in 16,621 stays (49.8%). Of the other stays, 16,258 had neither enzyme in the first 24 hours, 192 had only AST, 35 had only ALT, and 295 had both values but fell outside AST 5–5,000 U/L or ALT 5–2,000 U/L. None of the analyzed stays had a death recorded on a calendar day before the enzyme measurement. Five stays had an enzyme chart time later than the hospital death time and were removed before analysis, leaving 16,616 patients. Table 1 gives their characteristics. Median age was 68.9 years (interquartile range 57.4–79.6), 9,400 (56.6%) were men, and the median SOFA score was 6 (4–9). The median ratio was 1.4 (1.0–2.0). From the lowest to the highest quartile, median AST was 35 U/L and then 60 U/L, and median ALT was 48 U/L and then 21 U/L. Median bilirubin stayed near 0.6–0.7 mg/dL. Crude mortality within 28 days of the enzyme measurement was 825/4,157 (19.8%), 938/4,170 (22.5%), 1,135/4,136 (27.4%), and 1,351/4,153 (32.5%).

Table 2 compares these 16,616 stays with the 16,780 eligible stays that had no valid pair. Using 28-day death from ICU admission for both groups, mortality was 4,239/16,616 (25.5%) when the enzymes were analyzed and 2,285/16,780 (13.6%) when no valid pair was available. Measured stays had a higher SOFA score (median 6 versus 4) and were less often in the cardiac surgical intensive care unit (7.8% versus 28.4%). Alcohol use disorder was recorded in 1,250 of the 16,616 analyzed stays (7.5%), viral hepatitis in 676, and fatty liver in 237.

### Unadjusted survival

Kaplan–Meier curves for all 16,616 patients, timed from the day of enzyme measurement, separated in quartile order (Figure 2A; log-rank P < 0.001). The highest quartile left the other curves during the first week.

### Adjusted association

In the imputed model, hazard ratios versus the lowest quartile were 1.04 (95% confidence interval 0.95–1.14), 1.17 (1.07–1.29), and 1.43 (1.31–1.57) (Figure 2C, Table 3). The complete-case model gave a highest-quartile hazard ratio of 1.36 (1.23–1.51). On the same complete cases, a restricted cubic spline of the log ratio rose across the observed range relative to the median (Figure 2E). The spline association was significant (P < 0.001), without a significant departure from a log-linear term (P = 0.20).

In the time-varying model the highest-quartile hazard ratio was 1.81 (1.60–2.05) during days 0–7 and 1.13 (0.99–1.28) during days 7–28. The later hazard was 0.62 times the first-week hazard (95% confidence interval 0.52–0.74; p < 0.001). When age, SOFA, lactate, and bilirubin were also allowed to change after day 7, the corresponding hazard ratios were 1.67 (1.48–1.89) and 1.23 (1.08–1.40), and their ratio was 0.74 (0.62–0.88). Using day 3 or day 14 as the cut yielded the same early pattern (Figure 4A, Table 5). With the cut at day 3, the highest-quartile hazard ratios were 2.04 (1.72–2.41) before the cut and 1.25 (1.12–1.38) after it. With the cut at day 14, they were 1.62 (1.46–1.79) and 0.99 (0.83–1.19). In absolute terms, 373/4,157 (9.0%) of the lowest quartile and 819/4,153 (19.7%) of the highest quartile died in the first week. Among patients still alive at day 7, the corresponding risks over the next 21 days were 452/3,784 (11.9%) and 532/3,334 (16.0%) (Figure 4B).

The first enzyme pair was recorded a median of 3.1 hours after ICU admission (interquartile range 1.1–8.1), and 11,085/16,616 pairs (66.7%) were recorded within 6 hours. Of the 4,249 deaths within 28 days, 3,363 (79.1%) occurred in hospital and had a clock time. On the calendar day of measurement, 269 patients died (36, 49, 74, and 110 from the lowest quartile to the highest). Removing those deaths left 16,347 patients. In that cohort the same time-varying model gave a highest-quartile hazard ratio of 1.81 (1.59–2.06) in the first week and 1.15 (1.01–1.30) from day 7 to day 28. The later hazard was 0.63 times the first-week hazard (0.53–0.76; p < 0.001) (Table 5). The first-week association and its excess over the later period therefore both remained after measurement-day deaths were set aside.

Of the 16,616 patients, 396 died within 24 hours after ICU admission. Excluding them left 16,220. The five stays removed for chart time were deaths in that same interval, so they are not part of this 16,220. Follow-up was still counted from the day of enzyme measurement. In that group, the same day-7 split, with quartiles included in the imputation, gave a highest-quartile hazard ratio of 1.82 (1.59–2.08) in the first week and 1.17 (1.03–1.32) from day 7 to day 28. Weighting all 16,616 analyzed stays for the probability of having a valid pair, in a model without lactate, gave a highest-quartile hazard ratio of 1.46 (1.32–1.61), compared with 1.55 (1.42–1.69) in the unweighted analysis of the same patients (effective sample size 12,993). With wider enzyme bounds of 1–10,000 U/L, the complete-case highest-quartile hazard ratio was 1.33 (1.20–1.47; 11,104 patients). With tighter bounds of 10–1,000 U/L, it was 1.41 (1.27–1.58; 9,716 patients).

On the 10,890 complete cases, a joint model entered log AST and log ALT separately. Per standard deviation, the hazard ratio was 1.38 (1.28–1.48) for AST and 0.81 (0.75–0.87) for ALT. The same joint model in the 30 imputed datasets on the 16,616 patients gave 1.42 (1.33–1.52) for AST and 0.77 (0.72–0.82) for ALT. Entering those imputed datasets as the log ratio and the mean of the two log enzymes gave 1.16 (1.13–1.20) per standard deviation of the ratio and 1.09 (1.06–1.13) per standard deviation of the mean level. On the complete cases, a likelihood-ratio test rejected a model that kept the ratio and dropped the mean level (P < 0.001). Adding the log ratio to the covariate model changed Harrell’s C-index from 0.714 to 0.717 at 28 days and from 0.748 to 0.752 for death within 7 days. Adding the Charlson index, weight, albumin, norepinephrine, invasive ventilation, and dialysis, again with the 30 imputed datasets, gave a highest-quartile hazard ratio of 1.33 (1.22–1.46). Creatine kinase was available in 38% of stays. In that complete-case subset, with log creatine kinase in the model, the highest-quartile hazard ratio was 1.27 (1.03–1.57). Excluding alcohol use disorder, viral hepatitis, and fatty liver in a complete-case analysis yielded a highest-quartile hazard ratio of 1.40 (1.25–1.56; 9,529 patients). Hospital death in MIMIC-IV was 670/4,157 (16.1%), 777/4,170 (18.6%), 975/4,136 (23.6%), and 1,156/4,153 (27.8%). Imputed odds ratios were 1.05 (0.93–1.18), 1.27 (1.12–1.43), and 1.49 (1.32–1.67).

### Subgroups

The association did not differ by sex (interaction hazard ratio 1.03, 95% confidence interval 0.84–1.27; p = 0.75) or by lactate (1.09, 0.85–1.41; p = 0.50) (Figure 3). Crude 28-day mortality rose across quartiles in both lactate groups (Figure 3A), and the interaction estimates are shown in Figure 3B and Table 4. It tended to be weaker at age 65 years or older (0.80, 0.64–1.00; p = 0.054) and in the cardiac surgical intensive care unit (0.61, 0.36–1.03; p = 0.067).

### Replication in eICU

eICU contained 22,405 eligible first stays after the liver-disease exclusion, and 12,374 had a valid ratio. Hospital mortality in these 12,374 patients was 424/3,881 (10.9%), 463/3,033 (15.3%), 538/2,715 (19.8%), and 862/2,745 (31.4%) (Figure 2B). Imputed odds ratios versus the lowest quartile, with hospital-clustered standard errors across 188 hospitals, were 1.15 (0.99–1.35), 1.30 (1.14–1.49), and 1.95 (1.69–2.26) (Figure 2D, Table 3). In the complete-case model (7,092 patients, 1,565 hospital deaths) the highest-quartile odds ratio was 1.96 (1.64–2.36); clustering by 157 hospitals gave 1.96 (1.63–2.36). The MIMIC-IV hospital-death odds ratio for the highest quartile was 1.49, so both databases show the same direction on one endpoint.

## Discussion

### Principal findings

In patients with sepsis who had both enzymes measured on ICU admission, a higher AST/ALT ratio was associated with higher 28-day mortality. Our study adds two findings to earlier work [15–17]. First, the association was strongest in the first week (highest-quartile hazard ratio 1.81) and weaker thereafter (1.13), a pattern that persisted when the cut-point was moved to day 3 or day 14. Second, the ratio did not contain all of the prognostic information in the two enzymes. Written as the log ratio and the log mean, the mean level retained a hazard ratio of 1.09 per standard deviation. The same model, written as log AST and log ALT, associated higher AST and lower ALT with mortality in opposite directions (hazard ratios 1.42 and 0.77 per standard deviation).

The association was consistent across sensitivity analyses, including further adjustment for comorbidity, nutritional status, and organ support, exclusion of recorded alcohol use disorder, viral hepatitis, and fatty liver, and weighting for the probability of enzyme measurement. It was also replicated for hospital mortality in eICU using cutpoints fixed in MIMIC-IV.

In absolute terms, 19.7% of patients in the highest quartile died in the first week, compared with 9.0% in the lowest; among patients alive at day 7, the gap narrowed (16.0% versus 11.9%), which may partly reflect depletion of high-risk patients after early deaths.

### What the enzyme pattern shows

The enzyme profile in the highest quartile was modest: median AST was 60 U/L and bilirubin was unchanged, far below the levels typical of hypoxic hepatitis [5–7,26]. A high ratio can be produced by a higher AST, a lower ALT, or both, so the quartile medians alone do not separate those contributions. In the joint model, each standard deviation of log AST was associated with a higher hazard and each standard deviation of log ALT with a lower hazard.

A low ALT can mark reduced physiological reserve. In older adults it tracks frailty, disability, and sarcopenia, and low ALT predicts higher long-term mortality [27,28]. The association here was not stronger at age 65 years or older (interaction hazard ratio 0.80, 0.64–1.00), so old-age frailty alone does not account for it. Albumin and weight, together with the Charlson index and markers of shock and respiratory failure, moved the highest-quartile hazard ratio from 1.43 to 1.33. ALT activity also depends on pyridoxal 5′-phosphate, and deficiency of that cofactor can lower the measured value [29]. We did not have vitamin B6 levels, so this remains a possible contributor rather than a demonstrated pathway.

### Relation to previous studies

Choi and colleagues reported an adjusted hazard ratio of 1.46 (95% CI 1.30–1.64) for an emergency-department ratio of at least 1.84 in septic shock [16]. Yang and colleagues described a nonlinear association with hospital death in older patients from MIMIC-IV [17]. Wang and Xu described the overall association and its nonlinear shape in a general adult intensive-care cohort drawn from the same two databases [15]. Chen, Luo, Yu, and Palmowski studied the ratio inside sepsis-associated kidney or liver injury [11–14]. None of those analyses tested whether the association is larger early in the first month than later, and none tested whether the ratio contains the prognostic information in both enzymes. These two questions are addressed by the present study. eICU repeats the hospital-death association under cutpoints fixed in MIMIC-IV.

### Implications

Adding the ratio to a model with severity covariates improved discrimination only marginally (C-index change 0.003 at 28 days and 0.004 at 7 days). The ratio is therefore better viewed as a routinely available marker of early mortality risk than as a stand-alone prediction tool. Because the overall enzyme level also carried prognostic information, the two values read together are more informative than the ratio alone.

### Limitations

**Selection and missing data.** Enzymes were measured in half of eligible stays. Most of the omitted stays had neither AST nor ALT ordered in the first 24 hours; a few hundred had only one enzyme or a value outside the prespecified bounds. Measured patients were more severely ill, and lactate was missing more often in the lowest quartile (40.0%) than in the highest (24.8%) (Table 6). The estimates apply to patients who had both enzymes measured within 24 hours of admission. Multiple imputation assumes that missing covariates are missing at random given the variables in the imputation model, and it does not remove bias from a decision to order the test. Weighting for the probability of measurement yielded a hazard ratio of 1.46.

**Measurement and confounding.** ICD codes miss some alcohol use disorder and some chronic liver disease. eICU identified a chronic liver-disease code in 1.9% of first sepsis stays, against 11.0% in MIMIC-IV, so some chronic liver disease likely remains, particularly in eICU, and could inflate the apparent association. Death after discharge is recorded as a date, without a clock time. Five stays with an enzyme chart time later than the hospital death time were removed before the models were fit. First-day SOFA, norepinephrine, ventilation, and dialysis can be recorded after the enzyme measurement. If these variables lie on the pathway from the enzyme pattern to death, adjusting for them may attenuate the estimate; if they act as colliders, the direction of bias is uncertain.

**Generalizability.** Both databases are from the United States, and the eICU stays date from 2014–2015. MIMIC-IV is a single center. The eICU analysis repeats the hospital-death association with a different sepsis definition. It could not assess the time-varying association, because eICU lacks post-discharge survival data. The cardiac surgical comparison was too imprecise to support a separate conclusion.

## Conclusion

Among patients with sepsis who had both enzymes measured, a higher admission AST/ALT ratio was associated with higher short-term mortality, most strongly in the first week. The ratio did not capture all of the prognostic information in AST and ALT. The association was replicated for hospital mortality in eICU, but the ratio added little discrimination beyond severity of illness.

## Acknowledgments

The authors thank the teams responsible for developing and maintaining MIMIC-IV, eICU-CRD, and PhysioNet.

## Data availability statement

MIMIC-IV v3.1 and eICU-CRD v2.0 are available through PhysioNet to credentialed users who complete the required training and data-use procedures. The code used for the cohorts, multiple imputation, time-split models, and figures in this manuscript is available at https://github.com/shupengqin/ast-alt-sepsis-mortality. Patient-level records are not included.

## Funding

The authors received no specific funding for this study.

## Conflict of interest

The authors declare that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.

## Author contributions

Xiaoye Xu, Tielong Chen, Pengqin Shu, and Zhao Liu contributed equally. Xiaoye Xu and Tielong Chen built the cohort and ran the analysis. Pengqin Shu interpreted the clinical findings and revised the manuscript. Zhao Liu prepared the figures. Fei Ying and Cheng Bao checked the analysis. Chao Chen contributed to data curation and clinical interpretation. Jing Sun supervised the study and led the writing. All authors read and approved the final manuscript.

## Figure legends

Figure 1. Cohort construction. Cardiac surgical and cardiothoracic units are retained in both databases. Quartile membership uses the MIMIC-IV cutpoints 1.04, 1.44, and 2.03 in both databases. In MIMIC-IV, 16,621 stays had a valid enzyme pair and 5 were then removed because the chart time followed the hospital death time, leaving 16,616. The eICU analysis count is 12,374.

Figure 2. (A) Kaplan–Meier survival in the 16,616 MIMIC-IV patients retained after removal of 5 post-death chart times, timed from the day of enzyme measurement. The log-rank test across quartiles gave P < 0.001. Numbers at risk count patients still under observation at the start of that day (follow-up time at least 0, 7, 14, 21, or 28 days). The day-7 survivors reported in the text are those followed beyond day 7, so they are fewer than the day-7 risk set by the deaths and censoring recorded on day 7. (B) Crude mortality by quartile. MIMIC-IV is 28-day death. eICU is hospital death. (C) Multiply imputed MIMIC-IV hazard ratios versus quartile 1. (D) Multiply imputed eICU odds ratios versus quartile 1. (E) Complete-case spline in MIMIC-IV, relative to the median ratio, among patients with complete covariates. The likelihood-ratio tests gave P < 0.001 for the spline association and P = 0.20 for nonlinearity.

Figure 3. (A) Crude 28-day mortality by admission lactate within each ratio quartile. (B) Interaction hazard ratios for the highest-versus-lowest contrast in one sample of 5,423 complete cases. A value below 1 means a smaller association in the named group. The cardiac surgical interval includes 1.

Figure 4. Alternative time cuts in MIMIC-IV, each fit to the same follow-up. (A) Imputed hazard ratios for the highest versus the lowest quartile before and after a cut at day 3, day 7, or day 14. (B) Percent dying in the first 7 days, and percent dying over the next 21 days among patients alive at day 7.

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

15. Wang Y, Xu Y. Association between aspartate aminotransferase to alanine aminotransferase ratio and 28-day mortality of ICU patients: a retrospective cohort study from MIMIC-IV database. PLoS One. 2025;20:e0324904. doi:10.1371/journal.pone.0324904.

16. Choi S, Nah S, Suh GJ, Choi SH, Chung SP, Kim WY, et al. Prognostic value of the AST/ALT ratio in patients with septic shock: a prospective, multicenter, registry-based observational study. Diagnostics (Basel). 2025;15:1773. doi:10.3390/diagnostics15141773.

17. Yang H, Meng L, Xin S, Chang C, Zhao X, Guo B, et al. Non-linear association between aspartate aminotransferase to alanine aminotransferase ratio and mortality in critically ill older patients: a retrospective cohort study. PLoS One. 2023;18:e0293749. doi:10.1371/journal.pone.0293749.

18. von Elm E, Altman DG, Egger M, Pocock SJ, Gøtzsche PC, Vandenbroucke JP, et al. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement: guidelines for reporting observational studies. Lancet. 2007;370:1453-1457. doi:10.1016/S0140-6736(07)61602-X.

19. Vandenbroucke JP, von Elm E, Altman DG, Gøtzsche PC, Mulrow CD, Pocock SJ, et al. Strengthening the Reporting of Observational Studies in Epidemiology (STROBE): explanation and elaboration. Epidemiology. 2007;18:805-835. doi:10.1097/EDE.0b013e3181577511.

20. Benchimol EI, Smeeth L, Guttmann A, Harron K, Moher D, Petersen I, et al. The REporting of studies Conducted using Observational Routinely-collected health Data (RECORD) statement. PLoS Med. 2015;12:e1001885. doi:10.1371/journal.pmed.1001885.

21. Johnson AEW, Bulgarelli L, Shen L, Gayles A, Shammout A, Horng S, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10:1. doi:10.1038/s41597-022-01899-x.

22. Johnson A, Bulgarelli L, Pollard T, Gow B, Moody B, Horng S, et al. MIMIC-IV (version 3.1) [dataset]. PhysioNet; 2024. doi:10.13026/kpb9-mt58.

23. Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU Collaborative Research Database, a freely available multi-center database for critical care research. Sci Data. 2018;5:180178. doi:10.1038/sdata.2018.178.

24. Pollard T, Johnson A, Raffa J, Celi LA, Badawi O, Mark R. eICU Collaborative Research Database (version 2.0) [dataset]. PhysioNet; 2019. doi:10.13026/C2WM1R.

25. Vincent JL, Moreno R, Takala J, Willatts S, De Mendonça A, Bruining H, et al. The SOFA (Sepsis-related Organ Failure Assessment) score to describe organ dysfunction/failure. Intensive Care Med. 1996;22:707-710. doi:10.1007/BF01709751.

26. Drolz A, Horvatits T, Michl B, Roedl K, Schellongowski P, Holzinger U, et al. Statin therapy is associated with reduced incidence of hypoxic hepatitis in critically ill patients. J Hepatol. 2014;60:1187-1193. doi:10.1016/j.jhep.2014.01.019.

27. Vespasiani-Gentilucci U, De Vincentis A, Ferrucci L, Bandinelli S, Antonelli Incalzi R, Picardi A. Low alanine aminotransferase levels in the elderly population: frailty, disability, sarcopenia, and reduced survival. J Gerontol A Biol Sci Med Sci. 2018;73:925-930. doi:10.1093/gerona/glx126.

28. Ramaty E, Maor E, Peltz-Sinvani N, Brom A, Grinfeld A, Kivity S, et al. Low ALT blood levels predict long-term all-cause mortality among adults. A historical prospective cohort study. Eur J Intern Med. 2014;25:919-921. doi:10.1016/j.ejim.2014.10.019.

29. Youssef EM, Wu GY. Subnormal serum liver enzyme levels: a review of pathophysiology and clinical significance. J Clin Transl Hepatol. 2024;12:428-435. doi:10.14218/JCTH.2023.00446.
