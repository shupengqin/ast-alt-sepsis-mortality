Admission AST/ALT ratio and 28-day mortality in sepsis: a two-database cohort study

Xiaoye Xu1#, Tielong Chen2#, Pengqin Shu2#, Zhao Liu2#, Fei Ying2, Cheng Bao2, Chao Chen2, Jing Sun2*

Xiaoye Xu, Tielong Chen, Pengqin Shu, and Zhao Liu contributed equally to this work and share first authorship.

1 Affiliated Mental Health Center & Hangzhou Seventh People's Hospital, Zhejiang University School of Medicine, Hangzhou, Zhejiang, China.

2 Hangzhou TCM Hospital Affiliated to Zhejiang Chinese Medical University, Hangzhou, Zhejiang, China.

*Corresponding author: Jing Sun, Hangzhou TCM Hospital Affiliated to Zhejiang Chinese Medical University, No. 453 Tiyuchang Road, Xihu District, Hangzhou City, Zhejiang Province, China. Email: sj524258755@163.com.

# Short title

AST/ALT ratio and 28-day mortality in sepsis

# Article type

Original Research

Main-text word count: 4,806; Abstract word count: 288; Figures: 4; Tables: 7.

# Abstract

Background: The aspartate aminotransferase (AST) to alanine aminotransferase (ALT) ratio is associated with mortality in sepsis. However, whether this association is stronger early, and whether the ratio fully captures the prognostic information contained in AST and ALT, remains unclear. Methods: Adults with sepsis in the Medical Information Mart for Intensive Care IV (MIMIC-IV) and the eICU Collaborative Research Database (eICU) were included after exclusion of a prespecified set of chronic liver-disease codes. Patients were grouped into quartiles of the first AST/ALT ratio within 24 hours (cutpoints 1.04, 1.44, and 2.03). Follow-up for 28-day death started on the day of enzyme measurement. Models included age, sex, the Sequential Organ Failure Assessment (SOFA) score, lactate, bilirubin, creatinine, care unit, and alcohol use disorder, with multiple imputation for missing covariates. The eICU outcome was hospital death. Results: Among 16,616 MIMIC-IV patients, 28-day mortality rose from 19.8% to 32.5%. The imputed highest-quartile hazard ratio was 1.43 (95% CI 1.31–1.57) over 28 days, 1.81 (1.60–2.05) in the first week, and 1.13 (0.99–1.28) from day 7 to day 28. In a 24-hour landmark analysis, the corresponding hazard ratios were 1.82 (1.58–2.09) and 1.20 (1.06–1.36). Hazard ratios per standard deviation were 1.16 (1.13–1.20) for the log ratio and 1.09 (1.06–1.13) for the mean log enzyme level. In 12,374 eICU patients, hospital mortality rose across quartiles from 10.9% to 31.4% (adjusted odds ratio for the highest quartile 1.95, 1.69–2.26). Conclusions: Among patients with sepsis who had both enzymes measured, a higher admission AST/ALT ratio was associated with higher short-term mortality; the association was strongest early and attenuated over follow-up. The ratio did not capture all of the prognostic information in the two enzymes and added little discrimination. The association with hospital mortality was replicated in eICU.

Keywords: sepsis; AST/ALT ratio; De Ritis ratio; mortality; intensive care; MIMIC-IV; eICU

# Introduction

Sepsis remains a leading cause of death in intensive care. Sepsis-3 defines it as life-threatening organ dysfunction caused by a dysregulated host response to infection [1]. The 2021 Surviving Sepsis Campaign guidelines concentrate early care on infection control and restoration of perfusion [2]. Liver dysfunction, however, is not the target of any specific treatment recommendation. The liver can amplify inflammation, alter drug handling, and fail as a target organ [3,4]. Hypoxic (ischemic) hepatitis represents its most severe circulatory form: aminotransferases rise sharply in cardiac, circulatory, or respiratory failure once other causes of hepatic necrosis are set aside [5–7]. In a medical intensive care series, hypoxic hepatitis was frequently accompanied by acute kidney injury and high short-term mortality [8].

Aspartate aminotransferase (AST) and alanine aminotransferase (ALT) differ in tissue distribution. AST is present in hepatocytes, skeletal muscle, and myocardium, whereas ALT is more concentrated in the liver [5]. Their ratio, often called the De Ritis ratio, has therefore been used as a severity marker beyond hepatology, including in COVID-19 and in traumatic brain injury [9,10]. In sepsis, most studies have examined the ratio within narrower syndromes, such as sepsis-associated acute kidney injury [11,12] and sepsis-associated liver injury [13,14].

The association between the ratio and mortality has also been examined beyond a single syndrome. In patients with sepsis or septic shock, Schupp and colleagues measured the ratio repeatedly during the first week and found that a ratio above 1.8 was associated with 30-day death [15]. In a general adult intensive care population from MIMIC-IV, not restricted to sepsis and with eICU as a sensitivity cohort, Wang and Xu reported a nonlinear association with 28-day death [16]. In a prospective multicenter registry of septic shock, Choi and colleagues found that an emergency-department ratio of at least 1.84 was associated with 28-day death and added information to the SOFA score [17]. Among older patients in MIMIC-IV, Yang and colleagues described a nonlinear association with hospital death that flattened above a ratio of about 1.8 [18]. These studies characterized the overall association, but none examined whether the association of a single admission ratio changes over follow-up, or whether the ratio captures the prognostic information carried by both enzymes. A ratio above 2 is also the classic pattern of alcoholic liver injury, so alcohol use must be accounted for before the pattern is attributed to sepsis.

We therefore aimed to examine whether the association between the admission AST/ALT ratio and 28-day mortality is stronger early than later, and whether the ratio fully captures the prognostic information contained in AST and ALT. Quartile cutpoints fixed in MIMIC-IV were applied to hospital mortality in eICU.

# Materials and Methods

## Study design and reporting

This was a retrospective cohort study in two databases. Reporting followed the STROBE statement and the RECORD statement for routinely collected health data [19–21]. The time-varying analysis was specified in MIMIC-IV. eICU was used to replicate the association with hospital death. No sample-size calculation was performed.

## Data sources

The Medical Information Mart for Intensive Care IV (MIMIC-IV) version 3.1 contains electronic health-record data from intensive care stays at one United States academic medical center [22,23]. The eICU Collaborative Research Database (eICU) version 2.0 contains stays from 208 United States hospitals in 2014–2015 [24,25]. Both databases were accessed under the PhysioNet credentialed-data agreement.

## Participants

MIMIC-IV eligibility required a first intensive care stay, age of at least 18 years, and Sepsis-3 status in the derived sepsis table [1]. The liver exclusion used ICD-9 codes 571.2, 571.5, 571.6, 571.8, and 572.2, and ICD-10 codes beginning with K70 or K74, or K71.7. Viral hepatitis and fatty liver were not in this list. Cardiac surgical intensive care stays were retained. In eICU, the first unit stay of a hospitalization was eligible when the diagnosis list contained “sepsis” or “septic.” We excluded records coded as cirrhosis or chronic liver disease, and ICD-9 571.2, 571.5, and 571.6. Cardiac, cardiothoracic, and cardiac surgical units were retained.

## Exposure

The exposure was the AST/ALT ratio from the first pair of results within 24 hours after ICU admission in MIMIC-IV or unit admission in eICU. This value is referred to as the admission ratio. In MIMIC-IV the two enzymes come from the same laboratory record. In eICU the first AST and the first ALT had the same time offset in 14,773 of 14,803 stays with both results (99.8%). Pairs were retained when AST was between 5 and 5,000 U/L and ALT was between 5 and 2,000 U/L. MIMIC-IV quartiles were formed at the 25th, 50th, and 75th percentiles and then fixed at 1.0357, 1.4444, and 2.0294. The same cutpoints were applied to the Kaplan–Meier curves and in eICU. Follow-up started on the calendar day of that enzyme pair. Patients whose death was recorded on an earlier calendar day were not included. The time from ICU admission to the enzyme chart time was summarized. For deaths in hospital, that chart time was compared with the recorded death time. Five stays in which the chart time followed the death time were removed before the models were fit. Deaths after discharge were recorded as dates without clock times.

## Covariates

The adjusted MIMIC-IV model included age, sex, the first-day Sequential Organ Failure Assessment (SOFA) score [26], the first lactate within 24 hours, the logarithm of the first bilirubin, the first creatinine, care unit, and alcohol use disorder. Alcohol use disorder was defined by ICD-10 F10 and ICD-9 303 or 305.0. Care unit was coded as medical, coronary, cardiac surgical, surgical or trauma, or other. A second model entered standardized log AST and standardized log ALT together. The same information was also entered as the standardized log ratio and the standardized mean of log AST and log ALT. These two forms are reparameterizations of the same model. The coefficient on the mean level is zero only when log AST and log ALT have coefficients of equal size and opposite sign, which is the case in which the ratio alone carries the association. A further model added the Charlson comorbidity index, first-day weight, albumin, norepinephrine in the first 24 hours, invasive ventilation in the first 24 hours, and dialysis in the first 24 hours. Creatine kinase was available for a minority of stays and was examined only in that subset. eICU models used the Acute Physiology and Chronic Health Evaluation (APACHE) score in place of SOFA and included age, sex, lactate, bilirubin, creatinine, care-unit group, and an alcohol-related diagnosis. That diagnosis was a string containing “alcohol” or ICD-9 291, 303, or 305.0.

## Outcomes

In MIMIC-IV, the outcome was death within 28 days of the enzyme measurement, using the recorded date of death. That date includes deaths matched from state records after hospital discharge. Patients alive at 28 days, or dying later, were censored at 28 days. Deaths on the day of measurement were assigned a time of half a day. eICU does not provide the same post-discharge vital status, so the replication outcome was death before hospital discharge. The same hospital-death endpoint was also fit in MIMIC-IV, so the two databases could be compared on one definition.

## Statistical analysis

Kaplan–Meier curves and crude rates included every patient with a valid ratio and used the same cutpoints in both databases. Lactate was missing in 31.9% of the 16,616 MIMIC-IV patients and was the main reason a complete-case model kept 10,890 stays. Adjusted hazard ratios were estimated from 30 imputed datasets, in line with the largest missing fraction of about 32%. Imputation used chained equations with posterior draws. The imputation model contained the quartile indicators, the analysis covariates, the 28-day death indicator, and the Nelson–Aalen cumulative hazard. The quartiles were complete, so they entered the imputation as predictors and were not themselves imputed. Alcohol use disorder came from diagnosis codes and had no missing values, so it was not imputed. The same imputed datasets were used for the 28-day model and for one time-varying Cox model. In that model, the quartile coefficients were allowed to change after day 7, whereas the coefficients of the other covariates were held constant. Complete-case estimates are shown beside the imputed ones. The creatine-kinase subset and the analysis that removed alcohol use disorder, viral hepatitis, and fatty liver are complete-case analyses. In the creatine-kinase subset, log creatine kinase was a covariate.

Hazard ratios for the second, third, and highest quartiles used the lowest quartile as the reference. A rank Schoenfeld test on the complete-case Cox model rejected proportional hazards for the highest quartile (P = 0.006), and also for SOFA, lactate, age, and bilirubin. For this reason, the time-varying model allowed the quartile coefficients to change after day 7. The day-7 cut-point was chosen after inspection of the survival curves and was not prespecified; cut-points at day 3 and day 14 were examined as robustness checks. To describe the change over time without a fixed cut-point, two complete-case models were fitted: one in which each quartile log hazard ratio changed linearly with log time, compared with a model of constant coefficients by a likelihood-ratio test, and one in which the highest-quartile log hazard ratio changed linearly with time, with a change in slope at day 7. Because these models are not nested, they were compared by the Akaike information criterion (AIC). Confidence bands were derived by the delta method. The coefficients of age, SOFA, lactate, and bilirubin were held constant in both models. A sensitivity analysis let age, SOFA, lactate, and bilirubin change after day 7 as well. Restricted cubic splines of the log ratio described the shape of the association. The plotted curve shows the hazard ratio relative to the median ratio.

Subgroup models contrasted the highest and lowest quartiles and included a product term so that an interaction P value could be reported for sex, age of at least 65 years, lactate of at least 4 mmol/L, and cardiac surgical intensive care. These comparisons were not registered in advance. A sensitivity analysis excluded alcohol use disorder, viral hepatitis (ICD-10 B15–B19, ICD-9 070), and fatty liver (ICD-10 K76.0 or K75.8).

eICU hospital death was analyzed with the same chained-equation imputation, including the quartile indicators, on all 12,374 patients with a ratio. The reported standard errors are robust and clustered by hospital. The complete-case logistic model is reported beside the imputed one. A sensitivity analysis excluded patients who had already died within 24 hours after ICU admission and refit the day-7 split in that group. Follow-up was still counted from the day of enzyme measurement. In a landmark analysis of the same group, follow-up restarted at 24 hours after ICU admission and ended at day 28 after ICU admission; in this analysis, day 7 and day 28 were counted from ICU admission, so the early interval spanned about 6 days. A further sensitivity analysis removed patients who died on the calendar day of measurement and refit the same day-7 split, including the later-to-early ratio. The time-varying model was repeated with the split placed at day 3 and at day 14. Inverse-probability weights came from a logistic model for having a valid enzyme pair, using age, sex, SOFA, a missing-SOFA indicator, alcohol use disorder, and care unit. Weights were stabilized, truncated at the 1st and 99th percentiles, and applied to a Cox model of the analyzed stays. Standardized mean differences compared the 16,616 analyzed stays, before and after weighting, with all 33,396 eligible stays. The denominator was the unweighted pooled standard deviation of the analyzed stays and the stays with no valid pair. Enzyme bounds of AST 5–5,000 U/L and ALT 5–2,000 U/L were compared with wider bounds (1–10,000) and tighter bounds (10–1,000), keeping the same quartile cutpoints. Discrimination was Harrell’s C-index for the covariate model, for that model plus the log ratio, and for that model plus log AST and log ALT together, averaged across the 30 imputed datasets, at 28 days and within the first 7 days. Percentile intervals for the 28-day increments were obtained from 120 bootstrap samples of the 10,890 complete cases, refitting all three models in each sample; these intervals describe the apparent increment and are not optimism-corrected. Analyses were run in Python 3.14.5 with pandas 3.0.4, patsy 1.0.2, SciPy 1.18.0, statsmodels 0.14.6, scikit-learn 1.9.0, and lifelines 0.30.0. Figures were drawn in Matplotlib.

# Results

## Cohort

After exclusion of the prespecified chronic liver-disease codes, MIMIC-IV contained 33,401 eligible sepsis stays, including cardiac surgical intensive care (Figure 1). A valid AST and ALT pair was available in 16,621 stays (49.8%). Of the other stays, 16,258 had neither enzyme in the first 24 hours, 192 had only AST, 35 had only ALT, and 295 had both values but fell outside AST 5–5,000 U/L or ALT 5–2,000 U/L. None of the analyzed stays had a death recorded on a calendar day before the enzyme measurement. Five stays had an enzyme chart time later than the hospital death time and were removed before analysis, leaving 16,616 patients. Table 1 gives their characteristics. Median age was 68.9 years (interquartile range 57.4–79.6), 9,400 (56.6%) were men, and the median SOFA score was 6 (4–9). The median ratio was 1.4 (1.0–2.0). From the lowest to the highest quartile, median AST increased from 35 to 60 U/L and median ALT decreased from 48 to 21 U/L, whereas median bilirubin remained at 0.6–0.7 mg/dL. Crude mortality within 28 days of the enzyme measurement was 825/4,157 (19.8%), 938/4,170 (22.5%), 1,135/4,136 (27.4%), and 1,351/4,153 (32.5%).

Table 2 compares these 16,616 stays with the 16,780 eligible stays that had no valid pair. Using 28-day death from ICU admission for both groups, mortality was 4,239/16,616 (25.5%) when the enzymes were analyzed and 2,285/16,780 (13.6%) when no valid pair was available. Measured stays had a higher SOFA score (median 6 versus 4) and were less often in the cardiac surgical intensive care unit (7.8% versus 28.4%). Alcohol use disorder was recorded in 1,250 of the 16,616 analyzed stays (7.5%), viral hepatitis in 676, and fatty liver in 237.

## Unadjusted survival

Kaplan–Meier curves for all 16,616 patients, timed from the day of enzyme measurement, separated in quartile order (Figure 2A; log-rank P < 0.001). The curve for the highest quartile separated from the others during the first week.

## Adjusted association

In the imputed model, hazard ratios versus the lowest quartile were 1.04 (95% confidence interval 0.95–1.14), 1.17 (1.07–1.29), and 1.43 (1.31–1.57) (Figure 2C, Table 3). The complete-case model gave a highest-quartile hazard ratio of 1.36 (1.23–1.51). On the same complete cases, the restricted cubic spline of the log ratio increased across the observed range relative to the median (Figure 2E). The overall association was significant (P < 0.001), with no evidence of departure from log-linearity (P = 0.20).

In the time-varying model the highest-quartile hazard ratio was 1.81 (1.60–2.05) during days 0–7 and 1.13 (0.99–1.28) during days 7–28. The later hazard was 0.62 times the first-week hazard (95% confidence interval 0.52–0.74; P < 0.001). When age, SOFA, lactate, and bilirubin were also allowed to change after day 7, the corresponding hazard ratios were 1.67 (1.48–1.89) and 1.23 (1.08–1.40), and their ratio was 0.74 (0.62–0.88). Moving the cut-point to day 3 or day 14 yielded the same early pattern (Table 5). With the cut at day 3, the highest-quartile hazard ratios were 2.04 (1.72–2.41) before the cut and 1.25 (1.12–1.38) after it. With the cut at day 14, they were 1.62 (1.46–1.79) and 0.99 (0.83–1.19). On the 10,890 complete cases, allowing the quartile coefficients to change with log time improved model fit (likelihood-ratio P < 0.001). In this model, the highest-quartile hazard ratio was 1.95 (1.61–2.35) at day 1, 1.28 (1.15–1.43) at day 7, 1.11 (0.96–1.27) at day 14, and 0.95 (0.79–1.15) at day 28, indicating a progressively attenuating association (Figure 4A, Table 5). In the piecewise-linear model, the corresponding hazard ratios were 1.65 (1.43–1.91) at day 1, 1.19 (1.01–1.40) at day 7, and 1.28 (0.97–1.69) at day 28; the log hazard ratio changed by −0.055 per day before day 7 (95% CI −0.091 to −0.020) and by 0.004 per day thereafter (−0.013 to 0.020). Both models therefore showed a stronger early association, whereas the course after day 7 depended on the model. The log-time model had a lower AIC (56,655 versus 56,666), although it also allowed the second and third quartiles to vary with time. In absolute terms, 373/4,157 (9.0%) of the lowest quartile and 819/4,153 (19.7%) of the highest quartile died in the first week. Among patients still alive at day 7, the corresponding risks over the next 21 days were 452/3,784 (11.9%) and 532/3,334 (16.0%) (Figure 4B).

The first enzyme pair was recorded a median of 3.1 hours after ICU admission (interquartile range 1.1–8.1), and 11,085/16,616 pairs (66.7%) were recorded within 6 hours. Of the 4,249 deaths within 28 days, 3,363 (79.1%) occurred in hospital and had a clock time. On the calendar day of measurement, 269 patients died (36, 49, 74, and 110 from the lowest quartile to the highest). Removing those deaths left 16,347 patients. In that cohort the same time-varying model gave a highest-quartile hazard ratio of 1.81 (1.59–2.06) in the first week and 1.15 (1.01–1.30) from day 7 to day 28. The later hazard was 0.63 times the first-week hazard (0.53–0.76; P < 0.001) (Table 5). Thus, both the first-week association and its excess over the later period persisted after measurement-day deaths were excluded.

Of the 16,616 patients, 396 died within 24 hours after ICU admission. Excluding them left 16,220. The five stays removed because of the chart-time inconsistency had also died within this interval and are therefore not among these 16,220. Follow-up was still counted from the day of enzyme measurement. In that group, the same day-7 split, with quartiles included in the imputation, gave a highest-quartile hazard ratio of 1.82 (1.59–2.08) in the first week and 1.17 (1.03–1.32) from day 7 to day 28. In the 24-hour landmark analysis of these 16,220 patients, the highest-quartile hazard ratio was 1.82 (1.58–2.09) from 24 hours to day 7 and 1.20 (1.06–1.36) from day 7 to day 28 after ICU admission, and the later hazard was 0.66 times the earlier hazard (0.55–0.80) (Table 5). In this group, 1,772 patients died by day 7 and 2,061 between day 7 and day 28; 20 deaths that occurred within 28 days of enzyme measurement but after day 28 from ICU admission were censored. Weighting all 16,616 analyzed stays for the probability of having a valid pair, in a model without lactate, gave a highest-quartile hazard ratio of 1.46 (1.32–1.61), compared with 1.55 (1.42–1.69) in the unweighted analysis of the same patients (effective sample size 12,993). After truncation, the weights ranged from 0.55 to 3.97 (median 0.88). Compared with all 33,396 eligible stays, the largest standardized mean differences were 0.26 for SOFA and −0.28 for cardiac surgical care before weighting, and −0.01 and 0.01 after weighting (Table 7). Truncating at the 5th and 95th percentiles, instead of the 1st and 99th, gave a highest-quartile hazard ratio of 1.49 (1.36–1.64); leaving the weights untruncated gave 1.46 (1.32–1.61). Predicted measurement probabilities in the cardiac surgical units overlapped between measured and unmeasured stays. With wider enzyme bounds of 1–10,000 U/L, the complete-case highest-quartile hazard ratio was 1.33 (1.20–1.47; 11,100 patients). With tighter bounds of 10–1,000 U/L, it was 1.41 (1.27–1.58; 9,713 patients).

On the 10,890 complete cases, a joint model entered log AST and log ALT separately. Per standard deviation, the hazard ratio was 1.38 (1.28–1.48) for AST and 0.81 (0.75–0.87) for ALT. The same joint model in the 30 imputed datasets on the 16,616 patients gave 1.42 (1.33–1.52) for AST and 0.77 (0.72–0.82) for ALT. Reparameterizing this model as the log ratio and the mean of the two log enzymes gave 1.16 (1.13–1.20) per standard deviation of the ratio and 1.09 (1.06–1.13) per standard deviation of the mean level. On the complete cases, a likelihood-ratio test rejected the model that retained the ratio but omitted the mean level (P < 0.001). Adding the log ratio to the covariate model changed Harrell’s C-index from 0.714 to 0.716 at 28 days and from 0.748 to 0.752 for death within 7 days. Adding log AST and log ALT together, instead of the ratio, changed these values from 0.714 to 0.718 and from 0.748 to 0.754. In bootstrap samples of the complete cases, the 95% intervals for the apparent 28-day increment were 0.0004–0.0033 for the ratio and 0.0024–0.0063 for both enzymes relative to the covariate model, and 0.0013–0.0041 for both enzymes relative to the ratio. Adding the Charlson index, weight, albumin, norepinephrine, invasive ventilation, and dialysis, again with the 30 imputed datasets, gave a highest-quartile hazard ratio of 1.33 (1.22–1.46). Creatine kinase was available in 38% of stays. In that complete-case subset, with log creatine kinase in the model, the highest-quartile hazard ratio was 1.27 (1.03–1.57). Excluding alcohol use disorder, viral hepatitis, and fatty liver in a complete-case analysis yielded a highest-quartile hazard ratio of 1.40 (1.25–1.56; 9,526 patients). Hospital death in MIMIC-IV was 670/4,157 (16.1%), 777/4,170 (18.6%), 975/4,136 (23.6%), and 1,156/4,153 (27.8%). Imputed odds ratios were 1.05 (0.93–1.18), 1.27 (1.12–1.43), and 1.49 (1.32–1.67).

## Subgroups

The association did not differ by sex (interaction hazard ratio 1.03, 95% confidence interval 0.84–1.27; P = 0.75) or by lactate (1.09, 0.85–1.41; P = 0.50), and it tended to be weaker at age 65 years or older (0.80, 0.64–1.00; P = 0.054) and in the cardiac surgical intensive care unit (0.61, 0.36–1.03; P = 0.067) (Figure 3B, Table 4). Crude 28-day mortality rose across quartiles in both lactate groups (Figure 3A).

## Replication in eICU

eICU contained 22,405 eligible first stays after the liver-disease exclusion, and 12,374 had a valid ratio. Hospital mortality in these 12,374 patients was 424/3,881 (10.9%), 463/3,033 (15.3%), 538/2,715 (19.8%), and 862/2,745 (31.4%) (Figure 2B). Imputed odds ratios versus the lowest quartile, with hospital-clustered standard errors across 188 hospitals, were 1.15 (0.99–1.35), 1.30 (1.14–1.49), and 1.95 (1.69–2.26) (Figure 2D, Table 3). In the complete-case model (7,092 patients, 1,565 hospital deaths) the highest-quartile odds ratio was 1.96 (1.64–2.36); clustering by 157 hospitals gave 1.96 (1.63–2.36). The MIMIC-IV hospital-death odds ratio for the highest quartile was 1.49, so the association was in the same direction in both databases for this common endpoint.

# Discussion

## Principal findings

In patients with sepsis who had both enzymes measured on ICU admission, a higher AST/ALT ratio was associated with higher 28-day mortality. Our study adds two findings to earlier work [15–18]. First, the association was strongest in the first week (highest-quartile hazard ratio 1.81) and weaker thereafter (1.13). This pattern persisted when the cut-point was moved to day 3 or day 14 and in a 24-hour landmark analysis, and models with continuous time suggested a gradual attenuation rather than a threshold at day 7. Second, the ratio did not capture all of the prognostic information in the two enzymes. When the model was expressed as the log ratio and the mean log enzyme level, the mean level remained associated with mortality (hazard ratio 1.09 per standard deviation). Equivalently, higher AST and lower ALT were associated with mortality in opposite directions (hazard ratios 1.42 and 0.77 per standard deviation).

The 28-day association persisted after further adjustment for comorbidity, nutritional status, and organ support, after exclusion of recorded alcohol use disorder, viral hepatitis, and fatty liver, and after weighting for the probability of enzyme measurement (hazard ratio 1.46). It was also replicated for hospital mortality in eICU using cutpoints fixed in MIMIC-IV.

In absolute terms, 19.7% of patients in the highest quartile died in the first week, compared with 9.0% in the lowest; among patients alive at day 7, the gap narrowed (16.0% versus 11.9%), which may partly reflect depletion of high-risk patients after early deaths.

## What the enzyme pattern shows

The enzyme profile in the highest quartile was modest: median AST was 60 U/L and bilirubin was unchanged, far below the levels typical of hypoxic hepatitis [5–7,27]. A high ratio can be produced by a higher AST, a lower ALT, or both, so the quartile medians alone do not separate those contributions. In the joint model, each standard deviation of log AST was associated with a higher hazard and each standard deviation of log ALT with a lower hazard.

A low ALT may indicate reduced physiological reserve. In older adults it is associated with frailty, disability, and sarcopenia, and low ALT predicts higher long-term mortality [28,29]. The association here was not stronger at age 65 years or older (interaction hazard ratio 0.80, 0.64–1.00), so old-age frailty alone does not appear to account for it. Albumin and weight, together with the Charlson index and markers of shock and respiratory failure, attenuated the highest-quartile hazard ratio from 1.43 to 1.33. ALT activity also depends on pyridoxal 5′-phosphate, and deficiency of that cofactor can lower the measured value [30]. We did not have vitamin B6 levels, so this remains a possible contributor rather than a demonstrated pathway.

## Relation to previous studies

Schupp and colleagues reported an adjusted hazard ratio of 1.70 (95% CI 1.18–2.45) for 30-day death with a ratio above 1.8, with measurements repeated through day 7 [15]. Choi and colleagues reported an adjusted hazard ratio of 1.46 (95% CI 1.30–1.64) for an emergency-department ratio of at least 1.84 in septic shock [17]. Yang and colleagues described a nonlinear association with hospital death in older patients from MIMIC-IV [18]. Wang and Xu described the overall association and its nonlinear shape in a general adult intensive-care cohort drawn from the same two databases [16]. Chen, Luo, Yu, and Palmowski studied the ratio within sepsis-associated kidney or liver injury [11–14]. None of these studies examined whether the association of a single admission ratio changes over the first month, or whether the ratio captures the prognostic information in both enzymes; both questions are addressed here.

## Implications

Adding the ratio to a model with severity covariates increased the 28-day C-index only marginally, from 0.714 to 0.716, and adding both enzymes instead raised it to 0.718. The ratio is therefore better viewed as a routinely available marker of early mortality risk than as a stand-alone prediction tool. Although the ratio is not a sufficient statistical summary of the two enzymes, using both values did not meaningfully improve discrimination.

## Limitations

Selection and missing data. Enzymes were measured in half of eligible stays. Most of the omitted stays had neither AST nor ALT ordered in the first 24 hours; a few hundred had only one enzyme or a value outside the prespecified bounds. Measured patients were more severely ill, and lactate was missing more often in the lowest quartile (40.0%) than in the highest (24.8%) (Table 6). The estimates apply to patients who had both enzymes measured within 24 hours of admission. Weighting for the probability of measurement yielded a similar hazard ratio (1.46). No weighting can recover the association among patients who were never tested. Multiple imputation assumes that missing covariates are missing at random given the variables in the imputation model, and it does not remove bias from a decision to order the test.

Measurement and confounding. ICD codes miss some alcohol use disorder and some chronic liver disease. eICU identified a chronic liver-disease code in 1.9% of first sepsis stays, against 11.0% in MIMIC-IV, so some chronic liver disease likely remains, particularly in eICU, and could inflate the apparent association. Death after discharge is recorded as a date, without a clock time. Five stays with an enzyme chart time later than the hospital death time were removed before the models were fit. In the primary model, first-day SOFA, norepinephrine, ventilation, and dialysis may have been recorded after the enzyme measurement; the 24-hour landmark analysis, in which these variables preceded the start of follow-up, showed a similar early-versus-late pattern, although it was restricted to patients alive at 24 hours. If these variables lie on the pathway from the enzyme pattern to death, adjusting for them may attenuate the estimate; if they act as colliders, the direction of bias is uncertain.

Generalizability. Both databases are from the United States, and the eICU stays date from 2014–2015. MIMIC-IV comes from a single center. The eICU analysis replicated the hospital-death association under a different sepsis definition. It could not assess the time-varying association, because eICU lacks post-discharge survival data. The cardiac surgical comparison was too imprecise to support a separate conclusion. Finally, the day-7 cut-point was not prespecified, so the early-versus-late contrast should be regarded as exploratory.

# Conclusion

Among patients with sepsis who had both enzymes measured within 24 hours of ICU admission, a higher AST/ALT ratio was associated with higher short-term mortality. This association was strongest early and attenuated over follow-up. The ratio did not fully capture the prognostic information in AST and ALT, and its gain in discrimination was small. The association with hospital mortality was replicated in eICU, where the time-varying pattern could not be assessed.

# Acknowledgments

The authors thank the teams responsible for developing and maintaining MIMIC-IV, eICU-CRD, and PhysioNet.

# Data availability statement

MIMIC-IV v3.1 and eICU-CRD v2.0 are available through PhysioNet to credentialed users who complete the required training and data-use procedures. The code used for the cohorts, multiple imputation, time-split models, and figures in this manuscript is available at https://github.com/shupengqin/ast-alt-sepsis-mortality. Patient-level records are not included.

# Funding

The authors received no specific funding for this study.

# Conflict of interest

The authors declare that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.

# Author contributions

Xiaoye Xu, Tielong Chen, Pengqin Shu, and Zhao Liu contributed equally. Xiaoye Xu and Tielong Chen built the cohort and ran the analysis. Pengqin Shu interpreted the clinical findings and revised the manuscript. Zhao Liu prepared the figures. Fei Ying and Cheng Bao checked the analysis. Chao Chen contributed to data curation and clinical interpretation. Jing Sun supervised the study and led the writing. All authors read and approved the final manuscript.

# Ethics statement

The collection of MIMIC-IV was approved by the Institutional Review Boards of Beth Israel Deaconess Medical Center and the Massachusetts Institute of Technology, and informed consent was waived because the records are de-identified [22,23]. The eICU Collaborative Research Database is a de-identified dataset released under the HIPAA safe-harbor provision; Privacert certified the re-identification risk (HIPAA Certification no. 1031219-2) [24,25]. This analysis used only those databases and involved no direct patient contact.

# References

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

15. Schupp T, Weidner K, Rusnak J, Jawhar S, Forner J, Dulatahu F, et al. Diagnostic and prognostic value of the AST/ALT ratio in patients with sepsis and septic shock. Scand J Gastroenterol. 2023;58:392-402. doi:10.1080/00365521.2022.2131331.

16. Wang Y, Xu Y. Association between aspartate aminotransferase to alanine aminotransferase ratio and 28-day mortality of ICU patients: a retrospective cohort study from MIMIC-IV database. PLoS One. 2025;20:e0324904. doi:10.1371/journal.pone.0324904.

17. Choi S, Nah S, Suh GJ, Choi SH, Chung SP, Kim WY, et al. Prognostic value of the AST/ALT ratio in patients with septic shock: a prospective, multicenter, registry-based observational study. Diagnostics (Basel). 2025;15:1773. doi:10.3390/diagnostics15141773.

18. Yang H, Meng L, Xin S, Chang C, Zhao X, Guo B, et al. Non-linear association between aspartate aminotransferase to alanine aminotransferase ratio and mortality in critically ill older patients: a retrospective cohort study. PLoS One. 2023;18:e0293749. doi:10.1371/journal.pone.0293749.

19. von Elm E, Altman DG, Egger M, Pocock SJ, Gøtzsche PC, Vandenbroucke JP, et al. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement: guidelines for reporting observational studies. Lancet. 2007;370:1453-1457. doi:10.1016/S0140-6736(07)61602-X.

20. Vandenbroucke JP, von Elm E, Altman DG, Gøtzsche PC, Mulrow CD, Pocock SJ, et al. Strengthening the Reporting of Observational Studies in Epidemiology (STROBE): explanation and elaboration. Epidemiology. 2007;18:805-835. doi:10.1097/EDE.0b013e3181577511.

21. Benchimol EI, Smeeth L, Guttmann A, Harron K, Moher D, Petersen I, et al. The REporting of studies Conducted using Observational Routinely-collected health Data (RECORD) statement. PLoS Med. 2015;12:e1001885. doi:10.1371/journal.pmed.1001885.

22. Johnson AEW, Bulgarelli L, Shen L, Gayles A, Shammout A, Horng S, et al. MIMIC-IV, a freely accessible electronic health record dataset. Sci Data. 2023;10:1. doi:10.1038/s41597-022-01899-x.

23. Johnson A, Bulgarelli L, Pollard T, Gow B, Moody B, Horng S, et al. MIMIC-IV (version 3.1) [dataset]. PhysioNet; 2024. doi:10.13026/kpb9-mt58.

24. Pollard TJ, Johnson AEW, Raffa JD, Celi LA, Mark RG, Badawi O. The eICU Collaborative Research Database, a freely available multi-center database for critical care research. Sci Data. 2018;5:180178. doi:10.1038/sdata.2018.178.

25. Pollard T, Johnson A, Raffa J, Celi LA, Badawi O, Mark R. eICU Collaborative Research Database (version 2.0) [dataset]. PhysioNet; 2019. doi:10.13026/C2WM1R.

26. Vincent JL, Moreno R, Takala J, Willatts S, De Mendonça A, Bruining H, et al. The SOFA (Sepsis-related Organ Failure Assessment) score to describe organ dysfunction/failure. Intensive Care Med. 1996;22:707-710. doi:10.1007/BF01709751.

27. Drolz A, Horvatits T, Michl B, Roedl K, Schellongowski P, Holzinger U, et al. Statin therapy is associated with reduced incidence of hypoxic hepatitis in critically ill patients. J Hepatol. 2014;60:1187-1193. doi:10.1016/j.jhep.2014.01.019.

28. Vespasiani-Gentilucci U, De Vincentis A, Ferrucci L, Bandinelli S, Antonelli Incalzi R, Picardi A. Low alanine aminotransferase levels in the elderly population: frailty, disability, sarcopenia, and reduced survival. J Gerontol A Biol Sci Med Sci. 2018;73:925-930. doi:10.1093/gerona/glx126.

29. Ramaty E, Maor E, Peltz-Sinvani N, Brom A, Grinfeld A, Kivity S, et al. Low ALT blood levels predict long-term all-cause mortality among adults. A historical prospective cohort study. Eur J Intern Med. 2014;25:919-921. doi:10.1016/j.ejim.2014.10.019.

30. Youssef EM, Wu GY. Subnormal serum liver enzyme levels: a review of pathophysiology and clinical significance. J Clin Transl Hepatol. 2024;12:428-435. doi:10.14218/JCTH.2023.00446.

# Figure legends

Figure 1. Cohort construction. Cardiac surgical and cardiothoracic units are retained in both databases. Quartile membership uses the MIMIC-IV cutpoints 1.04, 1.44, and 2.03 in both databases. In MIMIC-IV, 16,621 stays had a valid enzyme pair and 5 were then removed because the chart time followed the hospital death time, leaving 16,616. The eICU analysis count is 12,374.

Figure 2. (A) Kaplan–Meier survival in the 16,616 MIMIC-IV patients retained after removal of 5 post-death chart times, timed from the day of enzyme measurement. The log-rank test across quartiles gave P < 0.001. Numbers at risk count patients still under observation at the start of that day (follow-up time at least 0, 7, 14, 21, or 28 days). The day-7 survivors reported in the text are those followed beyond day 7, so they are fewer than the day-7 risk set by the deaths and censoring recorded on day 7. (B) Crude mortality by quartile. MIMIC-IV is 28-day death. eICU is hospital death. (C) Multiply imputed MIMIC-IV hazard ratios versus quartile 1. (D) Multiply imputed eICU odds ratios versus quartile 1. (E) Complete-case spline in MIMIC-IV, relative to the median ratio, among patients with complete covariates. The likelihood-ratio tests gave P < 0.001 for the spline association and P = 0.20 for nonlinearity.

Figure 3. (A) Crude 28-day mortality by admission lactate within each ratio quartile. (B) Interaction hazard ratios for the highest-versus-lowest contrast in one sample of 5,423 complete cases. A value below 1 means a smaller association in the named group. The cardiac surgical interval includes 1.

Figure 4. (A) Complete-case hazard ratio for the highest versus the lowest quartile. The blue curve lets all quartile contrasts change with log time. The purple curve lets the highest-quartile slope change at day 7. Shaded bands are 95% confidence intervals. Split-period estimates at day 3, day 7, and day 14 are given in Table 5. (B) Percent dying in the first 7 days, and percent dying over the next 21 days among patients alive at day 7.

# Main table

Table 1. MIMIC-IV patients with a valid AST/ALT ratio. Deaths are counted from the day of enzyme measurement.

|  | Overall | Q1 | Q2 | Q3 | Q4 |
| --- | --- | --- | --- | --- | --- |
| N | 16616 | 4157 | 4170 | 4136 | 4153 |
| Age, years | 68.9 (57.4-79.6) | 66.4 (54.7-77.1) | 69.3 (57.9-79.9) | 70.8 (60.0-80.8) | 69.5 (57.9-80.1) |
| Male | 9400 (56.6%) | 2542 (61.2%) | 2371 (56.9%) | 2221 (53.7%) | 2266 (54.6%) |
| SOFA | 6 (4-9) | 5 (3-8) | 6 (4-8) | 6 (4-9) | 7 (4-10) |
| Lactate, mmol/L | 1.9 (1.3-3.1) | 1.6 (1.2-2.6) | 1.8 (1.2-3.0) | 1.9 (1.3-3.3) | 2.2 (1.4-3.7) |
| Bilirubin, mg/dL | 0.6 (0.4-1.2) | 0.6 (0.4-1.3) | 0.6 (0.4-1.1) | 0.6 (0.4-1.1) | 0.7 (0.4-1.3) |
| Creatinine, mg/dL | 1.2 (0.8-2.0) | 1.1 (0.8-1.8) | 1.2 (0.8-1.9) | 1.2 (0.8-2.0) | 1.3 (0.9-2.3) |
| AST, U/L | 40 (23-92) | 35 (21-75) | 33 (21-75) | 38 (23-80) | 60 (33-145) |
| ALT, U/L | 28 (16-65) | 48 (27-111) | 27 (17-62) | 22 (13-47) | 21 (12-46) |
| AST/ALT | 1.4 (1.0-2.0) | 0.8 (0.6-0.9) | 1.2 (1.1-1.3) | 1.7 (1.6-1.8) | 2.7 (2.3-3.4) |
| 28-day death | 4249 (25.6%) | 825 (19.8%) | 938 (22.5%) | 1135 (27.4%) | 1351 (32.5%) |

Note. n = 16,616 after removal of 5 stays whose enzyme chart time followed the hospital death time. Death is within 28 days of the enzyme measurement. Continuous variables are median (interquartile range).

Table 2. Analysis cohort and eligible MIMIC-IV stays without a valid AST/ALT pair. Deaths are counted from ICU admission.

|  | Analysis cohort | No valid pair |
| --- | --- | --- |
| N | 16616 | 16780 |
| Age, years | 68.9 (57.4-79.6) | 70.0 (59.2-80.1) |
| Male | 9400 (56.6%) | 9591 (57.2%) |
| SOFA | 6 (4-9) | 4 (3-6) |
| Lactate, mmol/L | 1.9 (1.3-3.1) | 1.8 (1.3-2.7) |
| Creatinine, mg/dL | 1.2 (0.8-2.0) | 1.0 (0.7-1.5) |
| Cardiac surgical ICU | 1294 (7.8%) | 4760 (28.4%) |
| 28-day death from ICU admission | 4239 (25.5%) | 2285 (13.6%) |

Note. The first column is the analysis cohort of 16,616, after removal of 5 stays whose enzyme chart time followed hospital death. The second column is eligible stays with no valid enzyme pair. Death in both columns is counted from ICU admission.

Table 3. Adjusted associations.

| Contrast | Estimate | 95% CI | P |
| --- | --- | --- | --- |
| Q2 vs Q1, imputed | HR 1.04 | 0.95-1.14 | 0.41 |
| Q3 vs Q1, imputed | HR 1.17 | 1.07-1.29 | <0.001 |
| Q4 vs Q1, imputed | HR 1.43 | 1.31-1.57 | <0.001 |
| Q4 vs Q1, complete case | HR 1.36 | 1.23-1.51 | <0.001 |
| Q4 vs Q1, days 0-7 | HR 1.81 | 1.60-2.05 | <0.001 |
| Q4 vs Q1, days 7-28 | HR 1.13 | 0.99-1.28 | 0.065 |
| Days 7-28 relative to days 0-7 | HR 0.62 | 0.52-0.74 | <0.001 |
| Log AST per SD, imputed | HR 1.42 | 1.33-1.52 | <0.001 |
| Log ALT per SD, imputed | HR 0.77 | 0.72-0.82 | <0.001 |
| Q4 vs Q1, imputed, expanded covariates | HR 1.33 | 1.22-1.46 | <0.001 |
| Q4 vs Q1, excluding alcohol, hepatitis, fatty liver | HR 1.40 | 1.25-1.56 | <0.001 |
| eICU Q2 vs Q1, imputed, clustered SE | OR 1.15 | 0.99-1.35 | 0.068 |
| eICU Q3 vs Q1, imputed, clustered SE | OR 1.30 | 1.14-1.49 | <0.001 |
| eICU Q4 vs Q1, imputed, clustered SE | OR 1.95 | 1.69-2.26 | <0.001 |

Note. Quartile rows for 28-day death, and the log AST and log ALT rows, are multiply imputed on 16,616 patients. The day 0-7, day 7-28, and relative rows are from one time-varying model on the same 30 imputed datasets. The complete-case quartile row uses 10,890 patients. The complete-case joint enzyme estimates on these 10,890 patients are shown in Table 5. The eICU rows use all 12,374 patients and hospital-clustered standard errors. The complete-case eICU estimate is in Table 5.

Table 4. Interaction hazard ratios in one sample: 5,423 complete cases from the lowest and highest quartiles, with 1,629 deaths. Every row uses this sample.

| Contrast | HR | 95% CI |
| --- | --- | --- |
| Men versus women | 1.03 | 0.84-1.27 |
| Age >= 65 versus younger | 0.80 | 0.64-1.00 |
| Lactate >= 4 versus lower | 1.09 | 0.85-1.41 |
| Cardiac surgical ICU versus other units | 0.61 | 0.36-1.03 |

Note. N = 5,423 and deaths = 1,629 are the size of the Q1-plus-Q4 complete-case analysis. A value below 1 means a smaller association in the named group.

Table 5. Sensitivity analyses for timing, selection, enzyme bounds, and hospital death.

| Analysis | Contrast | Estimate | 95% CI |
| --- | --- | --- | --- |
| Alive at 24 hours, days 0-7 | Q4 vs Q1 | HR 1.82 | 1.59-2.08 |
| Alive at 24 hours, days 7-28 | Q4 vs Q1 | HR 1.17 | 1.03-1.32 |
| Landmark at 24 hours, until day 7 | Q4 vs Q1 | HR 1.82 | 1.58-2.09 |
| Landmark at 24 hours, days 7-28 | Q4 vs Q1 | HR 1.20 | 1.06-1.36 |
| Landmark at 24 hours, later relative to early | Q4 | HR 0.66 | 0.55-0.80 |
| Measurement-day deaths removed, days 0-7 | Q4 vs Q1 | HR 1.81 | 1.59-2.06 |
| Measurement-day deaths removed, days 7-28 | Q4 vs Q1 | HR 1.15 | 1.01-1.30 |
| Measurement-day deaths removed, later relative to early | Q4 | HR 0.63 | 0.53-0.76 |
| Log-time interaction, complete case, day 1 | Q4 vs Q1 | HR 1.95 | 1.61-2.35 |
| Log-time interaction, complete case, day 7 | Q4 vs Q1 | HR 1.28 | 1.15-1.43 |
| Log-time interaction, complete case, day 14 | Q4 vs Q1 | HR 1.11 | 0.96-1.27 |
| Log-time interaction, complete case, day 28 | Q4 vs Q1 | HR 0.95 | 0.79-1.15 |
| Piecewise-linear time, complete case, day 1 | Q4 vs Q1 | HR 1.65 | 1.43-1.91 |
| Piecewise-linear time, complete case, day 28 | Q4 vs Q1 | HR 1.28 | 0.97-1.69 |
| Days 0-3, imputed | Q4 vs Q1 | HR 2.04 | 1.72-2.41 |
| After day 3, imputed | Q4 vs Q1 | HR 1.25 | 1.12-1.38 |
| Days 0-14, imputed | Q4 vs Q1 | HR 1.62 | 1.46-1.79 |
| After day 14, imputed | Q4 vs Q1 | HR 0.99 | 0.83-1.19 |
| Inverse-probability weighted | Q4 vs Q1 | HR 1.46 | 1.32-1.61 |
| AST/ALT bounds 1-10,000, complete case | Q4 vs Q1 | HR 1.33 | 1.20-1.47 |
| AST/ALT bounds 10-1,000, complete case | Q4 vs Q1 | HR 1.41 | 1.27-1.58 |
| MIMIC-IV hospital death, imputed | Q4 vs Q1 | OR 1.49 | 1.32-1.67 |
| Log AST per SD, complete case | AST per SD | HR 1.38 | 1.28-1.48 |
| Log ALT per SD, complete case | ALT per SD | HR 0.81 | 0.75-0.87 |
| Log ratio per SD, imputed | Log AST minus log ALT | HR 1.16 | 1.13-1.20 |
| Log mean level per SD, imputed | Mean of log AST and log ALT | HR 1.09 | 1.06-1.13 |
| eICU complete case, hospital-clustered SE | Q4 vs Q1 | OR 1.96 | 1.63-2.36 |

Note. Quartile indicators were predictors in the imputation model. The rows labeled “Alive at 24 hours” exclude patients who died within 24 hours after ICU admission, and follow-up is still counted from the day of enzyme measurement. The landmark rows use those same patients but restart follow-up at 24 hours after ICU admission. Day 7 and day 28 are counted from ICU admission, so the early interval is about 6 days. Follow-up is censored at day 28 from ICU admission. Twenty deaths that were inside 28 days of the enzyme measurement fell after that censoring time. The measurement-day rows remove the 269 deaths on the calendar day of the enzyme pair from the 16,616-patient cohort and refit the day-7 split. The ratio row is the later hazard divided by the first-week hazard. The log ratio and log mean rows are the joint enzyme model rewritten in that form; the mean-level term is the test of whether the ratio alone is sufficient. The log-time rows are from one complete-case model (10,890 patients) in which the quartile coefficients change in proportion to log time. A likelihood-ratio test of those interactions gave P < 0.001. The piecewise-linear rows let the highest-quartile slope change at day 7. The day-3 and day-14 rows use the imputed time-varying model on the same patients. Inverse-probability weights used age, sex, SOFA, alcohol use disorder, and care unit, and were applied to all 16,616 analyzed stays; that model does not include lactate. Enzyme-bound rows keep the original quartile cutpoints. The hospital-death odds ratio uses the same covariates as the 28-day hazard model. The log AST and log ALT rows are the complete-case joint model (10,890 patients). The eICU row is the complete-case model with standard errors clustered by hospital (157 hospitals).

Table 6. Missing covariates by AST/ALT quartile in MIMIC-IV.

| Variable | Overall | Q1 | Q2 | Q3 | Q4 |
| --- | --- | --- | --- | --- | --- |
| Lactate missing | 31.9% | 40.0% | 33.1% | 29.9% | 24.8% |
| Bilirubin missing | 3.7% | 3.9% | 4.1% | 3.6% | 3.3% |
| Creatinine missing | 0.1% | 0.1% | 0.1% | 0.1% | 0.1% |
| SOFA missing | 0 | 0 | 0 | 0 | 0 |
| Albumin missing | 42.3% | 42.8% | 42.1% | 42.3% | 42.0% |

Note. Percent of the 16,616 analyzed stays. Lactate is the covariate that removes the most patients from a complete-case model, and it is missing more often in the lower quartiles. Alcohol use disorder, age, and sex had no missing values. Multiple imputation assumes that, given the variables in the imputation model, the missing values do not depend on the unobserved number.

Table 7. Standardized mean differences before and after inverse-probability weighting.

| Variable | Unweighted SMD | Weighted SMD |
| --- | --- | --- |
| Age | -0.04 | 0.01 |
| Male | -0.01 | -0.01 |
| SOFA | 0.26 | -0.01 |
| Alcohol use disorder | 0.06 | 0.00 |
| Medical ICU | 0.12 | 0.00 |
| Coronary care unit | 0.10 | 0.00 |
| Surgical or trauma ICU | 0.06 | -0.01 |
| Cardiac surgical ICU | -0.28 | 0.01 |
| Other unit | -0.01 | 0.00 |

Note. Each column compares the 16,616 analyzed stays with all 33,396 eligible stays. The denominator is the unweighted pooled standard deviation of the analyzed stays and the 16,780 stays with no valid pair, so weighting changes the mean but not the scale. Weights are stabilized and truncated at the 1st and 99th percentiles. Median weight 0.88 (interquartile range 0.73-1.07); truncated range 0.55-3.97; effective sample size 12,993. A positive value means a higher mean among analyzed stays. SOFA was complete.
