#  APDI: Assessing the Performance of Dietary Assessment Instruments
Food Frequency Questionnaires (FFQs) are extensively utilized in nutritional epidemiology to gauge dietary intake over specified periods. However, FFQs often encounter measurement errors. 

Given the availability of a reference instrument (F), APDI evaluates the performance of dietary assessment  instrument (Q) like FFQs using two common metrics:
1. The correlation (ρ) between true intake (T) and the assessment instrument (Q).
2. The attenuation factor (λ).

Building on the measurement error model by Thompson et al. (2008), our R package APDI employs the method of moments for estimation. This package estimates:
- The correlation between nutrient intake (Q_N) and true nutrient intake (T_N) for each nutrient (N).
- The correlation between energy intake (Q_E) and true energy intake (T_E).   
- The correlation between nutrient residuals (Q_R) and true residuals (T_R).
- The attenuation factors for Q_N, Q_E, and Q_R.

# Usage Notes:
1. We recommend transforming each variable (e.g., using a log transformation) to approximate a normal distribution before applying our functions.
2. Our calculations assume there is no missing data in the dataset.
3. Following Thompson et al. (2008), nutrients and energy are measured twice for each Q and F.
   
Reference: Thompson FE, Kipnis V, Midthune D, Freedman LS, Carroll RJ, Subar AF, Brown CC, Butcher MS, Mouw T, Leitzmann M, Schatzkin A. Performance of a food-frequency questionnaire in the US NIH-AARP (National Institutes of Health-American Association of Retired Persons) Diet and Health Study. Public Health Nutr. 2008 Feb;11(2):183-95. doi: 10.1017/S1368980007000419.
   
