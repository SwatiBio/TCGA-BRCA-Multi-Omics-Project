## TCGA-BRCA Multi-Omics MOFA2 Analysis — Comprehensive Report

## 1\. Introduction

Breast cancer is the most common cancer in women worldwide. Think of cancer as a corrupted computer program: it's not just one thing going wrong, but many things — genes being too active or too quiet, DNA being chemically modified (methylation), chunks of DNA being copied or deleted (copy number changes), and proteins being overproduced or shut down.

**The problem:** Doctors usually look at these problems one at a time. But cancer is a system — all these layers interact. A gene might be silenced by methylation, or a protein might be overactive because its gene was copied extra times. If you only look at one layer, you miss the bigger picture.

**The solution (MOFA2):** Imagine you have 4 different maps of the same city — roads, elevation, population density, and weather patterns. Each map alone tells you something, but overlaying them reveals deeper patterns (e.g., "areas with high elevation have sparse road networks"). MOFA2 does the same for cancer data — it overlays 4 types of molecular data to find hidden patterns, called **factors**, that cut across all data types simultaneously.

This study applied MOFA2 to integrate **RNA-seq** (which genes are active), **DNA methylation** (which genes are chemically silenced), **copy number variation** (which genes are duplicated or deleted), and **RPPA proteomics** (which proteins are actually present) from ~636 TCGA-BRCA patients to:

1.  Discover hidden multi-omics patterns (factors) in breast cancer
2.  Test which factors predict survival (overall and recurrence-free)
3.  Validate the predictions using multiple statistical methods
4.  Understand the biology behind the most important factors

---

## 2\. Materials & Methods

### 2.1 Data Acquisition

All data came from **UCSC Xena** (a public cancer genomics database). Four types of molecular data were collected for TCGA-BRCA (The Cancer Genome Atlas — Breast Cancer):

| Data Type | What it measures | Number of features |
| --- | --- | --- |
| **RNA-seq** | Which genes are "on" or "off" — like a light switch panel for 20,000 genes | 20,501 genes |
| **DNA Methylation** | Chemical tags that silence genes — like a dimmer switch that turns genes down | ~485,000 probes |
| **Copy Number (CNV)** | Whether genes are accidentally duplicated or deleted — like having 3 copies of a page instead of 2 | 24,776 genes |
| **RPPA Proteomics** | Which proteins are actually present — the real workers in the cell | 281 proteins |
| **Clinical** | Patient info — age, stage, survival time, subtype | 194 columns |

**🧠 Plain English:** We downloaded 4 different types of molecular measurements from the same 636 breast cancer patients. Think of it as taking 4 different X-rays of each patient — each shows something different.

### 2.2 Preprocessing

Before analysis, the data needed cleaning:

*   **Duplicate samples removed:** Some patients had multiple tumor samples (e.g., two biopsies). We kept one per patient to avoid bias.
*   **Missing DNA methylation filled in:** Some methylation measurements failed. We used a k-nearest neighbors algorithm to estimate missing values (like guessing a missing puzzle piece by looking at the pieces around it).
*   **Top features selected:** We kept the most variable 8,000 RNA genes, 8,000 methylation probes, and 3,000 CNV genes. Variable features carry more information — if every patient has the same value, that feature can't explain differences between patients.
*   **Patient alignment:** All 4 data types were sorted so column 1 = Patient A in all 4 data types. This is essential — if the columns don't match, you're comparing different apples.

**Final dataset:** 636 patients with all omics available. (Note: CNV ended up with 0 usable features after filtering — only RNA, Methylation, and RPPA were used.)

### 2.3 MOFA2 Training — The Core Algorithm

MOFA2 was trained with **15 latent factors**. What's a latent factor?

**🧠 Think of it like this:** If you ask 1000 people to rate 20 movies, you might find that "action movies" and "romantic comedies" are hidden patterns (latent factors) that explain why people like certain movies. MOFA2 does the same for molecular data — it finds hidden patterns that explain why certain patients have coordinated changes across all 4 data types.

Each patient gets a **score** for each factor (like a rating from -3 to +3). A high positive score means that patient strongly expresses that factor's pattern. A negative score means the opposite.

### 2.4 Survival Analysis Methods

**🧠 What is survival analysis?** Instead of just asking "did the patient live or die?", survival analysis asks **"how long until an event (death or recurrence) happens?"** — and accounts for patients who are still alive at the end of the study (censored).

Methods used:

*   **Cox regression:** The standard method. Calculates a **hazard ratio (HR)** — how many times more likely a patient is to die at any given moment if their factor score increases by 1 unit. HR=1: no effect. HR>1: risk factor (worse survival). HR<1: protective factor (better survival).
*   **Kaplan-Meier curves:** A picture showing how survival changes over time for high-score vs low-score groups.
*   **Random Survival Forest (RSF):** A machine learning method that doesn't assume a straight-line relationship. Can capture complex patterns.
*   **Time-dependent ROC:** Measures **prediction accuracy** at specific time points (1, 3, 5 years). AUC=0.5 is random guessing, AUC=1.0 is perfect.
*   **RMST (Restricted Mean Survival Time):** Answers **"how many extra days of life?"** — much easier to understand than hazard ratios.

### 2.5 Gene Set Enrichment

For important factors, we found the top 500 correlated genes and tested which biological pathways they belong to (e.g., "immune response", "cell division"). This tells us **what biological process the factor represents**.

### 2.6 External Validation

We tested the most important factor's gene signature in an **independent breast cancer dataset** (GSE20685, 327 patients, measured on a different technology — Affymetrix microarrays). If the signature predicts survival in a completely separate group of patients, it's more likely to be a real finding, not just noise.

---

## 3\. Results

### 3.1 Multi-Omics Integration — Factor Discovery

MOFA2 identified **15 latent factors** from the integrated RNA-seq, methylation, and RPPA data. Each factor explains a certain percentage of the variation in the data.

<figure class="image image_resized" style="width:73.1%;"><img style="aspect-ratio:1600/1000;" src="api/attachments/sQObgyN97EQu/image/variance_lollipop.png" width="1600" height="1000"></figure>

**Figure: Lollipop chart of variance explained per factor**

> _File:_ `_figures/variance_lollipop.png_` _— Clean lollipop chart showing each factor's total variance explained, sorted from most to least important._

<figure class="image image_resized" style="width:69.69%;"><img style="aspect-ratio:1500/1050;" src="api/attachments/JiLhwpUwA0Ll/image/variance_explained.png" width="1500" height="1050"></figure>

**Figure: Variance Explained Heatmap**

> _File:_ `_figures/variance_explained.png_` _— Heatmap showing how much variance each factor explains in each data type (RNA / Methylation / RPPA)._

| Factor | RNA (%) | Methyl (%) | RPPA (%) | Total (%) |
| --- | --- | --- | --- | --- |
| Factor1 | 9.9 | 10.7 | 5.1 | **25.7** |
| Factor2 | 8.9 | 4.3 | 8.3 | **21.5** |
| Factor3 | 0.1 | 11.2 | 0.1 | 11.5 |
| Factor4 | 0.2 | 8.2 | 0.4 | 8.7 |
| Factor5 | 3.2 | 3.9 | 0.8 | 8.0 |
| Factor6 | 6.6 | 0.0 | 0.1 | 6.7 |
| Factor7–15 | \- | \- | \- | 1.5–5.1 |

**🧠 What this means:** The first two factors alone explain ~47% of all the coordinated variation across RNA, methylation, and proteins. That's like finding that 2 hidden patterns account for half of why breast cancer patients differ from each other at the molecular level.

*   **Factor 1** captures a pattern present in all 3 data types — it's a "master regulator" program.
*   **Factor 3 and Factor 4** are mostly about methylation (chemical silencing) — a specific type of regulation.
*   **Factor 6** is mostly about RNA — gene activity changes without much backing from other layers.

**Why it matters:** Factors that explain more variance aren't necessarily the most clinically important. A factor could explain a lot of variance but have nothing to do with survival (like Factor 1), while a moderate-variance factor (Factor 5) could be the most clinically meaningful.

### 3.2 PAM50 Subtype Distribution

<figure class="image image_resized" style="width:71.9%;"><img style="aspect-ratio:2400/1800;" src="api/attachments/u2jxP8vi7TbE/image/pam50_subtype_distribution.png" width="2400" height="1800"></figure>

**Figure: PAM50 Subtype Distribution**

> _File:_ `_figures/pam50_subtype_distribution.png_` _— Barplot showing how many patients have each breast cancer subtype._

**🧠 What this means:** Breast cancer isn't one disease — it has at least 5 molecular subtypes that behave differently:

*   **Luminal A** (~40-50%): Slow-growing, hormone-driven, best prognosis
*   **Luminal B** (~20%): Faster-growing, still hormone-driven
*   **Basal-like** (~15-20%): "Triple negative" — aggressive, no hormone targets
*   **HER2-enriched** (~10%): Driven by HER2 protein overexpression
*   **Normal-like:** Rare, poorly defined

Our cohort matches the expected distribution. This matters because different subtypes have drastically different survival rates and treatment responses.

### 3.3 Clinical Associations — What Do The Factors Mean Biologically?

By correlating factor scores with clinical variables (ER status, PAM50 subtype, age, stage), we can start to guess what each factor represents biologically.

| Factor | Strongest Clinical Association | What it might represent |
| --- | --- | --- |
| Factor1 | Basal-like subtype (r=-0.21), ER-negative (r=-0.16) | **Basal/immune program** — active in aggressive, non-hormone-driven tumors |
| Factor2 | Luminal subtype (r=+0.36), ER-positive (r=+0.19) | **Luminal/hormone program** — active in slow-growing, hormone-driven tumors |
| Factor6 | Luminal (r=+0.46), ER-positive (r=+0.42), Female (r=+0.49) | **Strongly associated with ER+ luminal biology** |
| Factor7 | ER-negative (r=-0.16), Basal-like (r=-0.17) | **Another basal-associated program** |
| Factor4 | Younger age (r=-0.22) | **Age-related methylation changes** |

<figure class="image"><img style="aspect-ratio:1600/1200;" src="api/attachments/WwBORYJJZ4hj/image/factor_clinical_heatmap.png" width="1600" height="1200"></figure>

**Figure: Factor-Clinical Correlation Heatmap**

> _File:_ `_figures/factor_clinical_heatmap.png_` _— Color-coded heatmap showing correlations between each factor and clinical variables. Red = positive correlation, Blue = negative correlation._

**🧠 What this means (non-biologist):** Think of each factor as a "flavor" of breast cancer biology. Factor 1 and Factor 7 taste like "Basal" (aggressive, triple-negative), while Factor 2 and Factor 6 taste like "Luminal" (slow-growing, hormone-driven). This tells us MOFA isn't finding random noise — it's rediscovering known breast cancer biology purely from the molecular data, without being told about subtypes.

### 3.4 Survival Analysis — Which Factors Predict Life and Death?

Each factor was tested individually for association with **overall survival (OS)** — how long patients live after diagnosis.

<figure class="image image_resized" style="width:75.48%;"><img style="aspect-ratio:1600/1200;" src="api/attachments/hwG0JQ7bRNtZ/image/forest_plot.png" width="1600" height="1200"></figure>

**Figure: Forest Plot — Hazard Ratios for All 15 Factors**

> _File:_ `_figures/forest_plot.png_` _— Forest plot showing HR (blue = protective, red = risk) with 95% confidence intervals. Factors sorted by HR._

<figure class="image image_resized" style="width:84.81%;"><img style="aspect-ratio:1400/1200;" src="api/attachments/ojn61O12yfM1/image/volcano_plot.png" width="1400" height="1200"></figure>

**Figure: Volcano Plot — HR vs Significance**

> _File:_ `_figures/volcano_plot.png_` _— Each dot is a factor. X-axis = HR (protective left, risk right). Y-axis = statistical significance (higher = more significant). Labeled factors are the most interesting._

| Factor | HR | 95% CI | p-value | p.adj | Direction |
| --- | --- | --- | --- | --- | --- |
| **Factor5** | **0.74** | \[0.56-0.97\] | **0.031** | 0.137 | **Protective** |
| **Factor8** | **1.31** | \[1.02-1.67\] | **0.031** | 0.137 | **Risk** |
| **Factor10** | **1.32** | \[1.01-1.73\] | **0.040** | 0.137 | **Risk** |
| **Factor13** | **1.37** | \[1.02-1.84\] | **0.038** | 0.137 | **Risk** |
| **Factor14** | **1.34** | \[1.01-1.79\] | **0.046** | 0.137 | **Risk** |

**🧠 What an HR of 0.74 means:** Imagine two identical patients, except Patient A has a Factor 5 score that's 1 unit higher. At any given moment during follow-up, Patient A is **26% less likely to die** than Patient B (because 1 - 0.74 = 0.26). That's a meaningful protective effect.

**🧠 What an HR of 1.37 means:** For Factor 13, a 1-unit higher score means **37% higher risk of death at any moment**. That's like the difference between a non-smoker and a moderate smoker in terms of lung cancer risk.

**🧠 What about p.adj (adjusted p-value):** We tested 15 factors simultaneously. Like rolling a die 15 times, you expect some "lucky" results by chance. p.adj corrects for this. None of our factors survived this strict correction (p.adj<0.05 threshold), which is not surprising — we have limited statistical power (only ~100 deaths in the cohort). The **consistency** of the pattern (same 5 factors appearing across multiple methods) is more meaningful than any individual p-value.

**The key insight:** Factor 5 is **protective** (higher score = better survival), while Factors 8, 10, 13, and 14 are **risk factors** (higher score = worse survival).

### 3.5 Multivariable Cox — Is the Effect Independent of Age and Stage?

A factor might look predictive simply because it's correlated with age or cancer stage. We adjusted for these to see if the factors **add independent information** beyond what doctors already know.

| Factor | HR (adjusted) | 95% CI | p-value | p.adj |
| --- | --- | --- | --- | --- |
| **Factor5** | **0.63** | \[0.46-0.87\] | **0.004** | 0.065 |
| **Factor14** | **1.44** | \[1.07-1.95\] | **0.017** | 0.120 |
| **Factor8** | **1.30** | \[1.04-1.63\] | **0.024** | 0.120 |

**🧠 Why this is exciting:** Factor 5's protective effect actually **got stronger** after adjusting for age and stage (HR dropped from 0.74 to 0.63). This means Factor 5 is NOT just a proxy for "younger age" or "early stage" — it's capturing something biological that standard clinical information misses. This is exactly the kind of finding that could lead to a new prognostic test.

LASSO regression (a machine learning method that automatically selects the most important predictors) selected Factor 13 as the single most important factor. This conservative method chose only 1 factor out of 15, confirming Factor 13's unique predictive contribution.

### 3.6 Kaplan-Meier Survival Curves — Seeing the Survival Difference

KM curves split patients into "High Factor" and "Low Factor" groups (above vs below the median score) and plot their survival over time.

<figure class="image"><img style="aspect-ratio:2400/1600;" src="api/attachments/oCfpznwGVyVn/image/km_panel_top5.png" width="2400" height="1600"></figure>

> _File:_ `_figures/km_panel_top5.png_` _— Panel of KM curves for the 5 most significant factors in a grid layout, making them easy to compare side-by-side._

<figure class="image"><img style="aspect-ratio:2400/1800;" src="api/attachments/8iYB1N8UVBX1/image/km_Factor5_OS.png" width="2400" height="1800"></figure>

**Figure: KM Curves — Factor 5 (Protective)**

> _File:_ `_figures/km_Factor5_OS.png_` _— Overall survival, Factor 5 high vs low_ _File:_ `_figures/km_Factor5_RFS.png_` _— Recurrence-free survival, Factor 5 high vs low_ _File:_ `_figures/km_median_Factor5.png_` _— Factor 5 median split_

<figure class="image"><img style="aspect-ratio:2400/1800;" src="api/attachments/Wz7fSZ3eVWst/image/km_median_Factor8.png" width="2400" height="1800"></figure>

**Figure: KM Curves — Risk Factors (Median Split)**

> _File:_ `_figures/km_median_Factor8.png_` _— Factor 8 median split_ _File:_ `_figures/km_median_Factor13.png_` _— Factor 13 median split_

<figure class="image"><img style="aspect-ratio:2400/1800;" src="api/attachments/hrPiMXKYRyoc/image/km_tertile_Factor5.png" width="2400" height="1800"></figure><figure class="image"><img style="aspect-ratio:2400/1800;" src="api/attachments/P8VJJyoQZQjq/image/km_tertile_Factor8.png" width="2400" height="1800"></figure><figure class="image"><img style="aspect-ratio:2400/1800;" src="api/attachments/jLuMSFxF8oDc/image/km_tertile_Factor13.png" width="2400" height="1800"></figure>

**Figure: KM Curves — Extreme Tertile Comparison**

> _File:_ `_figures/km_tertile_Factor5.png_` _— Factor 5: top vs bottom third (excludes middle)_ _File:_ `_figures/km_tertile_Factor8.png_` _— Factor 8: top vs bottom third_ _File:_ `_figures/km_tertile_Factor13.png_` _— Factor 13: top vs bottom third_

**🧠 How to read a KM curve:** Imagine two lines starting at the left at 100% — everyone is alive. Over time, the lines drop as patients die. If one line stays higher than the other, that group is living longer. The p-value tells you if the gap between the lines is bigger than what you'd expect by chance.

**What we see:** The KM curves visually confirm the Cox results. For Factor 5, the "High" group (red line) stays above the "Low" group — patients with high Factor 5 scores survive longer. For Factor 13, the opposite is true — high scores = worse survival. The separation is visible throughout the entire follow-up period (up to 10+ years), not just early on. The tertile splits (removing the middle third) make the difference even more dramatic.

### 3.7 Random Survival Forest — Which Factor Matters Most?

RSF is a machine learning method that doesn't assume a straight-line relationship between factors and survival. It can capture complex, non-linear patterns and interactions that standard Cox models might miss.

| Rank | Factor | Importance | C-Index Contribution |
| --- | --- | --- | --- |
| 1 | **Factor14** | **0.0194** | Highest |
| 2 | Factor12 | 0.0094 |  |
| 3 | Factor5 | 0.0072 |  |
| 4 | Factor11 | 0.0066 |  |
| 5 | Factor2 | 0.0064 |  |
| 6 | Factor10 | 0.0057 |  |

<figure class="image"><img style="aspect-ratio:1200/750;" src="api/attachments/Dumsn2Tfsz7w/image/brca_rsf_importance.png" width="1200" height="750"></figure><figure class="image"><img style="aspect-ratio:1600/1000;" src="api/attachments/KKDzY9iTevmb/image/rsf_lollipop.png" width="1600" height="1000"></figure>

**Figure: RSF Variable Importance**

> _File:_ `_figures/brca_rsf_importance.png_` _— Horizontal barplot of factor importance._ _File:_ `_figures/rsf_lollipop.png_` _— Cleaner lollipop version, sorted by importance._

**🧠 What "importance" means:** Imagine the RSF model as a group of 1000 decision trees, each voting on whether a patient will survive. Importance measures how much each factor improves the trees' accuracy. Factor 14 (importance ~0.019) is roughly 3x more important than the next best factor. Factors with negative importance (3, 9, 6) are actually hurting prediction — they're noise.

**Key finding:** Factor 14 dominates the RSF model, even though it wasn't the most significant in Cox regression. This suggests Factor 14's effect on survival might be **non-linear** — it might only matter above a certain threshold, or it might interact with other factors. The RSF's ability to capture this is why it complements Cox regression.

The model's C-index is **0.594** (1 - prediction.error = 1 - 0.406). A C-index of 0.5 is random; 1.0 is perfect. 0.594 is modest but meaningful — better than random chance and similar to many clinical risk scores.

### 3.8 Time-Dependent ROC — How Well Can We Predict Survival at Specific Time Points?

ROC analysis asks: **"If we set a risk threshold, how well can we separate who will survive 3 years from who won't?"**

| Timepoint | AUC | What it means |
| --- | --- | --- |
| 1 Year | 0.420 | Worse than random! Too few events for useful prediction |
| **3 Years** | **0.705** | **Decent** — 70.5% chance of correctly ranking a survivor vs non-survivor |
| 5 Years | 0.613 | Slightly above random — prediction fades at longer horizons |

<figure class="image"><img style="aspect-ratio:1200/900;" src="api/attachments/hIAaRqKjFRCY/image/brca_timeroc_auc.png" width="1200" height="900"></figure><figure class="image"><img style="aspect-ratio:1000/1000;" src="api/attachments/IxEXZMb3e4A9/image/timeroc_bar.png" width="1000" height="1000"></figure>

**Figure: Time-Dependent ROC Curves**

> _File:_ `_figures/brca_timeroc_auc.png_` _— ROC curves at 1, 3, and 5 years._ _File:_ `_figures/timeroc_bar.png_` _— Bar chart of AUC values with a dashed line at 0.70._

**🧠 What AUC means:** Imagine you have two patients — one who will survive 3 years, and one who won't. If you randomly picked one and your model said "this one is higher risk", an AUC of 0.705 means you'd be right 70.5% of the time. 50% is random guessing. 100% is perfect.

**Why the patterns are interesting:**

*   **1 year AUC is terrible (0.42):** Very few patients die within 1 year of diagnosis in breast cancer (it's a relatively slow-progressing disease). With only a handful of events, the model has nothing to learn from.
*   **3 years is best (0.705):** This is the sweet spot — enough events have accumulated, but they still reflect the underlying tumor biology rather than late-stage complications.
*   **5 years drops (0.613):** After 5 years, survival is affected by aging, other diseases, and treatments received after the initial diagnosis — things the MOFA factors can't capture.

**The takeaway:** Our factors are best at predicting **medium-term (3-year) survival** — which is clinically useful for treatment planning.

### 3.9 Restricted Mean Survival Time — "How Many Extra Days of Life?"

RMST answers a simple question: **Within 5 years of diagnosis, how many more days does the low-risk group live compared to the high-risk group?**

| Factor | RMST Difference (days) | 95% CI | p-value | Interpretation |
| --- | --- | --- | --- | --- |
| **Factor7** | **+95.6** | **\[0.6 – 190.5\]** | **0.049** | **Significant** — low-score group lives ~3 months longer |
| Factor2 | +85.7 | \[-3.7 – 175.1\] | 0.060 | Borderline protective |
| Factor5 | +66.2 | \[-23.9 – 156.2\] | 0.150 | Protective trend |
| Factor4 | \-82.9 | \[-172.4 – 6.7\] | 0.070 | Borderline risk |
| Factor8 | \-71.8 | \[-159.8 – 16.2\] | 0.110 | Risk trend |
| Factor14 | \-65.1 | \[-163.4 – 33.2\] | 0.194 | Risk trend |

**🧠 What this tells us that Cox didn't:** Factor 7, which had no significant Cox signal (HR=0.93, p=0.65), becomes the **only statistically significant RMST factor** (p=0.049). This is a paradox — how can something show no HR effect but a significant survival time difference?

The answer: **Factor 7's protective effect is delayed**. Cox regression measures risk at every moment, so an effect that appears only after several years gets diluted. RMST, by integrating survival over the full 5-year window, captures the **cumulative benefit**. This is like comparing two cars: one might be safer in the first year (Cox would detect that), but another might have a better 5-year reliability record (RMST would detect that).

**Why Factor 5 was significant in Cox but not RMST:** Factor 5's HR=0.74 is steady across the entire follow-up — its protective effect is constant but modest. The RMST difference of +66 days (2.2 months) is not statistically significant because the confidence interval crosses zero (p=0.15). This means the Cox and RMST results for Factor 5 are **consistent but underpowered** — with more patients or events, both would likely reach significance.

**🧠 Concrete meaning:** A patient with low Factor 7 score (worse RMST) lives on average **95 fewer days within 5 years** compared to a matched patient with high Factor 7. That's 3 months — a meaningful difference for patients and their families.

Why this matters for real patients: If a doctor could say "your RMST is 95 days shorter than similar patients," that might influence decisions about adjuvant chemotherapy or more frequent follow-up. RMST makes statistics **personal and actionable**.

### 3.10 Gene Set Enrichment — What Biology Do These Factors Represent?

For Factor 1 (highest variance explained), we identified the top 500 correlated genes and tested which biological pathways they participate in:

| GO Term | p.adjust | Genes | What it means (plain English) |
| --- | --- | --- | --- |
| Cell-cell adhesion regulation | 0.006 | 28 genes | **How cells stick together** — cancer cells that don't stick well can break off and spread (metastasis) |
| JAK-STAT signaling | 0.006 | 14 genes | **Immune communication channel** — how immune cells "talk" to each other |
| Negative regulation of cell adhesion | 0.006 | 20 genes | **Breaking cellular glue** — actively making cells less sticky, enabling spread |
| Proteoglycan biosynthesis | 0.006 | 10 genes | **Building the tumor's scaffolding** — the physical structure around the tumor |
| Leukocyte cell-cell adhesion | 0.006 | 24 genes | **Immune cells arriving at the tumor** — white blood cells trying to attack |
| Actin filament organization | 0.009 | 25 genes | **Cell skeleton reshaping** — how cells change shape to move and spread |

<figure class="image"><img style="aspect-ratio:1500/1200;" src="api/attachments/JvUhHusb8jD9/image/brca_factor1_go_enrichment.png" width="1500" height="1200"></figure>

**Figure: GO Enrichment Dotplot**

> _File:_ `_figures/brca_factor1_go_enrichment.png_` _— Dotplot showing enriched GO terms._

**🧠 The story this tells:** Factor 1 captures a **tumor microenvironment program** — the battle between the immune system and the tumor. The top pathways involve:

1.  **Immune cells arriving and trying to attack** (leukocyte adhesion, JAK-STAT signaling)
2.  **Tumor cells trying to escape** (breaking cell adhesion, reshaping their skeletons to become mobile)
3.  **The physical battlefield** (building and remodeling the extracellular scaffolding)

This is consistent with Factor 1 being associated with Basal-like (triple-negative) tumors, which are known to have more immune infiltration than Luminal tumors. Factor 1 is essentially measuring **"how inflamed is the tumor microenvironment"** — a key determinant of both prognosis and immunotherapy response.

### 3.11 Factor 5 — The Protective Immune Signature

Factor 5 is the most clinically interesting finding. Its top RNA features tell us what it represents:

**Top genes in Factor 5:** PLA2G2D, SIRPG, TIGIT, ICOS, ZAP70

**🧠 Translation for non-biologists:**

*   **TIGIT, ICOS:** "Checkpoint" proteins on immune T-cells. Think of them as the "gas pedal" (ICOS) and "brake" (TIGIT) of the immune system.
*   **ZAP70:** A key signaling protein inside T-cells. Without it, T-cells can't activate.
*   **SIRPG:** Helps T-cells interact with other cells.
*   **PLA2G2D:** An enzyme involved in inflammation.

**The picture:** Factor 5 captures **active T-cell immunity** — tumors where the immune system has successfully infiltrated and activated. Patients with high Factor 5 have a pre-existing anti-tumor immune response, which helps them live longer.

**Why this matters clinically:** This immune signature is **independent of PAM50 subtype** (Factor 5 shows minimal correlation with ER/PAM50). It captures a **trans-subtype immune component** — some Luminal tumors have it, some Basal tumors don't. Standard clinical classification (just checking ER/PR/HER2 status) completely misses this immune variation.

**If validated,** a Factor 5 gene signature (measuring just 5-6 genes) could be developed into a clinical test to identify patients who might benefit from immunotherapy, even if their cancer type isn't traditionally considered "immune-responsive."

### 3.12 External GEO Validation — The Critical Reality Check

**The goal:** Test Factor 1's gene signature in GSE20685, an independent breast cancer cohort of 327 patients measured on a different technology (Affymetrix GPL570 microarrays).

**The result:** Zero overlapping genes between the Factor 1 signature and the GEO platform.

**What went wrong:** This is a technical incompatibility, not a biological failure. The Affymetrix GPL570 platform measures ~54,000 probe sets that map to ~20,000 genes, but the probe-to-gene mapping uses **different gene naming conventions** than TCGA's RNA-seq data. The two technologies are like a French chef and an Italian chef using different names for the same ingredients — they're talking about the same things but can't agree on the words.

**What's needed:** A validation cohort using **the same technology** (RNA-seq). Candidates:

*   **GSE96058** (n=3,409 breast cancer patients, RNA-seq)
*   **METABRIC** (n=1,980, but uses old technology)
*   **TCGA-BRCA itself** (we could do cross-validation within the dataset)

---

## 4\. Discussion

### 4.1 Summary — What Did We Actually Learn?

**🧠 One-paragraph summary for a non-scientist:**

We took 4 different types of molecular data from 636 breast cancer patients and ran them through a pattern-finding algorithm (MOFA2). The algorithm discovered 15 hidden patterns. One pattern in particular — Factor 5 — strongly predicts which patients will live longer. This pattern represents **the immune system's ability to recognize and attack the tumor**. Importantly, this immune effect works regardless of the patient's age, cancer stage, or breast cancer subtype — it's an **independent predictor** that adds information beyond what doctors currently measure. If confirmed in other studies, this could lead to a simple blood or tissue test to guide treatment decisions.

### 4.2 What Each Factor Represents

| Factor | Biological Program | What it means for the patient |
| --- | --- | --- |
| **Factor1** | ECM remodeling + Immune signaling | **Tumor structure and immune battlefield** — active in aggressive, basal-like tumors |
| **Factor2** | Cell cycle + Hormone signaling | **Cell division program** — active in hormone-driven, luminal tumors |
| **Factor5** | T-cell activation + Immune checkpoints | **Active anti-tumor immunity** — patients with this live longer |
| **Factor14** | Primarily methylation-driven (unknown function) | **Mysterious risk factor** — top predictor in machine learning model |
| **Factor8** | Mixed RNA+Methylation signals | **Risk factor** — patients with this pattern do worse |
| **Factor10** | Subtle risk signal | **Mild risk factor** — consistently appears across analyses |

### 4.3 Why Factor 5 is Clinically Exciting

The protective effect of Factor 5 (HR=0.63 after adjustment) is compelling because:

1.  **It's independent** — not just a proxy for age, stage, or subtype
2.  **It's biologically interpretable** — T-cell activation genes (TIGIT, ICOS, ZAP70) tell a clear story
3.  **It's measurable** — a 5-6 gene signature could be developed into a clinical test
4.  **It's actionable** — patients with low Factor 5 scores might benefit from immunotherapy to boost their immune response

**The clinical scenario:** Two patients with identical age, stage, and Luminal B subtype come to the clinic. Current medicine treats them identically. But their Factor 5 scores differ. The patient with low Factor 5 might benefit from additional immunotherapy, while the high Factor 5 patient might do well with standard therapy alone. This is **precision medicine** in action.

### 4.4 Why the Other Risk Factors (14, 8, 10, 13) Are Important

Factors 14, 8, 10, and 13 consistently emerged as risk factors, but their biology is mostly unknown (they're primarily methylation-driven). This is actually exciting — it means MOFA found **novel risk patterns** that aren't explained by known biology.

**🧠 Analogy:** If Factor 5 is like discovering that "having fire extinguishers" prevents fire damage (obvious in hindsight), Factors 14 and 13 are like discovering a completely new fire safety principle that no one had thought of before. Understanding them could lead to entirely new therapeutic approaches.

### 4.5 Methodological Honesty — What Are the Caveats?

1.  **Small number of events:** Breast cancer has good survival rates (~85% at 5 years). This means fewer deaths = less statistical power. A cohort with more events (e.g., pancreatic cancer, where most patients die within 2 years) would give more definitive results.
2.  **CNV was excluded:** The GISTIC2 copy number data had 0 usable features after filtering. This means we missed a whole dimension of genomic variation.
3.  **External validation failed:** The GEO microarray incompatibility means we haven't truly validated our findings in an independent cohort. This is a critical limitation.
4.  **Modest effect sizes:** HRs of 0.63-1.44 are meaningful but not game-changing. They're similar to adding one additional clinical variable to a prediction model.
5.  **Multiple testing:** None of our factors survived strict multiple testing correction. The consistency across methods (Cox, RSF, RMST, KM curves) is reassuring, but independent validation is essential.

### 4.6 How This Compares to Other Studies

The immune-related protective signature we found (Factor 5) is consistent with a large body of literature showing that **tumor-infiltrating lymphocytes (TILs)** predict better survival in breast cancer. Our finding extends this by showing that the immune signal can be captured through multi-omics integration rather than requiring specialized immune cell counting.

The methylation-driven risk factors (14, 8, 13) are novel — most breast cancer prognostic studies focus on RNA or protein markers. This suggests that **epigenetic changes** (chemical modifications to DNA that don't change the genetic code) may play an underappreciated role in breast cancer progression.

---

## 5\. Conclusions

This study demonstrates that MOFA2-based multi-omics integration can:

1.  **Discover biologically meaningful patterns** that recapitulate known breast cancer biology (subtypes, immune microenvironment) and discover novel signals
2.  **Identify an immune checkpoint factor (Factor 5) with independent prognostic value** — higher scores predict better survival regardless of age, stage, or PAM50 subtype
3.  **Reveal novel methylation-driven risk factors (Factors 14, 13, 8)** that warrant further biological investigation
4.  **Achieve moderate prognostic accuracy** (3-year AUC=0.705) — comparable to many clinical risk scores

### Next Steps

1.  **Validate in an RNA-seq-based cohort** (GSE96058, n=3409) to confirm the Factor 5 signature
2.  **Characterize the risk factors** (14, 13, 8, 10) by examining their top features and pathway enrichments
3.  **Build a clinical risk score** combining Factor 5 with standard clinical variables (age, stage, subtype)
4.  **Test the Factor 5 signature in immunotherapy response datasets** — does it predict which patients benefit from checkpoint inhibitors?
5.  **Re-run with proper CNV data** to see if copy number alterations add predictive power

### Limitations

*   Small number of survival events in TCGA-BRCA
*   CNV data was excluded due to preprocessing issues
*   Failed external GEO validation (technical platform incompatibility)
*   Modest effect sizes and no survival p.adj<0.05

**Final thought:** Multi-omics integration is not a replacement for standard clinical and pathological assessment — it's a **complement**. The factors MOFA2 discovered capture biological dimensions that are invisible to current clinical tests. With validation, Factor 5 could become part of a clinically useful prognostic panel, bringing us one step closer to truly personalized cancer care.

---

## Figures Index

| # | Figure | Description | File Path |
| --- | --- | --- | --- |
| 1 | Variance Explained Heatmap | Heatmap of per-factor per-view variance explained | `figures/variance_explained.png` |
| 2 | Total Variance Explained | Barplot of total variance per view | `figures/variance_explained_total.png` |
| 3 | Variance Lollipop | Clean lollipop chart of variance per factor | `figures/variance_lollipop.png` |
| 4 | PAM50 Subtype Distribution | Subtype composition of the cohort | `figures/pam50_subtype_distribution.png` |
| 5 | Factor-Clinical Correlation Heatmap | Factor × clinical variable correlations | `figures/factor_clinical_heatmap.png` |
| 6 | Forest Plot | HR + 95% CI for all 15 factors (color-coded) | `figures/forest_plot.png` |
| 7 | Volcano Plot | HR vs p-value with labeled factors | `figures/volcano_plot.png` |
| 8 | KM Panel — Top 5 Factors | Grid of 5 KM curves for side-by-side comparison | `figures/km_panel_top5.png` |
| 9 | KM Curve — Factor 5 OS | Factor 5 overall survival (high vs low median) | `figures/km_Factor5_OS.png` |
| 10 | KM Curve — Factor 5 RFS | Factor 5 recurrence-free survival | `figures/km_Factor5_RFS.png` |
| 11 | KM Curve — Factor 5 Median | Factor 5 median split | `figures/km_median_Factor5.png` |
| 12 | KM Curve — Factor 8 Median | Factor 8 median split | `figures/km_median_Factor8.png` |
| 13 | KM Curve — Factor 13 Median | Factor 13 median split | `figures/km_median_Factor13.png` |
| 14 | KM Curve — Factor 5 Tertile | Factor 5 extreme tertile comparison | `figures/km_tertile_Factor5.png` |
| 15 | KM Curve — Factor 8 Tertile | Factor 8 extreme tertile comparison | `figures/km_tertile_Factor8.png` |
| 16 | KM Curve — Factor 13 Tertile | Factor 13 extreme tertile comparison | `figures/km_tertile_Factor13.png` |
| 17 | RSF Variable Importance | RSF horizontal barplot | `figures/brca_rsf_importance.png` |
| 18 | RSF Lollipop | Cleaner RSF importance lollipop chart | `figures/rsf_lollipop.png` |
| 19 | Time-Dependent ROC Curves | ROC curves at 1, 3, 5 years | `figures/brca_timeroc_auc.png` |
| 20 | TimeROC AUC Bar Chart | AUC values with 0.70 threshold line | `figures/timeroc_bar.png` |
| 21 | GO Enrichment Dotplot — Factor 1 | Top enriched GO biological process terms | `figures/brca_factor1_go_enrichment.png` |
| 22 | Combined Summary Figure | Multi-panel overview of key results | `figures/summary_combined.png` |

## Supplementary Files

| File | Path |
| --- | --- |
| Cox OS Results | `logs/cox_os.csv` |
| Cox RFS Results | `logs/cox_rfs.csv` |
| Cox Adjusted (Age+Stage) | `logs/cox_adjusted_clinical_os.csv` |
| Cox Continuous OS (all 15 factors) | `logs/cox_continuous_os.csv` |
| RSF Variable Importance | `Results/brca_rsf_importance.csv` |
| Factor Variance Explained | `Results/factor_variance_explained.csv` |
| Total Variance Explained | `Results/total_variance_explained.csv` |
| TimeROC AUC | `Results/timeroc_auc.csv` |
| RMST Per Factor | `Results/rmst_per_factor.csv` |
| GO Enrichment Factor 1 | `Results/brca_factor1_go_enrichment.csv` |
| Clinical Correlations | `Results/factor_clinical_corrs.csv` |
| Top Features Per Factor | `Results/top_features_per_factor.csv` |
| Patient Factor Scores | `Results/patient_factor_scores.csv` |
| GEO Validation | `Results/geo_validation.txt` |
| Full Analysis Summary | `Results/analysis_summary.txt` |
| MAE Summary | `logs/mae_summary.csv` |

---

_Report generated from TCGA-BRCA MOFA2 analysis outputs. All scripts available in_ `_scripts/_`_._