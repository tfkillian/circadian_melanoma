# circadian_melanoma

This repository contains a literature review of the role of circadian
rhythm-related genes in melanoma and an interrogation of these genes in the Depmap cancer data

### Relevance of circadian rhythm in the context of cancer

The circadian rhythm, sometimes referred to as the internal biological clock,
coordinates a vast array of physiological processes over a 24-hour cycle,
including hormone release, metabolism, and cell division. In recent years,
accumulating evidence has linked disruptions in circadian rhythms to the
development and progression of various cancers. When the normal timing of key
molecular clocks is altered—due to factors such as shift work, sleep deprivation,
or genetic mutations—cells may experience dysregulated gene expression, impaired
DNA repair mechanisms, and abnormal cell-cycle progression. These changes can
facilitate oncogenic transformation, promote tumor growth, and even influence
the efficacy of treatments. As a result, understanding the intricate relationship
between circadian regulation and cancer offers new avenues for both prevention
strategies and more precisely timed therapeutic interventions.

Emerging research has identified several circadian rhythm-related genes implicated in melanoma, influencing tumor progression, immune response, and patient outcomes. Key findings include:

    BMAL1 (ARNTL): Studies have shown that BMAL1 expression is significantly reduced in melanoma tissues compared to normal skin, indicating a disrupted circadian clock. Higher BMAL1 expression levels correlate with longer overall survival in melanoma patients, suggesting its potential role as a favorable prognostic marker.

    RORA (Retinoic Acid Receptor-Related Orphan Receptor-α): RORA expression is downregulated in melanoma patients, and higher RORA levels are associated with better prognosis following immunotherapy. This suggests that RORA may enhance anti-tumor immunity and could be a potential therapeutic target.

    PER1 and PER2 (Period Circadian Regulator 1 and 2): These genes are integral components of the circadian clock. Reduced expression of PER1 and PER2 has been observed in various cancers, including melanoma, and is associated with tumor progression. PER1 overexpression can induce apoptosis in cancer cells, while PER2 downregulation is linked to increased tumor growth.

    NR1D2 (Nuclear Receptor Subfamily 1 Group D Member 2): Higher expression levels of NR1D2 are associated with longer survival among skin cutaneous melanoma patients, indicating its potential protective role.

"RORA", "PER3", "PER1", "CRY2", "RORB", "PER2", "NR1D2", 
"RORC", "CRY1", "ARNTL", "NPAS2", "CLOCK", "NR1D1", "ARNTL2"

These findings underscore the significance of circadian rhythm-related genes in melanoma biology. Further research is warranted to explore their potential as biomarkers and therapeutic targets in melanoma treatment.

## Analyses

A RMarkdown HTML report details an exploration of the nature of circadian genes in melanoma
cell lines from the Broad Institute [Depmap](https://depmap.org/portal/data_page/?tab=allData)
datasets using the [depmap R package](https://www.bioconductor.org/packages/release/data/experiment/html/depmap.html),
cancer dependency data described by [Tsherniak, Aviad, et al. "Defining a cancer dependency map." Cell 170.3 (2017): 564-576.](https://www.ncbi.nlm.nih.gov/pubmed/28753430).

## Literature

[githib Rhythm_analyze](https://github.com/weijinchen01/Rhythm_analyze)

[Shafi AA, Knudsen KE. Cancer and the circadian clock. Cancer Res 2019;79:3806–14]()

[The Circadian Clock Component RORA Increases Immunosurveillance in Melanoma by Inhibiting PD-L1 Expression](https://pmc.ncbi.nlm.nih.gov/articles/PMC11247325/)

[A cellular hierarchy in melanoma uncouples growth and metastasis](https://www.nature.com/articles/s41586-022-05242-7)

[Decoding the interplay between genetic and non-genetic drivers of metastasis](https://www.nature.com/articles/s41586-024-07302-6)

[Secreted Apoe rewires melanoma cell state vulnerability to ferroptosis](https://www.science.org/doi/full/10.1126/sciadv.adp6164)
