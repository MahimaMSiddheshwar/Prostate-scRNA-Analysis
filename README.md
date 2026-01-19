# Single-Cell Atlas of Mouse Prostate Epithelium

### Basal–Luminal Differentiation Mapped Using Seurat

---

## 📖 Project Overview

This project presents a **complete, reproducible single-cell RNA sequencing (scRNA-seq) analysis pipeline** for mouse prostate epithelium using dataset **GSE111429**. The goal is to characterize **cellular heterogeneity**, reconstruct **basal–luminal differentiation trajectories**, and generate a high-quality **Seurat-based single-cell atlas**.

The analysis follows a standard and publication-aligned Seurat workflow, with biological interpretation grounded in **canonical keratin marker expression**.

---

## 🧰 Tech Stack & Tools

![R](https://img.shields.io/badge/R-276DC3?logo=r\&logoColor=white)
![Seurat](https://img.shields.io/badge/Seurat-scRNA--seq-blue)
![SingleCell](https://img.shields.io/badge/Single--Cell-RNA--seq-green)
![Linux](https://img.shields.io/badge/Linux-000?logo=linux\&logoColor=white)
![ggplot2](https://img.shields.io/badge/ggplot2-Visualization-orange)
![HPC](https://img.shields.io/badge/HPC-Cluster-red)

---

## 1. Dataset Summary

* **Dataset:** GSE111429 (10x Genomics)
* **Organism:** *Mus musculus*
* **Tissue:** Prostate epithelium
* **Assay:** Single-cell RNA-seq

**Input Files:**

* `barcodes.tsv`
* `genes.tsv`
* `matrix.mtx`

**Biological Objective:**
Identify and characterize **basal**, **luminal**, and **intermediate epithelial states** using established keratin gene signatures.

---

## 2. Quality Control (QC)

### Filters Applied

* Cells with **≥ 1000 detected genes** retained
* **Mitochondrial content (percent.mt)** calculated and monitored
* Removal of low-quality and dying cells

### QC Visualizations

* Violin plots (nFeature_RNA, nCount_RNA, percent.mt)
* Pre- vs post-filtering QC comparison panels

**QC Outcome:**
Post-filtering data showed stabilized gene/UMI distributions and removal of mitochondrial outliers, indicating improved biological signal.

---

## 3. Normalization, HVG Selection & Scaling

Performed using Seurat’s standard workflow:

* Log-normalization (`LogNormalize`)
* Highly variable gene (HVG) detection using **vst** method
* Optional regression of mitochondrial content
* Data scaling and centering

---

## 4. Dimensionality Reduction (PCA)

* PCA performed on **top 2000 HVGs**
* PC evaluation using:

  * Elbow plot
  * PC heatmaps (PC1–PC12)

**Selected PCs:** PC1–PC12

**Interpretation:**

* Higher PCs capture epithelial lineage variation
* Several PCs enriched for keratin genes (Krt5, Krt14, Krt8, Krt18), indicating true biological structure

---

## 5. Clustering & Low-Dimensional Embedding

### Clustering

* Louvain clustering applied on PCA-reduced space
* Resolution range tested: **0.25–0.4**
* Final clustering yielded **12 distinct clusters**

### Embedding

* **UMAP** for global structure
* **t-SNE (PC1–PC12)** for publication-style visualization (consistent with Karthaus et al. 2020)

**Result:**
Clear segregation of basal and luminal epithelial populations with stable cluster boundaries.

---

## 6. Marker Gene Analysis

* Differential expression using `FindAllMarkers` (Wilcoxon test)
* Top markers extracted per cluster
* Heatmap generated for top 10 markers per cluster

### Key Lineage Markers Identified

**Basal Markers:**

* Krt5
* Krt14
* Krt15
* Krt17
* Trp63

**Luminal Markers:**

* Krt8
* Krt18
* Krt19
* Krt4
* Krt7

These genes showed high log-fold change, strong cluster specificity, and match known epithelial lineage biology.

---

## 7. Cell Type Annotation

Cluster identities were assigned using marker expression patterns and FeaturePlot overlays:

* **Basal epithelial cells:** Krt5⁺ / Krt14⁺ / Trp63⁺
* **Luminal epithelial cells:** Krt8⁺ / Krt18⁺ / Epcam⁺
* **Intermediate progenitors:** Mixed basal–luminal keratin expression

This reveals a continuous basal-to-luminal differentiation spectrum.

---

## 8. Visualizations Generated

All plots were automatically saved, including:

* QC violin plots
* HVG scatter plots
* PCA elbow and loading plots
* PC heatmaps
* UMAP and t-SNE embeddings
* Marker gene heatmaps
* FeaturePlots for keratin panels
---

## 9. Key Biological Conclusions

1. Mouse prostate epithelium shows **distinct basal and luminal compartments**
2. Basal cells are enriched for **structural keratins (Krt5/Krt14/Krt17)**
3. Luminal cells express **secretory keratins (Krt8/Krt18/Krt19)**
4. Intermediate clusters represent **transitional differentiation states**
5. Cluster architecture is consistent with published prostate epithelial biology
6. The workflow is fully reproducible and aligned with **Karthaus et al. 2020**

---

## 10. Conclusion

This analysis establishes a robust single-cell atlas of mouse prostate epithelium and demonstrates clear basal–luminal differentiation using a standard Seurat pipeline. The results highlight the power of scRNA-seq in resolving epithelial lineage structure and provide a reusable framework for future prostate and epithelial tissue studies.
