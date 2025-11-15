# Single-Cell Atlas of Mouse Prostate Epithelium: Basal–Luminal Differentiation Mapped Using Seurat

## 📖 Overview
This project presents a complete, reproducible **single-cell RNA sequencing (scRNA-seq)** analysis pipeline using the **mouse prostate epithelial dataset (GSE111429)**.  
The primary objective is to explore **cellular heterogeneity**, reconstruct **basal vs. luminal epithelial lineages**, and generate a clean **Seurat-based atlas** of the tissue.

The workflow includes:
- Raw data loading (10x Genomics matrix)
- Quality control & filtering
- Normalization and HVG detection
- PCA and dimensionality reduction
- Clustering (Louvain)
- t-SNE & UMAP visualization
- Marker gene analysis
- Basal/Luminal identity assignment using canonical keratin genes

This README summarizes the analysis, results, and biological interpretation.

---

# 1. Dataset Summary

**Dataset:** GSE111429 (10x Genomics)  
**Organism:** *Mus musculus*  
**Tissue:** Prostate epithelium  
**Assay:** Single-cell RNA-seq  
**Files Used:**  
- `barcodes.tsv`  
- `genes.tsv`  
- `matrix.mtx`

**Biological Goal:**  
Identify and characterize **basal**, **luminal**, and **intermediate** cell states by leveraging keratin gene expression patterns.

---

# 2. Quality Control (QC)

### **Filters Applied**
- **≥ 1000 genes** per cell  
- Mitochondrial genes: **percent.mt recorded and monitored**  
- Removal of low-quality / dead cells  
- Before vs. After QC visualized using:
  - Violin plots  
  - Combined QC panel across nGenes, nUMI, percent.mt  

### **QC Insights**
- Post-filtering: Data became cleaner and biologically consistent  
- nUMI and nGene distributions stabilized  
- percent.mt outliers removed

---

# 3. Normalization, HVG Selection & Scaling

Performed using Seurat’s pipeline:
- Log-normalization (`LogNormalize`)
- Highly variable genes detected using **vst**
- Optionally regressed out mitochondrial content
- Data scaled and centered

---

# 4. Dimensionality Reduction (PCA)

- PCA performed on the top 2000 HVGs
- PC inspection done via:
  - PC heatmaps (PC1–PC12)
  - Elbow plot  
- PC1–PC12 selected for downstream clustering (based on variance knee point)

### **Interpretation**
- Higher PCs separate epithelial subtypes  
- Some PCs dominated by keratin signatures (Krt5, Krt14, Krt8, Krt18) indicating real lineage differences

---

# 5. Clustering & Embedding

## **Clustering**
- Louvain algorithm applied on PCA-reduced data  
- Selected resolution: **0.25–0.4**  
- Yielded **12 clusters**

## **Embedding**
Generated both:
- **UMAP**
- **t-SNE (PC1–PC12)** — similar to Karthaus 2020 paper

### **Result**
- Well-separated structure of epithelial states
- Compact and stable clusters
- Clear basal vs. luminal segregation

---

# 6. Marker Gene Analysis

Used:
- `FindAllMarkers` (Wilcoxon test)
- Top markers extracted per cluster
- Heatmap generated for top 10 markers per cluster

**Key markers identified (from your results):**

### **Basal Markers**
- **Krt5**
- **Krt14**
- **Krt15**
- **Krt17**
- **Trp63**

### **Luminal Markers**
- **Krt8**
- **Krt18**
- **Krt19**
- **Krt4**
- **Krt7**

### **Why these genes?**
These genes appeared at the **top of the FindAllMarkers output**, meaning:
- Strong log2 fold change  
- High pct.1 (cluster-specific expression)  
- Known epithelial lineage markers  

Therefore, they accurately label your clusters.

---

# 7. Cell Type Annotation

Based on marker expression:

### **Basal clusters (Krt5/Krt14/Trp63+)**
- High in structural keratins  
- Stem/progenitor features  
- Match classical basal epithelial lineage

### **Luminal clusters (Krt8/Krt18/Epcam+)**
- Secretory cell signatures  
- Represent differentiated luminal compartment  

### **Intermediate / Progenitor clusters**
- Mixed Krt5/Krt8 expression  
- Transitional states between basal ↔ luminal

Cluster IDs were manually annotated using FeaturePlot overlays.

---

# 8. Visualizations Included

All plots were saved automatically:
- QC violin plots  
- PCA loading plots  
- HVG scatter  
- Elbow plot  
- PC heatmaps  
- UMAP  
- t-SNE  
- Marker heatmap  
- FeaturePlot for keratin panels  
- Cluster-wise violin plots  

File location:  
`C:/Users/mmsid/Documents/My Project/scRNA Analysis/GSE111429_RAW/`

---

# 9. Key Biological Conclusions

1. **Mouse prostate epithelium contains clear basal and luminal compartments.**  
2. **Basal cells show overexpression of structural keratins Krt5/Krt14/Krt17.**  
3. **Luminal cells express secretory keratins Krt8/Krt18/Krt19.**  
4. **Intermediate populations indicate differentiation transitions.**  
5. **Cluster architecture matches known prostate epithelial biology.**  
6. Workflow is fully reproducible and aligned with Karthaus et al. 2020 pipeline style.

---


