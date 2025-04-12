# scDCABC

**A Deep Biclustering Method Integrating Denoising, Dimensionality Reduction, and Clustering for scRNA-seq Data**

## 📘 Introduction

**scDCABC** is a deep biclustering algorithm designed for single-cell RNA sequencing (scRNA-seq) data. It integrates:
- **Deep Count Autoencoder (DCA)** for denoising and dimensionality reduction,
- **Walktrap algorithm** for cell clustering on a shared nearest neighbor (SNN) graph,
- and **gene clustering** to discover gene modules associated with specific cell subpopulations.

This method is robust to the zero-inflation problem in scRNA-seq data and provides enhanced bio-interpretability.

---

## 📄 Citation

If you use this code, please cite:

**Tang, Xiaoqi** and **Lan, Chaowang** (2024):  
*scDCABC: A deep biclustering method integrating denoising, dimensionality reduction, and clustering for scRNA-seq data*.  
Proceedings of the 2024 International Conference on Bioinformatics and Biomedical Science (ICBBS 2024).  
DOI: [10.1145/3704198.3704204](https://doi.org/10.1145/3704198.3704204)

---


## 🛠 Requirements

This implementation is tested with:

- Python 3.8+
- TensorFlow 2.4.4
- Keras 2.4.3
- NumPy 1.19.5
- Pandas 1.4.4
- SciPy 1.10.1
- scikit-learn 1.3.2
- igraph 0.11.6
- scanpy 1.8.2
- anndata 0.7.8
