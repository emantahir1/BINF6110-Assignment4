# BINF*6110 Assignment 4: scRNA-seq Analysis of IAV Nasal Mucosa

**Course:** BINF*6110 - Applied Bioinformatics   
**Dataset:** Kazer et al. 2024 - Mouse nasal mucosa scRNA-seq following IAV infection     
**BioProject:** [PMC11324402](https://pmc.ncbi.nlm.nih.gov/articles/PMC11324402/)

---

## Table of Contents
1. [Introduction](#introduction)
2. [Methods](#methods)
3. [Results](#results)
4. [Discussion](#discussion)
5. [References](#references)

---

## Introduction

Influenza A virus (IAV) is a negative-sense, single-stranded RNA virus and a major respiratory pathogen affecting both humans and animals worldwide (Krammer et al., 2018). It causes an estimated three to five million cases of severe illness each year and remains a significant public health concern due to its ability to rapidly evolve and cause pandemics (WHO, 2023). Infection begins in the nasal mucosa, which serves as the primary site of viral entry and early replication (Oslund & Bhatt, 2013). The nasal cavity contains several distinct regions, including the respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LNG), each with different cell types and roles in immune defense and viral detection (Kazer et al., 2024).

The innate immune response to IAV is triggered within hours of infection and is largely driven by immune cells in the nasal mucosa (Iwasaki & Pillai, 2014). Among these, macrophages play a key role as both early detectors of infection and regulators of the inflammatory response. These cells recognize viral components through pattern recognition receptors, which activate signaling pathways that lead to the production of type I and type II interferons, as well as interferon-stimulated genes (ISGs) (Schneider et al., 2014). This interferon response is one of the earliest and most important antiviral defenses, and both its timing and intensity can strongly influence disease outcome (Iwasaki & Pillai, 2014).

Studying how different cell types respond to IAV infection requires methods that can capture variation at the single-cell level. Bulk RNA sequencing averages gene expression across all cells in a sample, which can hide important differences between cell types and overlook rare populations (Jovic et al., 2022). In contrast, single-cell RNA sequencing (scRNA-seq) profiles gene expression in individual cells, allowing for the identification of distinct cell types, their marker genes, and how each population responds to infection (Jovic et al., 2022). This is especially useful in complex tissues like the nasal mucosa, where many different cell types coexist and respond differently to viral infection.

For analyzing scRNA-seq data, Seurat is a widely used framework that provides tools for quality control, normalization, dimensionality reduction, clustering, and visualization (Hao et al., 2021). One common challenge in these datasets is batch effects caused by technical differences between samples, which can interfere with biological interpretation (Luecken et al., 2022). Harmony is often used to address this issue by correcting batch effects in PCA space, allowing cells from different samples to be more accurately compared (Korsunsky et al., 2019). After integration and clustering, differential expression analysis is best performed using pseudobulk approaches, which group cells by sample and apply bulk RNA-seq methods such as DESeq2. This reduces false positives compared to analyzing individual cells directly (Squair et al., 2021). To interpret these results, gene set enrichment analysis (GSEA) can be used to identify biological pathways that are enriched among differentially expressed genes (Subramanian et al., 2005).

In this study, we reanalyze the dataset from Kazer et al. (2024) to examine how nasal mucosa cell populations respond to IAV infection in mice. We focus specifically on the respiratory mucosa, comparing uninfected mice to those five days post-infection, a time point associated with peak myeloid cell recruitment. Using Seurat, Harmony, and pseudobulk DESeq2, we identify 29 distinct cell populations, annotate their identities, and investigate how macrophages respond transcriptionally to infection.
