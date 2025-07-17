# CADD Score Optimization for ClinVar Variants

## Overview
This notebook implements a comprehensive analysis pipeline for optimizing CADD score thresholds across different genomic regions using ClinVar variant classifications. The analysis includes:

1. **Data preprocessing and standardization**
2. **Statistical analysis of CADD scores by genomic region**
3. **Threshold optimization using ROC and PR curves**
4. **Visualization and validation of results**

## Key Features
- ANNOVAR functional annotation standardization
- Enhanced Youden's J statistic normalization
- Category-specific threshold optimization
- Comprehensive statistical testing
- Production-ready filtering functions

## Reusable component

- CADD_ClinVar_threshold_optimization_data.parquet serves as the input dataset for this analysis. To reproduce the same results, refer to section 4.1 of the Entrypoint using the corresponding Parquet file available on GitHub.

- This script is designed to be generalizable and can be reapplied to other datasets, provided they conform to the same data structure.