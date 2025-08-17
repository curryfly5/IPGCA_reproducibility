# IPGCA_reproducibility
This repository contains the code that was used for the Integrated Pig Cell Atlas (IPGCA) project.
## Reproduce:
You can find the main code used in the IPGCA project in the notebooks of this repository. We subdivided the notebooks into nine main categories:
| Folder                  | Content                                                                    |
| ----------------------- | -------------------------------------------------------------------------- |
| `00_preprocess`         | Download public datasets, quality filtering, basic normalisation           |
| `01_batch_effect`       | Batch-correction benchmarking (Harmony, Scanorama, scVI...)                |
| `02_annotation`         | Cell-type annotation & marker gene selection                               |
| `03_cov_variance`       | Covariate-variance decomposition of technical and biological factors       |
| `04_biological_finding` | Core biological results (e.g. rare cell discovery, trajectory inference)   |
| `05_deconvolution`      | Generate cell-type signature matrices & deconvolve bulk RNA-seq            |
| `06_IPGCA_extend`       | Extension analyses: mapping new datasets into IPGCA                        |
| `07_colocalization`     | Co-localisation of eQTLs with GWAS hits using IPGCA                        |
| `08_cross_species`      | Cross-species alignment (human–mouse–pig) using IPGCA                      |

and wellcome to single cell gut database based on IPGCA (https://alphaindex.zju.edu.cn/scgut)
