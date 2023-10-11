### Data

#### Expression data 
Expression values (gene counts and CPM) for all quads at all timepoints are in: 
- [`counts_strand_comb.Rdata`](counts_strand_comb.Rdata)
- [`cpm_strand_comb.Rdata`](cpm_strand_comb.Rdata)

Allele specific expression data can be found in: 
- [`exprs_all.Rdata`](exprs_all.Rdata)

#### Phenotypes and metadata 
- [`armadillo_hem.Rdata`](armadillo_hem.Rdata)
- [`metadata.Rdata`]()

#### Genotype information
- [`snps_overlaps.Rdata`](snps_overlaps.Rdata)


### Analysis outputs
#### Testing/training splits for identity analysis 
- [`ase_ratios.test_train.Rdata`](ase_ratios.test_train.Rdata)
- [`cpm_test_train.Rdata`](cpm_test_train.Rdata)
- 
#### XCI skewing estimates 
- [`skew.est.max.genes.Rdata`](skew.est.max.genes.Rdata)
- [`estimated_ratios.Rdata`](estimated_ratios.Rdata)

#### Models for identity analysis
- [`ase_rand.null.ratios2.xsub.Rdata`](ase_rand.null.ratios2.xsub.Rdata)
- [`ase_rand.ratios.noX.Rdata`](ase_rand.ratios.noX.Rdata)
- [`final.model.Rdata`](final.model.Rdata)


#### Perfect predictor genes 
Genes that were able to predict armadillos across time, per quad can be found in this file: 
- [`perfect.pred.arm.Rdata`](perfect.pred.arm.Rdata)

#### SNP binomial modeling analyses
- [`black_list_homo_snps.rds`](black_list_homo_snps.rds)
- [`ase_cov.rds`](ase_cov.rds)
The remaining files are too large to be uploaded here, but can be located at: 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8429576.svg)](https://doi.org/10.5281/zenodo.8429576)

 
#### Other data 
Genome annotation
- [`gene_annotations_v0_95_mod.Rdata`](gene_annotations_v0_95_mod.Rdata)
Functional gene sets:
- [`GO.arm.Rdata`](GO.arm.Rdata)
Homologs: 
- [`homologs.Rdata`](homologs.Rdata)
X chromosome chains and objects
- [`chrX.chainHg38.txt.gz`](chrX.chainHg38.txt.gz)
- [`chrX.chainMm10.txt.gz`](chrX.chainMm10.txt.gz)
- [`chromInfo.txt.gz`](chromInfo.txt.gz)
- [`cytoBandIdeo.txt.gz`](cytoBandIdeo.txt.gz)
  
### Code
See [notebooks](../notebooks/) and the [bin](../notebooks/bin/) for all the analyses and scripts. 
 
