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
- [``]()

#### Perfect predictor genes 
Genes that were able to predict armadillos across time, per quad can be found in this file: 
- [`perfect.pred.arm.Rdata`](perfect.pred.arm.Rdata)

#### SNP binomial modeling analyses
- [`black_list_homo_snps.rds`](black_list_homo_snps.rds)
- [`ase_cov.rds`](ase_cov.rds)
- [`cumsum___emp_p_adj_cor_All.tsv`](cumsum___emp_p_adj_cor_All.tsv)
- [`cumsum___emp_p_adj_cor_No X.tsv`](cumsum___emp_p_adj_cor_No X.tsv)
- [`cumsum___emp_p_adj_cor_X.tsv`](cumsum___emp_p_adj_cor_X.tsv)
- [`qrank_data___all_emp_p_adj_cor.rds`](qrank_data___all_emp_p_adj_cor.rds)

 
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
 
