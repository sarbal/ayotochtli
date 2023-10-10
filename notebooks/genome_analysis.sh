


```
cat Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.gtf |  awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,".",$7}}' | tr -d '";'  > Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.bed

bedtools flank -i Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.bed -g chrNameLength.txt -l 2000 -r 0 -s > genes.2kb.promoters.bed

bedtools getfasta -fi dasNov3_spike.fa -bed genes.2kb.promoters.bed -fo genes.2kb.promoters.bed.fa

bedtools getfasta -fi dasNov3_spike.fa -bed Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.bed -fo genes.fa

cat Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.gtf |  awk 'OFS="\t" {if ($3=="CDS") {print $1,$4-1,$5,$10,".",$7}}' | tr -d '";'  > Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.CDS.bed

cat Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.gtf |  awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,".",$7}}' | tr -d '";'  > Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.exons.bed

bedtools nuc -fi dasNov3_spike.fa -bed Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.exons.bed > Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.genome.content

bedtools nuc -fi dasNov3_spike.fa -bed Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.CDS.bed > Dasypus_novemcinctus.Dasnov3.0.95_mod.spike.genome.content.CDS
```
