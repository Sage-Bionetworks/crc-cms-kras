#!/bin/bash

rm -rf figures/
mkdir figures/

# Figure 1
# Figure 1a: CIRC vs KRAS in TCGA
cp output/circ-tcga-all-kras-MT-vs-WT.pdf figures/fig1a.pdf
## output/tcga-all-kras-tbl.xls 
# Figure 1b: CIRC vs KRAS in KFSYSCC
cp output/circ-kfs-all-kras-MT-vs-WT.pdf figures/fig1b.pdf
## output/kfs-all-kras-tbl.xls 

## more run.out | grep 'CIRC vs KRAS' | more | grep '\-all:'
## tcga-all: Analyzing CIRC vs KRAS using Wilcoxon rank sum test with continuity correction: W = 16382, n_WT = 219, n_MT = 125, p = 0.00239184964757301
## kfs-all: Analyzing CIRC vs KRAS using Wilcoxon rank sum test with continuity correction: W = 12084, n_WT = 175, n_MT = 115, p = 0.00381629339682518

# Supplemental Figure 1: 
# Supplemental Figure 1a: # of genes overlapping CIRC and Bindea in TCGA
cp output/gene-set-heatmap-tcga-all.pdf figures/supp-fig1a.pdf
# Supplemental Figure 1b: # of genes overlapping CIRC and Bindea in KFS
cp output/gene-set-heatmap-kfs-all.pdf figures/supp-fig1b.pdf

# Supplemental Figure 2:
# Supplemental Figure 2a: CIRC vs Bindea sets correlation in TCGA
cp output/tcga-bindea-correlation.pdf figures/supp-fig2a.pdf
# Supplemental Figure 2b: CIRC vs Bindea sets correlation in KFS
cp output/kfs-bindea-correlation.pdf figures/supp-fig2b.pdf

# Figure 2A: CIRC vs Bindea sets correlation dendrogram in TCGA
cp output/tcga-bindea-dendrogram.pdf figures/fig2a.pdf

# Supplemental Figure 3: CIRC vs Bindea sets correlation dendrogram in KFS
cp output/kfs-bindea-dendrogram.pdf figures/supp-fig3.pdf

# Figure 2B: Gene sets vs KRAS in TCGA and KFS
cp output/gene-set-heatmap.pdf figures/fig2b.pdf
more run.out | grep -i tcga-gene | grep KRAS | grep -i cytotoxic
## tcga-gene-sets: Analyzing Cytotoxic cells vs KRAS using Wilcoxon rank sum test with continuity correction: W = 15515, n_WT = 219, n_MT = 125, p = 0.0394550958525481
more run.out | grep -i kfs-gene | grep KRAS | grep -i cytotoxic
## kfs-gene-sets: Analyzing Cytotoxic cells vs KRAS using Wilcoxon rank sum test with continuity correction: W = 11662, n_WT = 175, n_MT = 115, p = 0.0220858283559412

more run.out | grep -i tcga-gene | grep KRAS | grep -i neutro   
## tcga-gene-sets: Analyzing Neutrophils vs KRAS using Wilcoxon rank sum test with continuity correction: W = 15981, n_WT = 219, n_MT = 125, p = 0.00974672529331536
more run.out | grep -i kfs-gene | grep KRAS | grep -i neutro   
## kfs-gene-sets: Analyzing Neutrophils vs KRAS using Wilcoxon rank sum test with continuity correction: W = 11974, n_WT = 175, n_MT = 115, p = 0.00622854051911139

# Supplemental Figure 4: CIRC vs KRAS in TCGA and KFS
cp output/circ-gene-heatmap.pdf figures/supp-fig4.pdf

## Supplemental Figure 5
## Supplemental Figure 5a: STAT1 vs KRAS in TCGA
cp output/STAT1-tcga-circ-kras-MT-vs-WT.pdf figures/supp-fig5a.pdf
## output/tcga-circ-kras-tbl.xls
## Supplemental Figure 5b: STAT1 vs KRAS in TCGA
cp output/STAT1-kfs-circ-kras-MT-vs-WT.pdf figures/supp-fig5b.pdf
## output/kfs-circ-kras-tbl.xls 

more run.out | grep -i KRAS | grep '\-circ' | grep kfs | grep STAT1
## kfs-circ: Analyzing STAT1 vs KRAS using Wilcoxon rank sum test with continuity correction: W = 13029, n_WT = 175, n_MT = 115, p = 2.17930782370407e-05

more run.out | grep -i KRAS | grep '\-circ' | grep tcga  | grep STAT1
## tcga-circ: Analyzing STAT1 vs KRAS using Wilcoxon rank sum test with continuity correction: W = 16069, n_WT = 219, n_MT = 125, p = 0.0072772800538876

## Supplemental Figure 6
## Supplemental Figure 6a: CXCL10 vs KRAS in TCGA
cp output/CXCL10-tcga-circ-kras-MT-vs-WT.pdf figures/supp-fig6a.pdf
## output/tcga-circ-kras-tbl.xls
## Supplemental Figure 6b: CXCL10 vs KRAS in KFS
cp output/CXCL10-kfs-circ-kras-MT-vs-WT.pdf figures/supp-fig6b.pdf
## output/kfs-circ-kras-tbl.xls 

more run.out | grep -i KRAS | grep '\-circ' | grep kfs | grep CXCL10
## kfs-circ: Analyzing CXCL10 vs KRAS using Wilcoxon rank sum test with continuity correction: W = 12348, n_WT = 175, n_MT = 115, p = 0.00107219740701228

more run.out | grep -i KRAS | grep '\-circ' | grep tcga  | grep CXCL10
## tcga-circ: Analyzing CXCL10 vs KRAS using Wilcoxon rank sum test with continuity correction: W = 16218, n_WT = 219, n_MT = 125, p = 0.0043467533177559

more run.out | grep -i KRAS | grep '\-circ' | grep tcga  | grep CIITA 
## tcga-circ: Analyzing CIITA vs KRAS using Wilcoxon rank sum test with continuity correction: W = 15715, n_WT = 219, n_MT = 125, p = 0.0223216042395995

more run.out | grep -i KRAS | grep '\-circ' | grep kfs  | grep CIITA
## kfs-circ: Analyzing CIITA vs KRAS using Wilcoxon rank sum test with continuity correction: W = 10654, n_WT = 175, n_MT = 115, p = 0.397560023143556

more run.out | grep -i KRAS | grep '\-circ' | grep kfs  | grep CD247
## kfs-circ: Analyzing CD247 vs KRAS using Wilcoxon rank sum test with continuity correction: W = 11319, n_WT = 175, n_MT = 115, p = 0.0721926652054417

more run.out | grep -i KRAS | grep '\-circ' | grep tcga  | grep CD247
## tcga-circ: Analyzing CD247 vs KRAS using Wilcoxon rank sum test with continuity correction: W = 15908, n_WT = 219, n_MT = 125, p = 0.0123353892395922


## Figure 3
## Figure 3A: CIRC vs CMS in TCGA
cp output/CIRC-tcga-all-cms.pdf figures/fig3a.pdf
## Figure 3B: CIRC vs neoantigen in TCGA
cp output/scatter-tcga-all-CIRC-vs-neoantigens.pdf figures/fig3b.pdf
## Figure 3C: CIRC vs tumor site in TCGA
cp output/CIRC-tcga-all-site.pdf figures/fig3c.pdf
## Figure 3D: CIRC vs msi status in TCGA
cp output/CIRC-tcga-all-status.pdf figures/fig3d.pdf

more run.out | grep -i '\-all:' | grep -i circ | grep -i cms | grep 'CMS2 vs' | grep -v 'KRAS MT only' | grep -E 'CMS1|CMS4' | grep -i tcga
## tcga-all: Analyzing CIRC CMS2 vs CMS1 using Wilcoxon rank sum test with continuity correction: W = 6373, n_CMS2 = 136, n_CMS1 = 51, p = 1.23959499498166e-18
## tcga-all: Analyzing CIRC CMS2 vs CMS4 using Wilcoxon rank sum test with continuity correction: W = 2243, n_CMS2 = 136, n_CMS4 = 87, p = 5.51962791445972e-15

more run.out | grep 'CIRC vs neoantigens' | grep 'all CIRC'
## tcga-all CIRC vs neoantigens: t = 7.18856145436229 df = 340 nrow = 342 cor = 0.363227593412596 pval = 4.1846192883292e-12

more run.out | grep -i '\-all:' | grep -i circ | grep -i site | grep -i tcga
## tcga-all: Analyzing CIRC site: rectum vs leftWilcoxon rank sum test with continuity correction: W = 3359, n_rectum = 60, n_left = 133, p = 0.0791927710817998
## tcga-all: Analyzing CIRC site: right vs leftWilcoxon rank sum test with continuity correction: W = 6860, n_right = 148, n_left = 133, p = 1.16662152475699e-05

more run.out | grep -i '\-all:' | grep -i circ | grep MSI   
## tcga-all: Analyzing CIRC status: MSS vs MSIWilcoxon rank sum test with continuity correction: W = 12383, n_MSS = 289, n_MSI = 53, p = 9.36844862356913e-13

## For 74 of 76 CMS1 samples imputed to be MSI, search run.out for Overlap
## Overlap of MSI and CMS1 in TCGA inferred:
##    
##      CMS1 CMS2 CMS3 CMS4 NOLBL
##  msi   74    1   42   33    19
##  mss    2  214   27  105    42

## For 78 of 81 MSI samples correctly identified in TCGA, search run.out for:
## Overlap of MSI inferred and annotated in TCGA:
##     
##      msi mss
##  msi  78   3
##  mss  91 387


## Supplemental Figure 7
## Supplemental Figure 7A: CIRC vs CMS in KFS
cp output/CIRC-kfs-all-cms.pdf figures/supp-fig7a.pdf
## Figure 7B: CIRC vs tumor site in KFS
cp output/CIRC-kfs-all-site.pdf figures/supp-fig7b.pdf
## Figure 7C: CIRC vs msi status in KFS
cp output/CIRC-kfs-all-msi-inferred-status.pdf figures/supp-fig7c.pdf

more run.out | grep -i '\-all:' | grep -i circ | grep -i cms | grep 'CMS2 vs' | grep -v 'KRAS MT only' | grep -E 'CMS1|CMS4' | grep -i kfs
## kfs-all: Analyzing CIRC CMS2 vs CMS1 using Wilcoxon rank sum test with continuity correction: W = 2070, n_CMS2 = 106, n_CMS1 = 26, p = 7.60378857272565e-05
## kfs-all: Analyzing CIRC CMS2 vs CMS4 using Wilcoxon rank sum test with continuity correction: W = 2166, n_CMS2 = 106, n_CMS4 = 78, p = 3.57087557384974e-08

more run.out | grep -i '\-all:' | grep -i circ | grep -i site | grep -i kfs 
## kfs-all: Analyzing CIRC site: rectum vs leftWilcoxon rank sum test with continuity correction: W = 2624, n_rectum = 188, n_left = 28, p = 0.980606268877483
## kfs-all: Analyzing CIRC site: right vs leftWilcoxon rank sum test with continuity correction: W = 886, n_right = 74, n_left = 28, p = 0.262273892076771

more run.out | grep -i kfs | grep -i infer | grep -i circ | grep MSI
## kfs-all-msi-inferred: Analyzing CIRC status: MSS vs MSIWilcoxon rank sum test with continuity correction: W = 11063, n_MSS = 212, n_MSI = 78, p = 1.019523849464e-05


## Supplemental Figure 8
cp output/neoantigens-tcga-all-cms.pdf figures/supp-fig8.pdf

## Supplemental Figure 9
## Supplemental Figure 9a: TCGA MSI imputation heatmap
cp output/tcga-msi.pdf figures/supp-fig9a.pdf
## Supplemental Figure 9b: TCGA TYMS expression vs MSI/MSS
cp output/tcga-msi-tyms-vs-msi.pdf figures/supp-fig9b.pdf

more run.out | grep -i tyms | grep Status | grep tcga-msi | grep Imputed
## tcga-msi: Analyzing TYMS vs Imputed Status using Wilcoxon rank sum test with continuity correction: W = 54712 p = 2.44358357254431e-35 n_MSS = 390, n_MSI = 169

## Supplemental Figure 10
## Supplemental Figure 10a: KFS MSI imputation heatmap
cp output/kfs-msi.pdf figures/supp-fig10a.pdf
## Supplemental Figure 10b: KFS TYMS expression vs MSI/MSS
cp output/kfs-msi-tyms-vs-annotated-msi.pdf figures/supp-fig10b.pdf

## For 29 of 30 CMS1 samples imputed to be MSI in KFS, search run.out for Overlap

## Overlap of MSI and CMS1 in KFS inferred:
##    
##      CMS1 CMS2 CMS3 CMS4 NOLBL
##  msi   29    1   13   38     7
##  mss    1  109   36   47    26

more run.out | grep -i tyms | grep Status | grep -v tcga
## Beeswarm: TYMS vs Inferred.Status W = 11946 p = 0.00102436427595825 n_MSI = 88, n_MSS = 219

## Supplemental Figure 11
## Supplemental Figure 11A: CIRC vs CMS in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-cms.pdf figures/supp-fig11a.pdf
## Supplemental Figure 11B: CIRC vs CMS in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-cms.pdf figures/supp-fig11b.pdf 

more run.out | grep -i '\-mss-msi-annotated-no-cms1' | grep -i tcga | grep -i circ | grep 'CMS2 vs CMS4' | grep -v KRAS
## tcga-mss-msi-annotated-no-cms1: Analyzing CIRC CMS2 vs CMS4 using Wilcoxon rank sum test with continuity correction: W = 2117, n_CMS2 = 134, n_CMS4 = 81, p = 7.02496068606718e-14

more run.out | grep -i '\-mss-msi-inferred-no-cms1' | grep -i kfs | grep -i circ | grep 'CMS2 vs CMS4' | grep -v KRAS
## kfs-mss-msi-inferred-no-cms1: Analyzing CIRC CMS2 vs CMS4 using Wilcoxon rank sum test with continuity correction: W = 1576, n_CMS2 = 105, n_CMS4 = 46, p = 0.000698953453612531


## Figure 4
## Figure 4A: multivariate analysis in TCGA
cp output/CIRC-tcga-all-full-model-forest.pdf figures/fig4a.pdf
## more output/CIRC-tcga-all-full-model-sum.tsv
## Figure 4B: multivariate analysis in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-full-model-forest.pdf figures/fig4b.pdf
## more output/CIRC-tcga-mss-msi-annotated-no-cms1-full-model-sum.tsv

## Supplemental Figure 12
## Figure 12A: multivariate analysis in KFS
cp output/CIRC-kfs-all-msi-inferred-full-model-forest.pdf figures/supp-fig12a.pdf
## More output/CIRC-kfs-all-msi-inferred-full-model-sum.tsv
## Figure 12B: multivariate analysis in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-full-model-forest.pdf figures/supp-fig12b.pdf
## more output/CIRC-kfs-mss-msi-inferred-no-cms1-full-model-sum.tsv

## Figure 5
## Figure 5A:  CIRC vs CMS2 MT vs others in TCGA (all)
cp output/CIRC-tcga-all-kras-cms.pdf figures/fig5a.pdf
## Figure 5B:  CIRC vs CMS2 MT vs others in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-kras-cms.pdf figures/fig5b.pdf

## Supplemental Figure 13
## Supplemental Figure 13A:  CIRC vs CMS2 MT vs others in KFS (all)
cp output/CIRC-kfs-all-kras-cms.pdf figures/supp-fig13a.pdf
## Supplemental Figure 13B:  CIRC vs CMS2 MT vs others in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-kras-cms.pdf figures/supp-fig13b.pdf

## Supplemental Figure 14
## Supplemental Figure 14A:  CIRC vs CMS MT vs WT in TCGA (all)
cp output/CIRC-tcga-all-kras-cms-faceted.pdf figures/supp-fig14a.pdf
## Supplemental Figure 14B:  CIRC vs CMS MT vs WT in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-kras-cms-faceted.pdf figures/supp-fig14b.pdf

## Supplemental Figure 15
## Supplemental Figure 15A:  CIRC vs CMS MT vs WT in KFS (all)
cp output/CIRC-kfs-all-kras-cms-faceted.pdf figures/supp-fig15a.pdf
## Supplemental Figure 15B:  CIRC vs CMS MT vs WT in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-kras-cms-faceted.pdf figures/supp-fig15b.pdf


