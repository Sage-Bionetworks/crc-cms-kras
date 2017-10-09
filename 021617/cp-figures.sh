#!/bin/bash

rm -rf figures/
mkdir figures/

# Figure 1
# Figure 1a: CIRC vs KRAS in TCGA
cp output/circ-tcga-all-kras-MT-vs-WT.pdf figures/
## output/tcga-all-kras-tbl.xls 
# Figure 1b: CIRC vs KRAS in KFSYSCC
cp output/circ-kfs-all-kras-MT-vs-WT.pdf figures/
## output/kfs-all-kras-tbl.xls 

# Supplemental Figure 1
# Supplemental 1a: MYC pathway vs KRAS in TCGA
cp output/myc.sig-tcga-pathways-kras-MT-vs-WT.pdf figures/
## output/tcga-pathways-kras-tbl.xls 
# Supplemental 1b: MYC pathway vs KRAS in KFSYSCC
cp output/myc.sig-kfs-pathways-kras-MT-vs-WT.pdf figures/
## output/kfs-pathways-kras-tbl.xls

# Supplemental Figure 2
# Supplemental 2a: WNT pathway vs KRAS in TCGA
cp output/wnt-tcga-pathways-kras-MT-vs-WT.pdf figures/
## output/tcga-pathways-kras-tbl.xls 
# Supplemental 2b: WNT pathway vs KRAS in KFSYSCC
cp output/wnt-kfs-pathways-kras-MT-vs-WT.pdf figures/
## output/kfs-pathways-kras-tbl.xls

# Supplemental Figure 1: 
# Supplemental Figure 1a: # of genes overlapping CIRC and Bindea in TCGA
cp output/gene-set-heatmap-tcga-all.pdf figures/
# Supplemental Figure 1b: # of genes overlapping CIRC and Bindea in KFS
cp output/gene-set-heatmap-kfs-all.pdf figures/

# Supplemental Figure 2:
# Supplemental Figure 2a: CIRC vs Bindea sets correlation in TCGA
cp output/tcga-bindea-correlation.pdf figures/
# Supplemental Figure 2a: CIRC vs Bindea sets correlation in KFS
cp output/kfs-bindea-correlation.pdf figures/

# Figure 2: CIRC vs Bindea sets correlation dendrogram in TCGA
cp output/tcga-bindea-dendrogram.pdf figures/

# Supplemental Figure 3: CIRC vs Bindea sets correlation dendrogram in KFS
cp output/kfs-bindea-dendrogram.pdf figures/

# Figure 3: Gene sets vs KRAS in TCGA and KFS
cp output/gene-set-heatmap.pdf figures/

# Figure 4: Gene sets vs KRAS in TCGA and KFS
cp output/circ-gene-heatmap.pdf figures/

## Supplemental Figure 4
## Supplemental Figure 4a: STAT1 vs KRAS in TCGA
cp output/STAT1-tcga-circ-kras-MT-vs-WT.pdf figures/
## output/tcga-circ-kras-tbl.xls
## Supplemental Figure 4b: STAT1 vs KRAS in TCGA
cp output/STAT1-kfs-circ-kras-MT-vs-WT.pdf figures/
## output/kfs-circ-kras-tbl.xls 

## Supplemental Figure 5
## Supplemental Figure 5a: CXCL10 vs KRAS in TCGA
cp output/CXCL10-tcga-circ-kras-MT-vs-WT.pdf figures/
## output/tcga-circ-kras-tbl.xls
## Supplemental Figure 5b: CXCL10 vs KRAS in KFS
cp output/CXCL10-kfs-circ-kras-MT-vs-WT.pdf figures/
## output/kfs-circ-kras-tbl.xls 

## Figure 5
## Figure 5A: CIRC vs CMS in TCGA
cp output/CIRC-tcga-all-cms.pdf figures/
cp output/scatter-tcga-all-CIRC-vs-neoantigens.pdf figures/
cp output/CIRC-tcga-all-site.pdf figures/
cp output/CIRC-tcga-all-status.pdf figures/

## Supplemental Figure 8
## Supplemental Figure 8A: CIRC vs CMS in KFS
cp output/CIRC-kfs-all-cms.pdf figures/
cp output/CIRC-kfs-all-site.pdf figures/
cp output/CIRC-kfs-all-msi-inferred-status.pdf figures/

## Supplemental Figure 9
## Supplemental Figure 9a: TCGA MSI imputation heatmap
cp output/tcga-msi.pdf figures/
## Supplemental Figure 9b: TCGA TYMS expression vs MSI/MSS
cp output/tcga-msi-tyms-vs-msi.pdf figures/

## Supplemental Figure 10
## Supplemental Figure 10a: KFS MSI imputation heatmap
cp output/kfs-msi.pdf figures/
## Supplemental Figure 10b: KFS TYMS expression vs MSI/MSS
cp output/kfs-msi-tyms-vs-annotated-msi.pdf figures/

## Figure 6
## Figure 6A: CIRC vs CMS in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-cms.pdf figures/
## cp output/CIRC-tcga-mss-msi-annotated-no-cms1-site.pdf figures/
## cp output/CIRC-kfs-mss-msi-inferred-no-cms1-site.pdf figures/
## Figure 6B: CIRC vs CMS in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-cms.pdf figures/ 

## Figure 7
## Figure 7A: multivariate analysis in TCGA
cp output/CIRC-tcga-all-full-model-forest.pdf figures/
## more output/CIRC-tcga-all-full-model-sum.tsv
## Figure 7B: multivariate analysis in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-full-model-forest.pdf figures/
## more output/CIRC-tcga-mss-msi-annotated-no-cms1-full-model-sum.tsv

## Supplemental Figure 11
## Figure 11A: multivariate analysis in KFS
cp output/CIRC-kfs-all-msi-inferred-full-model-forest.pdf figures/
## More output/CIRC-kfs-all-msi-inferred-full-model-sum.tsv
## Figure 11B: multivariate analysis in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-full-model-forest.pdf figures/
## more output/CIRC-kfs-mss-msi-inferred-no-cms1-full-model-sum.tsv

## Figure 8
## Figure 8A:  CIRC vs CMS2 MT vs others in TCGA (all)
cp output/CIRC-tcga-all-kras-cms.pdf figures/
## Figure 8B:  CIRC vs CMS2 MT vs others in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-kras-cms.pdf figures/

## Supplemental Figure 12
## Supplemental Figure 12A:  CIRC vs CMS2 MT vs others in KFS (all)
cp output/CIRC-kfs-all-kras-cms.pdf figures/
## Supplemental Figure 12B:  CIRC vs CMS2 MT vs others in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-kras-cms.pdf figures/

## Supplemental Figure 13
## Supplemental Figure 13A:  CIRC vs CMS MT vs WT in TCGA (all)
cp output/CIRC-tcga-all-kras-cms-faceted.pdf figures/
## Supplemental Figure 13B:  CIRC vs CMS MT vs WT in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-kras-cms-faceted.pdf figures/

## Supplemental Figure 14
## Supplemental Figure 14A:  CIRC vs CMS MT vs WT in KFS (all)
cp output/CIRC-kfs-all-kras-cms-faceted.pdf figures/
## Supplemental Figure 14B:  CIRC vs CMS MT vs WT in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-kras-cms-faceted.pdf figures/

## Anova
more output/CIRC-tcga-all-kras-cms-interaction-anova.tsv
more output/CIRC-tcga-mss-msi-annotated-no-cms1-kras-cms-interaction-anova.tsv 

more output/CIRC-kfs-all-kras-cms-interaction-anova.tsv
more output/CIRC-kfs-mss-msi-inferred-no-cms1-kras-cms-interaction-anova.tsv 

## Supplemental Figure 15
## Supplemental Figure 15A:  CIRC vs CMS2 MT vs CMS other MT in TCGA (all)
cp output/CIRC-tcga-all-cms-mt.pdf figures/
## Supplemental Figure 15B:  CIRC vs CMS2 MT vs CMS other MT in TCGA (MSS only)
cp output/CIRC-tcga-mss-msi-annotated-no-cms1-cms-mt.pdf figures/

## Supplemental Figure 16
## Supplemental Figure 16A:  CIRC vs CMS2 MT vs CMS other MT in KFS (all)
cp output/CIRC-kfs-all-cms-mt.pdf figures/
## Supplemental Figure 16B:  CIRC vs CMS2 MT vs CMS other MT in KFS (MSS only)
cp output/CIRC-kfs-mss-msi-inferred-no-cms1-cms-mt.pdf figures/
