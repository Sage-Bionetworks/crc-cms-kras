more tcga-immune-hallmark-kras.tsv | grep -E 'pval|Cyto'

more kfs-immune-hallmark-kras.tsv | grep -E 'pval|Cyto'

more tcga-immune-hallmark-kras.tsv | grep -E 'pval|Neutro'

more kfs-immune-hallmark-kras.tsv | grep -E 'pval|Neutro'

more tcga-immune-hallmark-kras.tsv | grep -E 'pval|Th1'

more kfs-immune-hallmark-kras.tsv | grep -E 'pval|Th1'

more kfs-immune-hallmark-cms.tsv  | grep -E 'gene.set|CIRC' | cut -f1,2,4,6,8

more tcga-immune-hallmark-cms.tsv  | grep -E 'gene.set|CIRC' | cut -f1,2,4,6,8

more output/CIRC-kfs-kras-cms-interaction-anova.tsv 

more output/CIRC-tcga-kras-cms-interaction-anova.tsv 

more kfs-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'pval|CIRC|cytot|neutro|Th1' | cut -f1,2,3,5,6

more tcga-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'pval|CIRC|cytot|neutro|Th1' | cut -f1,2,3,4,5


more kfs-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'pval|gamma' 

more tcga-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'pval|gamma'

more kfs-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'myc|wnt|pval' | cut -f1,4

more tcga-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'myc|wnt|pval' | cut -f1,6

##more output/STAT1-kfs-full-model-sum.tsv
##more output/CXCL10-kfs-full-model-sum.tsv
##more output/CIITA-kfs-full-model-sum.tsv
##more output/STAT1-tcga-full-model-sum.tsv
##more output/CXCL10-tcga-full-model-sum.tsv
##more output/CIITA-tcga-full-model-sum.tsv
##more output/HALLMARK_INTERFERON_GAMMA_RESPONSE-tcga-full-model-sum.tsv
##more output/HALLMARK_INTERFERON_GAMMA_RESPONSE-kfs-full-model-sum.tsv

##grep -E 'KRASWT|CMS1|CMS4' output/*kfs-full-model-sum.tsv | grep -v INTER

grep -E 'KRASWT' output/*-full-model-sum.tsv | grep -E 'STAT1|CXCL10|CIITA'
grep -E 'CMS1' output/*-full-model-sum.tsv | grep -E 'STAT1|CXCL10|CIITA'
grep -E 'CMS4' output/*-full-model-sum.tsv | grep -E 'STAT1|CXCL10|CIITA'

more kfs-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'pval|stat1|cxcl10|ciita' | cut -f1,4
more tcga-immune-hallmark-kras-cms.tsv  | cut -f1,2,6,10,14,18,22,26 | grep -E -i 'pval|stat1|cxcl10|ciita' | cut -f1,6


