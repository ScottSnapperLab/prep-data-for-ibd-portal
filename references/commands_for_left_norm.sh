bcftools annotate -x INFO/MLEAC,INFO/MLEAF \
-Ou /home/gus/MEGAsync/zim/main/BCH/Projects/IBD_portal/data_preprocessing/test_fams/VCFs/CHB099.vcf.bgz \
| \
bcftools annotate -a /home/gus/MEGAsync/zim/main/BCH/Projects/IBD_portal/data_preprocessing/pipeline_output/test_fams/make_anno_tables/CHB099.csv.gz \
-h /home/gus/MEGAsync/zim/main/BCH/Projects/IBD_portal/data_preprocessing/pipeline_output/test_fams/make_anno_tables/CHB099.hdr \
-c CHROM,POS,AGE_MO,BAM_OK,FAM_ID,MOI \
-Ov \
| \
sed -e 's/AD,Number=./AD,Number=R/g' -e 's/PL,Number=./PL,Number=G/g' | \
bcftools norm \
-f /run/media/gus/Storage/BCH/data/g1k/reference_genome/human_g1k_v37.fasta \
-m- \
-Ov - \
| \
vt normalize - \
-r /run/media/gus/Storage/BCH/data/g1k/reference_genome/human_g1k_v37.fasta \
-o  CHB099.vt.norm.none.vcf



bcftools annotate -a /home/gus/MEGAsync/zim/main/BCH/Projects/IBD_portal/data_preprocessing/pipeline_output/test_fams/make_anno_tables/CHB099.csv.gz \
-h /home/gus/MEGAsync/zim/main/BCH/Projects/IBD_portal/data_preprocessing/pipeline_output/test_fams/make_anno_tables/CHB099.hdr \
-c CHROM,POS,AGE_MO,BAM_OK,FAM_ID,MOI \
-Ov \
/home/gus/MEGAsync/zim/main/BCH/Projects/IBD_portal/data_preprocessing/test_fams/VCFs/CHB099.vcf.bgz \
| \
sed -e 's/AD,Number=./AD,Number=R/g' -e 's/PL,Number=./PL,Number=G/g' | \
bcftools norm \
-f /run/media/gus/Storage/BCH/data/g1k/reference_genome/human_g1k_v37.fasta \
-m- \
-Ov - \
| \
vt normalize - \
-r /run/media/gus/Storage/BCH/data/g1k/reference_genome/human_g1k_v37.fasta \
-o  CHB099.vt.norm.none.vcf
