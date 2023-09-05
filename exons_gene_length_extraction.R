####=====================================================================================================================####
####=====================================================================================================================####
#### 1. ��ΰѴ� TCGA ���ص� HTseq count ת���� FPKM ####


#### 1�����ȿ�һ��ʲô��FPKM��
# FPKM��fragments per kilobase million����fragment per kilobase of transcript per million mapped reads����ÿһ�����map�ϵ�reads��map�������ӵ�ÿһǧ������ϵ�fragments������

# ��ʽ���£�
# FPKM = read counts / (mapped reads(Millions)* exon length) #�����exon length��λ�� kb��

# ����Ҫ�����ص���ÿ�� exon �ĳ��ȣ�Ϊ���󳤶ȣ�����Ҫ����ÿ��������ע���ļ����������������أ�
# https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files������ GDC ��������ַ����ȥ�Ժ�ѡ������ļ��������أ�
# Annotation Files --- GDC.h38 Flattened GENCODE v22 GFF (used in RNA-Seq alignment and by HTSeq)
#                          genecode.v22.annotation.gtf.gz
#                            md5: 291330bdcff1094bc4d5645de35e0871
#                            file size: 39.0 MB

# ���غ��ѹ���á�


# chooseBioCmirror() #һ��ѡ�񱱾�����
# chooseCRANmirror() #һ��ѡ�����ݾ���
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor") #�廪����


# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("GenomicFeatures")

# install.packages("tidyverse")



#### 2������exon���ȣ����ұ����dataframe ####

library(GenomicFeatures) #��ȡexon��������
setwd("E:/BaiduNetdiskDownload/R_data_Analysis/2.Papers/Analysis/AML/Analysis_Whole/1.dataPreProcess")
txdb <- makeTxDbFromGFF("gencode.v22.annotation.gtf", format="gtf") #��ȡ�������ص�ע���ļ�

#ͨ�� exonsBy ��ȡÿ�������ϵ����������ӵ���ʼλ�����ֹλ�㣬Ȼ���� reduce ȥ�����ص�����Ĳ��֡�
#�����㳤��
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene, function(x){sum(width(reduce(x)))})
#����������̫���ˣ�ֱ����ȡע���ļ���Ļ�������ensemble����exon������Ϣ

View(exons_gene_lens) # �õ�������Ⱥ󣬴򿪿�һ�����ݼ���ʲô����
head(exons_gene_lens) # �� view ������ѡһ

#�鿴�ж��ٻ�������
length(exons_gene_lens)

#��һ����������
class(exons_gene_lens) #�鿴���֣��Ǹ��б����͵����ݣ�[1] "list"����Ҫת����dataframe��������

#��������ת����list ---> dataframe
exons_gene_lens1 <- as.data.frame(exons_gene_lens)
dim(exons_gene_lens1)#�鿴ת���������ά��
str(exons_gene_lens1)#�鿴���ݽṹ
class(exons_gene_lens1)#�鿴��������

#ת��֮����������Ϊ������������Ҫ��һ��ת����
exons_gene_lens1 <- t(exons_gene_lens1)
dim(exons_gene_lens1)#�鿴����ά��
head(exons_gene_lens1)#�鿴���ݳ�ʲô������ʱ����Ϊ��������

#ɾ����һ�У�
exons_gene_lens2 <- exons_gene_lens1
View(exons_gene_lens2)
exons_gene_lens2 <- as.data.frame(exons_gene_lens2)
rownames(exons_gene_lens2) #�鿴������ʲô����
head(exons_gene_lens2)
dim(exons_gene_lens2)


#---------------------------------------------------------------------------------------------------#
#�ѻ���С������������ȥ����

rownames(exons_gene_lens2) #��ȡ����
xc <- gsub("\\.(\\.?\\d*)", "", rownames(exons_gene_lens2)) #ɾ��С����
xc #�鿴����
rownames(exons_gene_lens2) = xc #��exon �б�������ͬ�� xc����ֱ���滻exon�������������С���㡣
#---------------------------------------------------------------------------------------------------#


#�� exon ���ȱ�����������Ϊ��Length��
colnames(exons_gene_lens2) = "Length"
View(exons_gene_lens2) #�鿴���ݡ�
write.table(exons_gene_lens2,file = "1.exons_gene_length.txt",sep = "\t",quote = F)
##=======================================================================================================##
#�����ļ�����Excel���ļ�������ͷ����һ�£���һ�м��ϡ�id������Length��������2�У�����������������󱣴档
##=======================================================================================================##

rm(list=ls()) #һ�����

####==========================================================================================================####
####==========================================================================================================####