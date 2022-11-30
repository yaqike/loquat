################################################################
################################################################
#####################                      #####################
#####################   枇杷实验室服务器     #####################
#####################                      #####################
################################################################
################################################################
#############                                      #############
#############   文中出现的____均为需要修改的参数      #############
#############                                      #############
################################################################
################################################################


### 一、QIIME 2_pipeline

##  0.准备工作

# 新建工作目录
mkdir /storage/wudi/qiime2/____

# 准备样本元信息metadata、原始数据seq/*.gz和manifest(需要自己手动编辑)
# 编辑完成后将metadata.txt、seq和manifest放入____文件夹

# 激活QIIME2工作环境
sudo su root
cd /storage/wudi/qiime2
docker run -t -i -v $(pwd):/data quay.io/qiime2/core:2022.2
cd ____

# 数据导入qiime2，格式为单端33格式
qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path manifest \
--output-path demux.qza \
--input-format SingleEndFastqManifestPhred33V2

##  1. 生成特征表和代表序列

# DADA2
# 支持多线程加速
# threads为需要的线程数，0则使用所有可用的内核
# trunc-len为截取的所需序列长度,0则不截取直接取全长
time qiime dada2 denoise-single \
--i-demultiplexed-seqs demux.qza \
--p-n-threads 0 \
--p-trunc-len 0 \
--o-representative-sequences rep-seqs.qza \
--o-table table.qza \
--o-denoising-stats dada2-stats.qza

# 特征表、代表序列、dada2去噪效果可视化
qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv \
--m-sample-metadata-file metadata.txt
	  
qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv
	  
qiime metadata tabulate \
--m-input-file dada2-stats.qza \
--o-visualization dada2-stats.qzv
    
##  2. 多样性分析

# 构建进化树用于多样性分析
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

# Alpha稀疏和深度选择
# p-max-depth参考table.qzv文件中Feature Count值(Sampling Depth的设置尽量包含大多数样本)
qiime diversity alpha-rarefaction \
--i-table table.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth ____ \
--m-metadata-file metadata.txt \
--o-visualization alpha-rarefaction.qzv

# 计算核心多样性 
# p-max-depth抽平参考table.qzv文件中Feature Count值(尽量包含大多数样本)
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth ____ \
--m-metadata-file metadata.txt \
--output-dir core-metrics-results

# Alpha多样性组间显著性分析和可视化
# 可选的index有faith_pd、shannon、observed_features、evenness
index=____
qiime diversity alpha-group-significance \
--i-alpha-diversity core-metrics-results/${index}_vector.qza \
--m-metadata-file metadata.txt \
--o-visualization core-metrics-results/${index}-group-significance.qzv

# Beta多样性组间显著性分析和可视化
# 可选的distance有unweighted_unifrac、bray_curtis、weighted_unifrac、jaccard
# 指定分组column可减少计算量(置换检验较耗时)
distance=____
column=____
qiime diversity beta-group-significance \
--i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
--m-metadata-file metadata.txt \
--m-metadata-column ${column} \
--o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
--p-pairwise

##  3. 物种组成分析

# 物种注释(classifier_gg_13_8_99_V5-V7.qza是V5-V7的训练文件)
qiime feature-classifier classify-sklearn \
--i-classifier /data/classifier_gg_13_8_99_V5-V7.qza \
--i-reads rep-seqs.qza \
--o-classification taxonomy.qza

# 可视化物种注释
qiime metadata tabulate \
--m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv

# 堆叠柱状图展示
qiime taxa barplot \
--i-table table.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file metadata.txt \
--o-visualization taxa-bar-plots.qzv

##  4. 差异分析ancom
# 为了保证分析的准确性，建议仅分析具体分类的某一对象(e.g.哪种序列变体在我们两个受试者的肠道样本中丰度存在差异)
# 确定分析的具体分类和对象
object=____
object_group=____
qiime feature-table filter-samples \
--i-table table.qza \
--m-metadata-file metadata.txt \
--p-where "[${object_group}]='${object}'" \
--o-filtered-table ${object}-table.qza

# 格式化特征表，添加伪计数
qiime composition add-pseudocount \
--i-table ${object}-table.qza \
--o-composition-table comp-${object}-table.qza

# 计算差异特征，指定需要分析的分组
column=____
time qiime composition ancom \
--i-table comp-${object}-table.qza \
--m-metadata-file metadata.txt \
--m-metadata-column ${column} \
--o-visualization ancom-${column}.qzv

# 按属水平合并
qiime taxa collapse \
--i-table ${object}-table.qza \
--i-taxonomy taxonomy.qza \
--p-level 6 \
--o-collapsed-table ${object}-table-l6.qza

# 格式化特征表，添加伪计数
qiime composition add-pseudocount \
--i-table ${object}-table-l6.qza \
--o-composition-table comp-${object}-table-l6.qza
  
# 计算差异属，指定分组类型比较
qiime composition ancom \
--i-table comp-${object}-table-l6.qza \
--m-metadata-file metadata.txt \
--m-metadata-column ${column} \
--o-visualization l6-ancom-${column}.qzv


### 二、PICRUSt2_pipeline

##  0.准备工作

# 新建工作目录并进入
# 将dada2降噪后rep-seqs.qza和table.qza后缀修改为rar
# 提取data文件夹中dna-sequences.fasta和feature-table.biom放入工作目录____
mkdir /storage/wudi/picrust2/____
cd /storage/wudi/picrust2/____

# 运行PICRUSt2
sudo su root
conda activate picrust2

# 计算每条序列的最近序列物种索引(NTSI)
picrust2_pipeline.py -s dna-sequences.fasta \
-i feature-table.biom \
-o picrust2_out_pipeline \
-p 16

# 结果注释
cd picrust2_out_pipeline
  
# 添加KO注释
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
-o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
  
# KEGG通路层级汇总(KEGG.Pathway.raw.txt文件为三级通路)
cd ..
zcat picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz > KEGG.KO.txt
python3 /storage/wudi/picrust2-2.5.0/summarizeAbundance.py \
-i KEGG.KO.txt \
-m /storage/wudi/picrust2-2.5.0/KO1-4.txt \
-c 2,3,4 -s ',+,+,' -n raw \
-o KEGG

# 提取KEGG.Pathway.raw.txt到本地并自行创建group.txt文件
# STAMP作图
# 作图前需除去PATH:备注（用excel处理即可，公式如下）
=IF(ISERROR(FIND(":ko",A2)),A2,LEFT(A2,LEN(A2)-14))



### 三、MicrobiomeAnalyst_pipeline

##  0.准备工作

# 新建工作目录
mkdir /storage/wudi/biom/____

# 提取qiime2中table.qza文件并修改后缀为rar
# 提取data文件夹中feature-table.biom放入工作目录____

# 运行虚拟环境
cd /storage/wudi/biom
source venv/bin/activate
cd /storage/wudi/biom/____

# 转换biom为经典格式
biom convert -i feature-table.biom -o table.txt --to-tsv

##  1.MicrobiomeAnalyst常规操作

# 根据Data Upload要求修改ASV table、Metadata file、Taxonomy table和phylogenetic tree
# 编辑table.txt标题行(加#NAME)得到asv_table.txt
# 编辑metadata.txt第一行(加#NAME)得到meta.txt
# 修改taxonomy.qza后缀为rar，提取data文件夹中taxonomy.tsv，利用Excel修改taxonomy.tsv得到taxa.txt(去掉 x__;去掉Confidence;第一行数据需要精确到Species)
# 修改rooted-tree.qza后缀为rar，提取data文件夹中tree.nwk

# 登录https://www.microbiomeanalyst.ca
# 详细教程可参考https://mp.weixin.qq.com/s/M1Pdr06Yo3APtCMq4w1hfQ

##  2.MicrobiomeAnalyst(Correlation Analysis)

# 根据Data Upload要求修改ASV table、Metadata file和Taxonomy table
# 可参考Example Datasets(Infected和UnInfected)



### 四、ggClusterNet_pipeline(仅在PC端进行测试，未在服务器测试)

##  0.准备工作

# 下载所需R包，下载后解压到R语言安装路径(下载链接：https://pan.baidu.com/s/1szwCjdvab2Tn9EXeU2BJnQ 提取码：peja)
# 查看R语言安装路径
.libPaths()

# 设置工作目录
setwd("E:ggClusterNet")
getwd()

# 加载所需R包
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)

##  1. 数据导入(依照Example Datasets修改标题行)
metadata = read.delim("2x_meta.tsv",row.names = 1)
otutab = read.delim("2x_asv_table.txt", row.names=1)
taxonomy = read.table("taxa.txt", row.names=1)

ps = phyloseq(sample_data(metadata),
              otu_table(as.matrix(otutab), taxa_are_rows=TRUE),
              tax_table(as.matrix(taxonomy))
              )

# 提取丰度最高的指定数量的otu进行构建网络
# N和r.threshold参数均可修改(此处采用N = 100r.threshold=0.90)
result = corMicro (ps = ps,
                   N = 100,
                   method.scale = "TMM",
                   r.threshold=0.90,
                   p.threshold=0.05,
                   method = "spearman"
)

# 提取相关矩阵
cor = result[[1]]
ps_net = result[[3]]

# 导出otu表格
otu_table = ps_net %>% 
  vegan_otu() %>%
  t() %>%
  as.data.frame()

tax = ps_net %>% vegan_tax() %>%
  as.data.frame()
tax$filed = tax$Phylum

##  2. 按照微生物分类不同设定分组(‘V3’默认为‘Phylum’，可修改)
group1 <- data.frame(ID = row.names(tax),group = tax$V3)
group1$group  =as.factor(group1$group)

# 计算布局
result1 = PolygonRrClusterG (cor = cor,nodeGroup =group1 )
node = result1[[1]]

# 计算节点并注释
tax_table = ps_net %>%
  vegan_tax() %>%
  as.data.frame()
nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
names(nodes)[names(nodes) == 'V3'] <- 'Phylum'

# 计算边
edge = edgeBuild(cor = cor,node = node)

##  3. 出图
pnet <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = as.factor(cor)),
                                data = edge, size = 0.5) +
geom_point(aes(X1, X2,fill = Phylum,size = mean),pch = 21, data = nodes) +
scale_colour_brewer(palette = "Set1") +
scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
theme(panel.background = element_blank()) +
theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
theme(legend.background = element_rect(colour = NA)) +
theme(panel.background = element_rect(fill = "white",  colour = NA)) +
theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

pnet

############################################################################
#################################          #################################
#################################   附录    #################################
#################################          #################################
############################################################################

### 〇、conda

##  基础操作

# 查看电脑上所有已创建的环境
conda info --env

# 删除不需要的环境
conda remove -n 环境名 --all

### 一、QIIME 2

##  0.环境搭建

# 安装docker
sudo apt-get install docker-ce=18.06.3~ce~3-0~ubuntu

# 查看版本号
docker version

# 安装qiime2
docker pull quay.io/qiime2/core:2022.2

# 退出docker
exit 

##  1.训练分类器

# 导入数据库文件并解压
tar -zxvf gg_13_8_otus.tar.gz

# 导入参考序列
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path gg_13_8_otus/rep_set/99_otus.fasta \
--output-path 99_otus.qza

# 导入物种分类信息
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path gg_13_8_otus/taxonomy/99_otu_taxonomy.txt \
--output-path ref-taxonomy.qza

# 引物提取参考序列的扩增区段(引物序列可自行修改，此处为V5-V7序列)
time qiime feature-classifier extract-reads \
--i-sequences 99_otus.qza \
--p-f-primer AACMGGATTAGATACCCKG \
--p-r-primer ACGTCATCCCCACCTTCC \
--o-reads ref-seqs.qza

# 基于筛选的指定区段，生成实验特异的分类器
time qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier classifier_gg_13_8_99_V5-V7.qza



### 二、PICRUSt2

##  0.环境搭建

# mamba安装
conda install mamba -n base -c conda-forge

# 安装python3.9以下版本并运行
conda create --name py37 python=3.7
conda activate py37

# 解压picrust2-2.5.0.tar.gz并进入文件夹
tar xvzf picrust2-2.5.0.tar.gz
cd picrust2-2.5.0/

# 创建并激活环境
mamba env create -f picrust2-env.yaml
conda activate picrust2
sudo pip install --editable .

##  1.PICRUSt2分析可选项

# 添加EC注释
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
-o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

# 添加pathway注释(基于MetaCyc数据库)
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
-o pathways_out/path_abun_unstrat_descrip.tsv.gz



### 三、MicrobiomeAnalyst

##  0.环境搭建

# 进入工作目录
cd /storage/wudi/biom

# 创建虚拟环境
python3 -m venv venv

# 安装biom包
pip install biom-format

# 安装biom2.0格式支持
pip install h5py


