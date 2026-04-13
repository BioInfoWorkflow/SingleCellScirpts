#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ scRNA.py   2025/4/28/-12:04 
╰───────────────────────────────────────╯ 
│ Description:
    用于RNA数据  预处理  
""" # [By: HuYw]

# region |- Import -|
import hdf5plugin
import scanpy as sc
import pandas as pd
import numpy as np
import os
# endregion

import argparse

exam=['-i', 'path/to/Megre.h5ad', 
'-o', 'path/to/analysis/Gex',  # 输出路径
'-sp', 'hg', # 物种
# '-BK', 'library_id', # 如果没有 不写此参数, 不进行 hvg的批次矫正
# '-MK', 'celltype_l1',# 如果umap想绘制更多 obs key, 在这续写
 ] 

parser = argparse.ArgumentParser(description="处理输入目录和输出文件路径")
parser.add_argument("-i", "--input", type=str, help="输入Merge Raw h5ad")
parser.add_argument("-o", "--outdir", type=str, required=True, help="输出Scanpy Run的总目录")
parser.add_argument("-BK", "--batchKey", type=str, default=None, help="obs.Batch key")
parser.add_argument("-MK", "--moreKey", type=str, nargs='+', default=[], required=False, help="想绘制的更多 more obs key")
parser.add_argument("-n", "--nhvg", type=int, default=3000, help="hvg数目, 3000")
parser.add_argument("-mg", "--mingene", type=int, default=300, help="cell最小基因测量数")
parser.add_argument("-mc", "--mincell", type=int, default=3, help="gene最小细胞表达数")
parser.add_argument("-mt", "--mtCutoff", type=int, default=10, help="mt cutoff 百分比, 默认是10")
parser.add_argument("-sp", "--species", type=str, default="hg", help="物种 (hg/mm) (人类/小鼠)")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

INFO=pd.DataFrame(index=['n.Cell', 'Median.Gene', 'Median.UMI', 'Median.MT', 'n.Cluster',
                         'n.Phase-G1', 'n.Phase-G2M', 'n.Phase-S'], columns=['Value'])
INFO['Value'] = -1

adata = sc.read(args.input)
if not args.species.lower().startswith("mix") or (adata.n_vars > 50000 or adata.var_names.str.contains(f'{args.species}_').sum() > 3000):
    adata = adata[:, adata.var_names.str.startswith(f'{args.species}_')]
    adata.var_names = adata.var_names.str.replace(f'{args.species}_', '')
    if 'gene_ids' in adata.var:
        adata.var['gene_ids'] = adata.var['gene_ids'].str.replace(f'{args.species}_', '')


n_obs = adata.n_obs


mtprefix = 'MT-' if args.species.lower().startswith('hg') else 'mt-'
adata.var["mt"] = adata.var_names.str.startswith(mtprefix)  # 线粒体基因
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, log1p=True
)
adata = adata[adata.obs['pct_counts_mt'] < args.mtCutoff, :]
print(f"Filter MT rate from {n_obs} -> {adata.n_obs}")

sc.pp.filter_cells(adata, min_genes=args.mingene)
sc.pp.filter_genes(adata, min_cells=args.mincell)
print(adata)
INFO.loc['n.Cell', 'Value'] = adata.n_obs
INFO.loc['Median.Gene', 'Value'] = round(adata.obs['n_genes'].median())
INFO.loc['Median.UMI', 'Value'] = round(adata.obs['total_counts'].median())
INFO.loc['Mean.MT', 'Value'] = round(adata.obs['pct_counts_mt'].mean())

# 必须是原始矩阵才可以
assert adata.X.max().is_integer(), f"[{args.input}] H5ad file not raw counts, job quit..."



# 保留raw counts
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata)

# 筛选HVG
sc.pp.highly_variable_genes(
	adata,
	n_top_genes=args.nhvg,
	layer='counts',
	flavor="seurat_v3",
	batch_key=args.batchKey,
)

# 细胞周期 基因集来自Seurat
s_genes = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1",
    "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"
]
g2m_genes = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "SMC4",
    "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B",
    "GTSE1", "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C",
    "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
    "G2E3", "GAS2L3", "CBX5", "CENPA"
]
# 细胞周期小鼠化
if args.species.lower().startswith("mm"):
    s_genes = [g[0].upper() + g[1:].lower() for g in s_genes]
    g2m_genes = [g[0].upper() + g[1:].lower() for g in g2m_genes]


sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
phase_df=adata.obs['phase'].value_counts().astype(int)
phase_df.index = 'n.Phase-' + phase_df.index 
INFO.loc[phase_df.index, 'Value']=phase_df.values
INFO['Value']=INFO['Value'].astype(int)

# 自动确定 PCA 的主成分数量
def calculate_pca_components(adata, fast=True):
    sc.tl.pca(adata, n_comps=50) if fast else sc.tl.pca(adata, n_comps=200)
    explained_variance_ratio = adata.uns['pca']['variance_ratio']
    # 使用“肘部”方法确定主成分数量
    n_components = 0
    # 如果没有找到合适的主成分数量，使用默认值
    if fast:
        n_components = min(adata.n_obs - 1, adata.n_vars - 1, 15)
    else:
        for n in range(1, len(explained_variance_ratio) - 1):
            if (explained_variance_ratio[n + 1] - explained_variance_ratio[n]) / (explained_variance_ratio[0] - explained_variance_ratio[-1]) < explained_variance_ratio[-1]:
                n_components = n + 1
                break
    return n_components


n_pcs=calculate_pca_components(adata, fast=True)
sc.pl.pca_variance_ratio(adata, log=False)
from matplotlib import pyplot as plt
plt.savefig(f"{args.outdir}/PCA.VarRatio.png", bbox_inches='tight', dpi=300)
# 自动确定聚类参数
def calculate_clustering_params(n_cells):
    if n_cells < 200:
        resolution = 2.5
        neighbors = 10  # 默认值，可以根据需要调整
    else:
        neighbors = int(2.7 * np.log(n_cells) - 10.21)
        resolution = neighbors / (4 * np.log(n_cells) - 19.25)
    return neighbors, resolution


neighbors, resolution = calculate_clustering_params(adata.n_obs)
print(f"Calculated neighbors: {neighbors}, resolution: {resolution}")

# 2. 构建 KNN 图
sc.pp.neighbors(adata, n_neighbors=neighbors, n_pcs=n_pcs)
# 3. Leiden 聚类
sc.tl.leiden(adata, resolution=resolution)
sc.tl.umap(adata)


INFO.loc['n.Cluster', 'Value'] = len(adata.obs['leiden'].unique())
INFO.to_csv(f"{args.outdir}/Metrics.tsv", sep='\t')
adata.obs.to_csv(f"{args.outdir}/obs.tsv", sep='\t')

# QC plots umap
obs_qc = ["n_genes", "total_counts", "pct_counts_mt"]
umap = sc.pl.umap(
    adata,
    color=obs_qc,
    # Setting a smaller point size to get prevent overlap
    size=5,
    ncols=3,
    return_fig=True,
    wspace=0.5,
)
umap.savefig(f"{args.outdir}/UMAP.QC.png", bbox_inches='tight', dpi=300)

# cate plost umap
obs_cate=["leiden", "phase"] + args.moreKey
if args.batchKey:
    obs_cate += [args.batchKey]


umap = sc.pl.umap(
    adata,
    color=obs_cate,
    # Setting a smaller point size to get prevent overlap
    size=5,
    return_fig=True,
    ncols=3,
    wspace=0.5,
)
umap.savefig(f"{args.outdir}/UMAP.Cate.png", bbox_inches='tight', dpi=300)


# deg
import pandas as pd
obs_key='leiden'
# 计算差异基因
sc.tl.rank_genes_groups(adata, groupby=obs_key, method="wilcoxon")
# standard_scale=None 才显示真实的norm值, 否则是rank百分数
fig_dotplot = sc.pl.rank_genes_groups_dotplot(
    adata, groupby=obs_key, standard_scale=None, n_genes=5, return_fig=True
)
# categories_order 指定dotplot顺序
# dendrogram=False 顺序才生效
fig_dotplot.savefig(f"{args.outdir}/DEG.{obs_key}.png", dpi=300, bbox_inches='tight')
# 获取所有组高变基因数据框
deg_groups = adata.uns["rank_genes_groups"]["names"].dtype.names
deg_df = []
for g in deg_groups:
    tmpdf = sc.get.rank_genes_groups_df(adata, group=g)
    tmpdf["group"] = g
    tmpdf=tmpdf[['names', 'group', 'scores', 'logfoldchanges', 'pvals_adj']]
    deg_df.append(tmpdf)


# 导出deg
deg_df = pd.concat(deg_df, ignore_index=True)
deg_df.to_csv(f"{args.outdir}/DEG.{obs_key}.tsv", sep='\t')
# 输出head5
print("Head 5 DEG for each cluster")
for g in deg_df['group'].unique():
    print(f'{g}: {deg_df[deg_df["group"]==g].sort_values("scores").head(5)["names"].tolist()}')



deg_df.sort_values(['scores'], ascending=False).groupby('group').head(5)

# 保存
adata.write(f"{args.outdir}/adata.h5ad")