#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# # # # # # # # # # # # 
"""
╭───────────────────────────────────────╮ 
│ scanpyMerge.py   2025/04/27/-13:56 
╰───────────────────────────────────────╯ 
│ Description:
    用于UDA 96孔合并
""" # [By: HuYw]

# region |- Import -|
import scanpy as sc
from distutils.version import LooseVersion
from scipy import sparse
import os
import scrublet
import pandas as pd
# endregion

import argparse


exam=['-i', 'path/to/cellranger.outs/', 
'-o', 'path/to/Merge.h5ad']

CHOOSE_MATRIX='outs/raw_feature_bc_matrix' # 使用raw 未filter的原始矩阵
CHOOSE_BAM='outs/possorted_genome_bam.bam' # 如果使用SNP拆分, 这里多输出一个 barcodes list -> bam 的映射表, 不需要不用管
VLN_COL_W = 5   # vln plot 单个item宽度, 在这没啥用

parser = argparse.ArgumentParser(description="处理输入目录和输出文件路径")
parser.add_argument("-i", "--input", type=str, help="输入的10x输出主目录, 其下有各个Sub输出子目录")
parser.add_argument("-o", "--output", type=str, required=True, help="输出Merge的h5ad文件路")
parser.add_argument("-mt", "--mtprefix", type=str, default="MT-", help="MT 基因开头")
parser.add_argument("-mg", "--mingene", type=int, default=100, help="最小基因测量数")
# parser.add_argument("-s", "--species", type=str, nargs="+", default=None, help="输出Merge的h5ad, 根据物种基因组mapping率 分发")
args = parser.parse_args(exam)

subdir = sorted(os.listdir(args.input), key=LooseVersion)
print(f"Find {len(subdir)} Subs in dir: {args.input}")

adata = []
for sub in subdir:
    dir_10x = f"{args.input}/{sub}/{CHOOSE_MATRIX}"
    if not os.path.isdir(dir_10x): 
        print(f"Warning: SubDir {sub} is not a standard 10x outputs! jump this sub dir! [{dir_10x}]")
        continue
    sub_data = sc.read_10x_mtx(dir_10x)
    # 根据设定筛除低质量细胞
    sc.pp.filter_cells(sub_data, min_genes=args.mingene)
    # 计算线粒体表达比例
    sub_data.var["mt"] = sub_data.var_names.str.startswith(args.mtprefix)  # 线粒体基因
    sc.pp.calculate_qc_metrics(
        sub_data, qc_vars=["mt", ], inplace=True, log1p=False
    )
    sub_data.obs['SubLib'] = sub
    sub_data.obs['Raw.Barcode'] = sub_data.obs_names
    sub_data.obs_names = sub + "." + sub_data.obs_names
    # doublet
    sub_data.X = sparse.csr_matrix(sub_data.X, dtype=int)
    # safe scrublet
    scrub = scrublet.Scrublet(sub_data.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    sub_data.obs['Doublet.Score'] = doublet_scores
    sub_data.obs['Doublet.Flag'] = predicted_doublets
    print(f"SubDir {sub} Read {sub_data.n_obs} cells")
    adata.append(sub_data)


adata = sc.concat(adata, axis=0, join='outer', label='concat_id', keys=subdir, merge='first')
adata.write(args.output)
# export barcodes
bdir=os.path.dirname(args.output)
jobname=os.path.basename(args.output).rsplit(".", 1)[0]
os.makedirs(f"{bdir}/barcodes/{jobname}", exist_ok=True)
# pair tsv
pairtab = pd.DataFrame(columns=['SubID','Bam','Barcode','BarcodePrefix'])
from tqdm import tqdm
for sub in tqdm(adata.obs['SubLib'].unique()):
    tmp = adata[adata.obs['SubLib']==sub]
    bampath=f"{args.input}/{sub}/{CHOOSE_BAM}"
    bcpath=f"{bdir}/barcodes/{jobname}/{sub}.tsv"
    tmp.obs['Raw.Barcode'].to_csv(bcpath, sep='\t', index=False, header=False)
    pairtab.loc[sub] = [sub, bampath, bcpath, sub + "."]
    if not os.path.isfile(bampath): print(f"Warning Bam No existed: {sub} [{bampath}]\n")
    if not os.path.isfile(bcpath): print(f"Warning BCtsv No existed: {sub} [{bcpath}]\n")


pairtab.to_csv(f"{bdir}/barcodes/{jobname}.pair.tsv", sep='\t', index=False, header=True)



# QC
import pandas as pd
# 质控报告 梯度
qcs=[
[0, 0],
[200, 300],
[300, 500],
[500, 800],
[600, 1000],
[800, 1200],
[1000, 1500]
]
QC=pd.DataFrame(qcs, columns=('min.Gene', 'min.UMI'))
QC['nCell']=0
QC['nCell.withDoublet']=0
QC['median.Gene']=0
QC['median.UMI']=0
for i, r in QC.iterrows():
	tmp = adata[(adata.obs['n_genes']>=r['min.Gene']) & (adata.obs['total_counts']>=r['min.UMI'])]
	QC.loc[i, 'nCell']=tmp[~tmp.obs['Doublet.Flag']].n_obs
	QC.loc[i, 'nCell.withDoublet']=tmp.n_obs
	QC.loc[i, 'median.Gene']=tmp.obs['n_genes'].median()
	QC.loc[i, 'median.UMI']=tmp.obs['total_counts'].median()


# 输出
QC.to_csv(f"{bdir}/{jobname}.QC_grad.csv", index=False)
# 基本过滤 Violin 图
from matplotlib import pyplot as plt
def vlnQCPlot(adata, prefix, groupby='SubLib'):
	sub_width= VLN_COL_W * adata.obs[groupby].unique().shape[0]
	plt.close('all')
	sc.pl.violin(
	adata,
	"n_genes_by_counts",
	groupby=groupby,
	size=0,
	show=False,
	rotation=45,
	inner='quart',
	save=f"{prefix}.{groupby}-Gene.png"
	)
	plt.close('all')
	sc.pl.violin(
	adata,
	"total_counts",
	groupby=groupby,
	size=0,
	show=False,
	rotation=45,
	inner='quart',
	save=f"{prefix}.{groupby}-UMI.png",
	)
	plt.close('all')
	sc.pl.violin(
	adata,
	"pct_counts_mt",
	groupby=groupby,
	size=0,
	show=False,
	rotation=45,
	inner='quart',
	save=f"{prefix}.{groupby}-MT.png"
	)
	plt.close('all')
	sc.pl.violin(
	adata,
	"Doublet.Score",
	groupby=groupby,
	size=0,
	show=False,
	rotation=45,
	inner='quart',
	save=f"{prefix}.{groupby}-Doublet.png",
	)
	plt.close('all')


sc.settings.set_figure_params(dpi=300, facecolor="white") # 美化
os.chdir(bdir)
adata.obs['Whole'] = 'Whole' # 临时元信息, 用于绘制QC全局图像
vlnQCPlot(adata, f'{jobname}.QC', groupby='Whole')
vlnQCPlot(adata, f'{jobname}.QC', groupby='SubLib') # 单样本展示