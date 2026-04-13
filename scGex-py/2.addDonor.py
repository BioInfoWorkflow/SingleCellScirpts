# 获取barcode, 进行拆donor

# region |- Import -|
from distutils.version import LooseVersion
import scanpy as sc
import pandas as pd
import os
# endregion

import argparse


exam=['-i', 'path/to/Merge.h5ad', 
    '-d', 'path/to/demuxlet/mappings/',] # 基于SNP的demuxlet 获取单个sub bam内 barcode->donor的映射表 目录, 内部的 sub名称 要和 10x原始dir subid名字一致

# 规定项目
BAR_COL_W = 0.8 # bar图的单个item宽度
QUANT_QCPLOT = 0.985 # 985分位数 作 QCplot的cutoff
VLN_COL_W = 2   # vln plot 单个item宽度

parser = argparse.ArgumentParser(description="处理输入目录和输出文件路径")
parser.add_argument("-i", "--input", type=str, help="Merge的h5ad文件")
parser.add_argument("-d", "--donordir", type=str, required=True, help="输出的barcodes - donor 根目录")
# parser.add_argument("-s", "--species", type=str, nargs="+", default=None, help="输出Merge的h5ad, 根据物种基因组mapping率 分发")
args = parser.parse_args(exam)

adata = sc.read(args.input)
# SNP merge
subgroups = tuple(f for f in os.listdir(args.donordir) if f.endswith('.tsv'))
snpMeta = []
for sub in subgroups:
    tmp = pd.read_csv(f"{args.donordir}/{sub}", sep='\t', index_col=0)
    # tmp.index = f"UDA-2-Gonad-1-{sub.split('.tsv')[0]}." + tmp.index
    snpMeta.append(tmp)
    print(f"\r{sub}", end="")

# merge
snpMeta = pd.concat(snpMeta)

# meta info
bcs = snpMeta.index.intersection(adata.obs_names)
snpMeta = snpMeta[[ 'NUM.SNPS', 'DROPLET.TYPE', 'SNG.POSTERIOR', 'Donor']]
for k in snpMeta.columns:
    adata.obs[k] = 'Filtered' if snpMeta[k].dtype.str == '|O' else -1
    adata.obs.loc[bcs, k] = snpMeta.loc[bcs, k].values.copy()


print(adata.obs['Donor'].value_counts())
adata=adata[adata.obs['Donor']!='Filtered'] # 移除咯~

adata.write(args.input)
# obs meta
obs_dir = os.path.dirname(args.input)
adata.obs.to_csv(f'{obs_dir}/obs.tsv', sep='\t')

donors = sorted(adata.obs['Donor'].unique().tolist())
# donors.remove("Filtered")

# plots
plotData = pd.DataFrame(index=donors, columns=['cell', 'median Gene', 'median UMI'])
for donor in donors:
    tmp = adata[adata.obs['Donor']==donor]
    plotData.loc[donor] = [tmp.n_obs, tmp.obs['n_genes_by_counts'].median(), tmp.obs['total_counts'].median()]


# plot
import os
import scanpy as sc
from matplotlib import pyplot as plt
import numpy as np

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


os.chdir(obs_dir)
os.makedirs("figures", exist_ok=True)
pdata = adata[(adata.obs['n_genes_by_counts']<np.quantile(adata.obs['n_genes_by_counts'], QUANT_QCPLOT)) & 
            (adata.obs['total_counts']<np.quantile(adata.obs['total_counts'], QUANT_QCPLOT)) ]
vlnQCPlot(pdata, 'Donor', groupby='Donor')
plt.close("all")


# plots
donors = adata.obs['Donor'].unique().tolist()
plotData = pd.DataFrame(index=sorted(donors, key=LooseVersion), columns=['cell', 'median Gene', 'median UMI'])
for donor in donors:
    tmp = adata[adata.obs['Donor']==donor]
    plotData.loc[donor] = [tmp.n_obs, tmp.obs['n_genes_by_counts'].median(), tmp.obs['total_counts'].median()]


plt.close("all")
# colors = plt.cm.gist_rainbow(np.linspace(0, 1, plotData.shape[0]))
# 获取颜色映射
cmap = plt.cm.get_cmap('tab20')
# 生成循环使用的颜色列表
colors = cmap(np.linspace(0, 1, plotData.shape[0]))

k='cell'
fig, ax = plt.subplots(figsize=(BAR_COL_W*len(donors), 6))
bars = ax.bar(plotData.index, plotData[k], color=colors,
            linewidth=1, edgecolor='black')
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{height}',
             ha='center', va='bottom',
             fontsize=8,rotation=45)

# 旋转x轴标签
plt.xticks(rotation=45, ha='right')
ax.set_title(k)
# plt.bar_label(ax, labels=plotData[k])
plt.savefig(f"figures/Homo.Donor.{k.replace(' ','_')}.png", bbox_inches='tight', dpi=300)
plt.close("all")

k='median Gene'
fig, ax = plt.subplots(figsize=(BAR_COL_W*len(donors), 6))
bars = ax.bar(plotData.index, plotData[k], color=colors,
            linewidth=1, edgecolor='black')
# 旋转x轴标签
plt.xticks(rotation=45, ha='right')
ax.set_title(k)
# plt.bar_label(ax, labels=plotData[k])
plt.savefig(f"figures/Homo.Donor.{k.replace(' ','_')}.png", bbox_inches='tight', dpi=300)
plt.close("all")

k='median UMI'
fig, ax = plt.subplots(figsize=(BAR_COL_W*len(donors), 6))
bars = ax.bar(plotData.index, plotData[k], color=colors,
            linewidth=1, edgecolor='black')
# 旋转x轴标签
plt.xticks(rotation=45, ha='right')
ax.set_title(k)
# plt.bar_label(ax, labels=plotData[k])
plt.savefig(f"figures/Homo.Donor.{k.replace(' ','_')}.png", bbox_inches='tight', dpi=300)
plt.close("all")