# 获取barcode, 进行拆donor

# region |- Import -|
import scanpy as sc
import os
import pandas as pd
# endregion

import argparse


exam=['-i', 'path/to/Megre.h5ad', 
    '-b', 'path/to/cellranger.outs/',
    '-o', 'path/to/sub.bcdir'] # 输出 每个子样品的 barcode list

parser = argparse.ArgumentParser(description="处理输入目录和输出文件路径")
parser.add_argument("-i", "--input", type=str, help="Merge的h5ad文件")
parser.add_argument("-o", "--outdir", type=str, required=True, help="输出的subdir barcodes主目录")
parser.add_argument("-b", "--bam10xdir", type=str, help="cellranger输出目录")
parser.add_argument("-k", "--key", type=str, default='SubLib', help="SubLib obs key名称")
parser.add_argument("-bc", "--barcodekey", type=str, default='Raw.Barcode', help="SubLib obs key名称")
# parser.add_argument("-s", "--species", type=str, nargs="+", default=None, help="输出Merge的h5ad, 根据物种基因组mapping率 分发")
args = parser.parse_args()

adata = sc.read(args.input)
subgroups = sorted(adata.obs[args.key].unique().tolist())
for sub in subgroups:
    barcode = adata[adata.obs[args.key]==sub].obs[args.barcodekey]
    barcode.to_csv(f'{args.outdir}/{sub}.txt', index=False, header=False)
    print(f"\r{sub}", end='')



# 生成 Bam Barcode 配对表
pairtab = pd.DataFrame(index=subgroups, columns=['Bam', 'Barcode', 'BarcodePrefix'])
pairtab.index.name = 'SubID'

pairtab['Bam'] = args.bam10xdir + '/' + pairtab.index + '/outs/gex_possorted_bam.bam'
pairtab['Barcode'] = args.outdir + '/' + pairtab.index + '.txt'
pairtab['BarcodePrefix'] = pairtab.index + '.'

pairtab.to_csv(f'{args.outdir}/../pair.tsv', sep='\t')