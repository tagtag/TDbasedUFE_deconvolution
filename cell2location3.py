import scanpy as sc
import scvi
import anndata
from cell2location.models import RegressionModel, Cell2location
import numpy as np
import pandas as pd
from scipy.io import mmread
import pickle

# ----------------------------
# 1. シングルセルデータの読み込み
# ----------------------------

# 遺伝子発現行列の読み込み（行が遺伝子、列が細胞）

adata_sc  = sc.read_csv("GSM8341772_control_3DE_genematrix.csv.gz")
adata_sc  = anndata.AnnData(X=adata_sc.X.T, obs = adata_sc.var, var=adata_sc.obs)


# ----------------------------
# 2. メタデータの統合
# ----------------------------
metadata_path = 'results3.tsv'
metadata = pd.read_csv(
  metadata_path, 
  sep='\t', 
  header=None, 
  index_col=0,
  names=['cell_label', 'cell_type']
)
adata_sc.obs['cell_type'] = metadata.loc[adata_sc.obs_names, 'cell_type']

# ----------------------------
# 3. シングルセルデータ前処理
# ----------------------------
# 基本的なフィルタリング
sc.pp.filter_genes(adata_sc, min_counts=1)
sc.pp.filter_cells(adata_sc, min_counts=1)

# 正規化前の生データを保存
adata_sc.layers['counts'] = adata_sc.X.copy()

# 正規化と対数変換
sc.pp.normalize_total(adata_sc, target_sum=1e4)
sc.pp.log1p(adata_sc)
adata_sc.raw = adata_sc  # 生データをバックアップ

# ----------------------------
# 4. リファレンスモデルのトレーニング
# ----------------------------
# モデル設定
RegressionModel.setup_anndata(
  adata=adata_sc,
  labels_key='cell_type',
  layer='counts'
)

# モデル初期化
mod = RegressionModel(adata_sc)

# トレーニング
mod.train(
  max_epochs=300,
  batch_size=2500,
  train_size=1,
  lr=0.002,
  accelerator="cpu"
)



# 結果の抽出
adata_sc = mod.export_posterior(
  adata_sc,
  sample_kwargs={
    'num_samples': 1000,
    'batch_size': 2500
  }
)
adata_sc.var_names_make_unique()


# リファレンスシグネチャの抽出
adata_ref = adata_sc.copy()
factor_names = adata_ref.uns['mod']['factor_names']
adata_ref.varm['means'] = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in factor_names]].copy()

# ----------------------------
# ----------------------------
# 5. Visiumデータの読み込み（H5ファイルなし版）
# ----------------------------

# 遺伝子発現行列の読み込み
vis_matrix_path = 'GSM8341553_sample3_matrix.mtx.gz'
vis_genes_path = 'GSM8341553_sample3_features.tsv.gz'
vis_barcodes_path = 'GSM8341553_sample3_barcodes.tsv.gz'

# 疎行列の読み込みと転置（スポット×遺伝子形式に変換）
vis_matrix = mmread(vis_matrix_path).T.tocsr()

# 遺伝子名とバーコードの読み込み
vis_genes = pd.read_csv(vis_genes_path, sep='\t', header=None)[0].values  # 第2列が遺伝子名
vis_barcodes = pd.read_csv(vis_barcodes_path, sep='\t', header=None)[0].values

# AnnDataオブジェクトの作成
adata_vis = anndata.AnnData(
  X=vis_matrix,
  obs=pd.DataFrame(index=vis_barcodes),
  var=pd.DataFrame(index=vis_genes)
)
adata_vis.var_names = vis_genes
adata_vis.var_names_make_unique()
# ----------------------------
# 6. 空間座標データの読み込み
# ----------------------------
tissue_positions_path = 'GSM8341553_sample3_tissue_positions_list.csv.gz'

# 座標データの読み込み（標準的な10X形式のCSV）
spatial_coords = pd.read_csv(
  tissue_positions_path,
  header=None,
  names=[
    'barcode', 'in_tissue', 'array_row', 'array_col',
    'pxl_row_in_fullres', 'pxl_col_in_fullres'
  ]
)

# 座標データをAnnDataに統合
spatial_coords.index = spatial_coords['barcode']
adata_vis.obsm['spatial'] = spatial_coords.loc[adata_vis.obs_names, [
  'pxl_col_in_fullres', 'pxl_row_in_fullres'
]].values.astype(int)

# ----------------------------
# 6. Visiumデータ前処理
# ----------------------------
# 基本的なフィルタリング
sc.pp.filter_genes(adata_vis, min_counts=10)
sc.pp.filter_cells(adata_vis, min_counts=5)

# 遺伝子の整合性確認
intersected_genes = np.intersect1d(adata_vis.var_names, adata_ref.var_names)
adata_vis = adata_vis[:, intersected_genes].copy()
adata_ref = adata_ref[:, intersected_genes].copy()

# ----------------------------
# 7. cell2locationモデルの適用
# ----------------------------
# モデル設定
Cell2location.setup_anndata(adata_vis)

# モデル初期化
mod_vis = Cell2location(
  adata_vis,
  cell_state_df=adata_ref.varm['means'],
  N_cells_per_location=30,
  detection_alpha=200
)

# トレーニング
mod_vis.train(
  max_epochs=30000,
  batch_size=adata_vis.n_obs,
  train_size=1
)

# 結果の抽出
adata_vis = mod_vis.export_posterior(
  adata_vis,
  sample_kwargs={
    'num_samples': 1000,
    'batch_size': mod_vis.adata.n_obs
  }
)

# ----------------------------
# 8. 結果の保存と可視化
# ----------------------------
# 細胞分布データの保存
cell_abundance = adata_vis.obsm['q05_cell_abundance_w_sf']
cell_abundance.to_csv('cell_abundances3.csv')

