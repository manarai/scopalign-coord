"""Redraw poster figures as editable vector SVG/PDF from saved arrays (no refits)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc, anndata as ad

REPO = Path('/Users/terooatt/Documents/scopalign-coord')
DATA = Path('/Users/terooatt/Documents/scopalign-manus/data')
OUT = REPO / 'poster' / 'figures'; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'svg.fonttype': 'none', 'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['IBM Plex Sans', 'Helvetica', 'Arial'], 'font.size': 16,
                     'axes.titlesize': 18, 'axes.spines.top': False, 'axes.spines.right': False})
BATCH = ('#1f5fa8', '#e07b00')
TYPES = ['B cells', 'DC', 'ILC', 'Megakaryocytes/platelets', 'Monocytes', 'T cells', 'pDC']
TCOL = dict(zip(TYPES, ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2']))

def save(fig, name):
    for ext in ('svg', 'pdf'):
        fig.savefig(OUT / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)

def scatter_pair(axes, xy, batch, types, labels, title):
    for (b, lab), col in zip(labels, BATCH):
        m = batch == b
        axes[0].scatter(xy[m, 0], xy[m, 1], s=3, alpha=.7, color=col, label=lab, linewidths=0)
    for t in sorted(set(types)):
        m = types == t
        axes[1].scatter(xy[m, 0], xy[m, 1], s=3, alpha=.7, color=TCOL.get(t, '#7f7f7f'), label=t, linewidths=0)
    axes[0].set_title(title + '\nby batch', loc='left'); axes[1].set_title(title + '\nby cell type', loc='left')
    for a in axes:
        a.set_xticks([]); a.set_yticks([]); a.spines[:].set_visible(False)

def legends(fig, axes_batch, axes_type, ncol):
    h, l = axes_batch.get_legend_handles_labels(); h2, l2 = axes_type.get_legend_handles_labels()
    fig.legend(h + h2, l + l2, loc='lower center', ncol=ncol, markerscale=4, frameon=False,
               bbox_to_anchor=(0.5, -0.02), fontsize=15)

def fig_t02c():
    v = np.load(REPO / 'reports/T02c/umap_views.npz', allow_pickle=True)
    g = __import__('json').loads((REPO / 'reports/T02c/report.json').read_text())['gated']
    batch, types = v['dataset'].astype(str), v['cell_type'].astype(str)
    fig, ax = plt.subplots(2, 2, figsize=(8.1, 8.6))
    scatter_pair(ax[0], v['joint'], batch, types, (('A', 'v2 chemistry'), ('B', 'v3 chemistry')),
                 f"Before (mix {g['mix_joint']:.2f})")
    scatter_pair(ax[1], v['sep_load'], batch, types, (('A', 'v2 chemistry'), ('B', 'v3 chemistry')),
                 f"After rotation (mix {g['mix_sep_load']:.2f})")
    fig.tight_layout(rect=(0, 0.1, 1, 1)); legends(fig, ax[0, 0], ax[0, 1], 3); save(fig, 'fig2_T02c_umap')

def fig_o1b():
    r = pd.read_csv(REPO / 'reports/O1b/views.csv').set_index('view')
    a = sc.read_h5ad(DATA / 'pbmc_10k_v3/pbmc_10k_v3.h5ad'); a = a[a.obs['qc_pass'].to_numpy(dtype=bool)]
    tab = pd.read_csv(DATA / 'pbmc_10k_v3/k0_halves.csv').set_index('barcode')['half']
    names = a.obs_names.astype(str); half = np.asarray(names.map(tab).fillna('').astype(str))
    t1 = a.obs['cell_type'].astype(str).to_numpy()[half == 'H1']; t2 = a.obs['cell_type'].astype(str).to_numpy()[half == 'H2']
    z1 = np.load(REPO / 'reports/K0/latent_pca30_h1.npy').astype(np.float64)
    z2 = np.load(REPO / 'reports/K0/latent_pca30_h2.npy').astype(np.float64)
    assert len(z1) == len(t1) and len(z2) == len(t2)
    qh = np.load(REPO / 'reports/O1b/Q_hat_F.npy'); ql = np.load(REPO / 'reports/O1b/Q_load_F.npy')
    batch = np.asarray(['H1'] * len(z1) + ['H2'] * len(z2)); types = np.concatenate([t1, t2])
    views = [('separate_frames', np.vstack([z1, z2]), 'No rotation'),
             ('E5_estimate', np.vstack([z1 @ qh.T, z2]), 'scOPAlign (operators only)'),
             ('known_answer', np.vstack([z1 @ ql.T, z2]), 'Known answer')]
    fig, ax = plt.subplots(3, 2, figsize=(8.1, 12.6))
    for i, (key, rep, title) in enumerate(views):
        x = ad.AnnData(X=np.zeros((len(rep), 1), dtype=np.float32)); x.obsm['view'] = rep
        sc.pp.neighbors(x, n_neighbors=15, use_rep='view', random_state=20261014)
        sc.tl.umap(x, min_dist=0.5, random_state=20261014)
        scatter_pair(ax[i], np.asarray(x.obsm['X_umap']), batch, types, (('H1', 'half 1'), ('H2', 'half 2')),
                     f"{title}\nmix {r.loc[key, 'mix']:.2f} · frame {r.loc[key, 'frame_ratio']:.2f}")
        for a_ in ax[i]: a_.set_title(a_.get_title(loc='left').replace('\nby', ' · by'), loc='left', fontsize=15)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); legends(fig, ax[0, 0], ax[0, 1], 3); save(fig, 'fig3_O1b_umap')

def fig_k2():
    o = {d: np.load(REPO / f'reports/K2/ops_{d}.npz', allow_pickle=True) for d in ('a', 'b')}
    ta, tb = o['a']['cell_type'].astype(str), o['b']['cell_type'].astype(str)
    types = sorted(t for t in set(ta) & set(tb) if (ta == t).sum() >= 50 and (tb == t).sum() >= 50)
    fig, ax = plt.subplots(1, len(types), figsize=(8.1, 2.9), sharey=False)
    for a_, t in zip(np.atleast_1d(ax), types):
        for d, lab, col, tt in (('a', 'v2', BATCH[0], ta), ('b', 'v3', BATCH[1], tb)):
            a_.hist(o[d]['lam_seeds'][0][tt == t].ravel(), bins=36, density=True, histtype='step', lw=2, color=col, label=lab)
        a_.set_title(t, fontsize=15); a_.set_yticks([]); a_.tick_params(labelsize=12)
    np.atleast_1d(ax)[0].set_ylabel('density'); h, l = np.atleast_1d(ax)[0].get_legend_handles_labels(); fig.legend(h, l, loc='upper right', ncol=2, fontsize=13, frameon=False, bbox_to_anchor=(1, 1.08))
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.supxlabel('operator eigenvalue', fontsize=15, y=0.02); save(fig, 'fig1_K2_spectra')

if __name__ == '__main__':
    for f in sys.argv[1:] or ('k2', 't02c', 'o1b'):
        globals()[f'fig_{f}']()
        print('done', f, flush=True)
