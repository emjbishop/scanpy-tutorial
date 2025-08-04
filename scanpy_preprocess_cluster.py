import marimo

__generated_with = "0.14.16"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Imports and setup""")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    # Core scverse libraries
    import scanpy as sc
    import anndata as ad

    # Data retrieval
    import pooch
    return ad, pooch, sc


@app.cell
def _(sc):
    sc.settings.set_figure_params(dpi=50, facecolor="white")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Get example data""")
    return


@app.cell
def _(pooch):
    EXAMPLE_DATA = pooch.create(
        path=pooch.os_cache("scverse_tutorials"),
        base_url="doi:10.6084/m9.figshare.22716739.v1/",
    )
    EXAMPLE_DATA.load_registry_from_doi()
    return (EXAMPLE_DATA,)


@app.cell
def _(EXAMPLE_DATA, ad, sc):
    samples = {
        "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
        "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
    }
    adatas = {}

    for sample_id, filename in samples.items():
        path = EXAMPLE_DATA.fetch(filename)
        sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata

    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    print(adata.obs["sample"].value_counts())
    adata
    return (adata,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## QC Basics

    Interesting that ribosomal and hemoglobin genes are considered (the Seurat tutorial just focuses on mitochondrial).
    """
    )
    return


@app.cell
def _(adata):
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""Interesting that the data is log-transformed at this step (`log1p=True`), which is not done at this stage of the Seurat QC visualization. It doesn't seem like the log-transformed values are what's being visualized below, however. I wonder where the transformed values go.""")
    return


@app.cell
def _(adata, sc):
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""I find this plot a little confusing. I guess you'd want to see mostly dark purple dots located in the upper left (too many counts per cell hints at doublets). This seems like a plot you could make for every dataset and after a while you'd get a sense of what it looks like for a high quality dataset.""")
    return


@app.cell
def _(adata, sc):
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    Intersting suggested cutoffs:
    ```
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    ```

    These are the Seurat tutorial ones:
    ```
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    ```

    I will use the Seurat minimums since those have worked well for me in the past. I also like the upper-limit for number of genes since this acts as a first pass at removing doublets.
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pp.filter_cells(adata, min_genes=200)  # Using Seurat suggested minimum
    sc.pp.filter_genes(adata, min_cells=3)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Removing doublets

    I like scDblFinder when working with Seurat. I've heard of scrublet.

    I would normally remove doublets before clustering but I'll leave them in for now since that's what the tutorial suggests. Seems that to remove I would filter by `predicted_doublet` or `doublet_score` in `.obs` (part of the AnnData object).
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pp.scrublet(adata, batch_key="sample")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Normalization

    No SCTransform here. Seems like there is a desire to have this added to scanpy (https://github.com/scverse/scanpy/issues/1643). For now I'll do the normalization suggested in the tutorial.
    """
    )
    return


@app.cell
def _(adata):
    # Saving count data
    adata.layers["counts"] = adata.X.copy()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""I get a warning saying the data has already been transformed, which is super nice! In Seurat you could normalize over and over and it wouldn't tell you. This is probably from the `sc.pp.calculate_qc_metrics()` call above.""")
    return


@app.cell
def _(adata, sc):
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Feature selection

    Nice that the Suerat 'flavor' is the default (for reproducibility).
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""An interesting plot, haven't seen one like this.""")
    return


@app.cell
def _(adata, sc):
    sc.pl.highly_variable_genes(adata)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Dimension reduction

    Pretty standard. I think the Seurat tutorial does a better job explaining what these steps do.
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.tl.pca(adata)
    return


@app.cell
def _(adata, sc):
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    return


@app.cell
def _(adata, sc):
    sc.pl.pca(
        adata,
        color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ### Making the UMAPs

    Very similar to Seurat.
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color="sample",
        # Setting a smaller point size to get prevent overlap
        size=2,
    )

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Clustering

    Also pretty standard. I like that they use leidan instead of louvain, which is considered better by some papers (although I find louvain better sometimes).
    """
    )
    return


@app.cell
def _(adata, sc):
    # Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)

    sc.pl.umap(adata, color=["leiden"])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Re-assess

    I love that this is an included step! I do this also to remove clusters of dying cells and sometimes tweak my QC values.
    """
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.umap(
        adata,
        color=["leiden", "predicted_doublet", "doublet_score"],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    )
    return


@app.cell
def _(adata, sc):
    sc.pl.umap(
        adata,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    ## Annotation

    I really like their bubble [plots](https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#marker-gene-set). 
    """
    )
    return


if __name__ == "__main__":
    app.run()
