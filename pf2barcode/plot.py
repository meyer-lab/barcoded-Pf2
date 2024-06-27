




def heatmap(
    adata: AnnData,
    var_names: _VarNames | Mapping[str, _VarNames],
    groupby: str | Sequence[str],
    *,
    use_raw: bool | None = None,
    log: bool = False,
    num_categories: int = 7,
    dendrogram: bool | str = False,
    gene_symbols: str | None = None,
    var_group_positions: Sequence[tuple[int, int]] | None = None,
    var_group_labels: Sequence[str] | None = None,
    var_group_rotation: float | None = None,
    layer: str | None = None,
    standard_scale: Literal["var", "obs"] | None = None,
    show_gene_labels: bool | None = None,
    show: bool | None = None,
    save: str | bool | None = None,
    figsize: tuple[float, float] | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    vcenter: float | None = None,
    norm: Normalize | None = None,
    **kwds,
) -> dict[str, Axes] | None:
    """\
    Heatmap of the expression values of genes.

    If `groupby` is given, the heatmap is ordered by the respective group. For
    example, a list of marker genes can be plotted, ordered by clustering. If
    the `groupby` observation annotation is not categorical the observation
    annotation is turned into a categorical by binning the data into the number
    specified in `num_categories`.

    Parameters
    ----------
    {common_plot_args}
    standard_scale
        Whether or not to standardize that dimension between 0 and 1, meaning for each variable or observation,
        subtract the minimum and divide each by its maximum.
    show_gene_labels
         By default gene labels are shown when there are 50 or less genes. Otherwise the labels are removed.
    {show_save_ax}
    {vminmax}
    **kwds
        Are passed to :func:`matplotlib.pyplot.imshow`.

    Returns
    -------
    Dict of :class:`~matplotlib.axes.Axes`

    Examples
    -------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        markers = ['C1QA', 'PSAP', 'CD79A', 'CD79B', 'CST3', 'LYZ']
        sc.pl.heatmap(adata, markers, groupby='bulk_labels', swap_axes=True)

    .. currentmodule:: scanpy

    See also
    --------
    pl.rank_genes_groups_heatmap
    tl.rank_genes_groups
    """
    var_names, var_group_labels, var_group_positions = _check_var_names_type(
        var_names, var_group_labels, var_group_positions
    )

    categories, obs_tidy = _prepare_dataframe(
        adata,
        var_names,
        groupby,
        use_raw=use_raw,
        log=log,
        num_categories=num_categories,
        gene_symbols=gene_symbols,
        layer=layer,
    )

    # check if var_group_labels are a subset of categories:
    if var_group_labels is not None:
        if set(var_group_labels).issubset(categories):
            var_groups_subset_of_groupby = True
        else:
            var_groups_subset_of_groupby = False

    if standard_scale == "obs":
        obs_tidy = obs_tidy.sub(obs_tidy.min(1), axis=0)
        obs_tidy = obs_tidy.div(obs_tidy.max(1), axis=0).fillna(0)
    elif standard_scale == "var":
        obs_tidy -= obs_tidy.min(0)
        obs_tidy = (obs_tidy / obs_tidy.max(0)).fillna(0)
    elif standard_scale is None:
        pass
    else:
        logg.warning("Unknown type for standard_scale, ignored")

    if groupby is None or len(categories) <= 1:
        categorical = False
        # dendrogram can only be computed  between groupby categories
        dendrogram = False
    else:
        categorical = True
        # get categories colors
        if isinstance(groupby, str) and isinstance(
            adata.obs[groupby].dtype, CategoricalDtype
        ):
            # saved category colors only work when groupby is valid adata.obs
            # categorical column. When groupby is a numerical column
            # or when groupby is a list of columns the colors are assigned on the fly,
            # which may create inconsistencies in multiple runs that require sorting
            # of the categories (eg. when dendrogram is plotted).
            if groupby + "_colors" not in adata.uns:
                # if colors are not found, assign a new palette
                # and save it using the same code for embeddings
                from ._tools.scatterplots import _get_palette

                _get_palette(adata, groupby)
            groupby_colors = adata.uns[groupby + "_colors"]
        else:
            # this case happen when adata.obs[groupby] is numeric
            # the values are converted into a category on the fly
            groupby_colors = None

    if dendrogram:
        dendro_data = _reorder_categories_after_dendrogram(
            adata,
            groupby,
            dendrogram,
            var_names=var_names,
            var_group_labels=var_group_labels,
            var_group_positions=var_group_positions,
            categories=categories,
        )

        var_group_labels = dendro_data["var_group_labels"]
        var_group_positions = dendro_data["var_group_positions"]

        # reorder obs_tidy
        if dendro_data["var_names_idx_ordered"] is not None:
            obs_tidy = obs_tidy.iloc[:, dendro_data["var_names_idx_ordered"]]
            var_names = [var_names[x] for x in dendro_data["var_names_idx_ordered"]]

        obs_tidy.index = obs_tidy.index.reorder_categories(
            [categories[x] for x in dendro_data["categories_idx_ordered"]],
            ordered=True,
        )

        # reorder groupby colors
        if groupby_colors is not None:
            groupby_colors = [
                groupby_colors[x] for x in dendro_data["categories_idx_ordered"]
            ]

    if show_gene_labels is None:
        if len(var_names) <= 50:
            show_gene_labels = True
        else:
            show_gene_labels = False
            logg.warning(
                "Gene labels are not shown when more than 50 genes are visualized. "
                "To show gene labels set `show_gene_labels=True`"
            )
    if categorical:
        obs_tidy = obs_tidy.sort_index()

    colorbar_width = 0.2
    norm = check_colornorm(vmin, vmax, vcenter, norm)

    # define a layout of 2 rows x 4 columns
    # first row is for 'brackets' (if no brackets needed, the height of this row
    # is zero) second row is for main content. This second row is divided into
    # three axes:
    #   first ax is for the categories defined by `groupby`
    #   second ax is for the heatmap
    #   third ax is for the dendrogram
    #   fourth ax is for colorbar

    dendro_width = 1 if dendrogram else 0
    groupby_width = 0.2 if categorical else 0
    if figsize is None:
        height = 6
        if show_gene_labels:
            heatmap_width = len(var_names) * 0.3
        else:
            heatmap_width = 8
        width = heatmap_width + dendro_width + groupby_width
    else:
        width, height = figsize
        heatmap_width = width - (dendro_width + groupby_width)

    if var_group_positions is not None and len(var_group_positions) > 0:
        # add some space in case 'brackets' want to be plotted on top of the image
        height_ratios = [0.15, height]
    else:
        height_ratios = [0, height]

    width_ratios = [
        groupby_width,
        heatmap_width,
        dendro_width,
        colorbar_width,
    ]
    fig = plt.figure(figsize=(width, height))

    axs = gridspec.GridSpec(
        nrows=2,
        ncols=4,
        width_ratios=width_ratios,
        wspace=0.15 / width,
        hspace=0.13 / height,
        height_ratios=height_ratios,
    )

    heatmap_ax = fig.add_subplot(axs[1, 1])
    kwds.setdefault("interpolation", "nearest")
    im = heatmap_ax.imshow(obs_tidy.values, aspect="auto", norm=norm, **kwds)

    heatmap_ax.set_ylim(obs_tidy.shape[0] - 0.5, -0.5)
    heatmap_ax.set_xlim(-0.5, obs_tidy.shape[1] - 0.5)
    heatmap_ax.tick_params(axis="y", left=False, labelleft=False)
    heatmap_ax.set_ylabel("")
    heatmap_ax.grid(False)

    if show_gene_labels:
        heatmap_ax.tick_params(axis="x", labelsize="small")
        heatmap_ax.set_xticks(np.arange(len(var_names)))
        heatmap_ax.set_xticklabels(var_names, rotation=90)
    else:
        heatmap_ax.tick_params(axis="x", labelbottom=False, bottom=False)
    # plot colorbar
    _plot_colorbar(im, fig, axs[1, 3])

    if categorical:
        groupby_ax = fig.add_subplot(axs[1, 0])
        (
            label2code,
            ticks,
            labels,
            groupby_cmap,
            norm,
        ) = _plot_categories_as_colorblocks(
            groupby_ax, obs_tidy, colors=groupby_colors, orientation="left"
        )

        # add lines to main heatmap
        line_positions = (
            np.cumsum(obs_tidy.index.value_counts(sort=False))[:-1] - 0.5
        )
        heatmap_ax.hlines(
            line_positions,
            -0.5,
            len(var_names) - 0.5,
            lw=1,
            color="black",
            zorder=10,
            clip_on=False,
        )

    if dendrogram:
        dendro_ax = fig.add_subplot(axs[1, 2], sharey=heatmap_ax)
        _plot_dendrogram(
            dendro_ax, adata, groupby, ticks=ticks, dendrogram_key=dendrogram
        )

    # plot group legends on top of heatmap_ax (if given)
    if var_group_positions is not None and len(var_group_positions) > 0:
        gene_groups_ax = fig.add_subplot(axs[0, 1], sharex=heatmap_ax)
        _plot_gene_groups_brackets(
            gene_groups_ax,
            group_positions=var_group_positions,
            group_labels=var_group_labels,
            rotation=var_group_rotation,
            left_adjustment=-0.3,
            right_adjustment=0.3,
        )

    return_ax_dict = {"heatmap_ax": heatmap_ax}
    if categorical:
        return_ax_dict["groupby_ax"] = groupby_ax
    if dendrogram:
        return_ax_dict["dendrogram_ax"] = dendro_ax
    if var_group_positions is not None and len(var_group_positions) > 0:
        return_ax_dict["gene_groups_ax"] = gene_groups_ax

    _utils.savefig_or_show("heatmap", show=show, save=save)
    show = settings.autoshow if show is None else show
    if show:
        return None
    return return_ax_dict