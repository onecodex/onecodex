from onecodex.lib.enums import Rank
from onecodex.viz._primitives import interleave_palette, prepare_props, get_base_classification_url
from onecodex.exceptions import OneCodexException, PlottingException
from onecodex.utils import is_continuous, has_missing_values


class VizPCAMixin(object):
    def plot_pca(
        self,
        rank=Rank.Auto,
        normalize="auto",
        org_vectors=0,
        org_vectors_scale=None,
        title=None,
        xlabel=None,
        ylabel=None,
        color=None,
        size=None,
        tooltip=None,
        return_chart=False,
        label=None,
        mark_size=100,
        width=None,
        height=None,
    ):
        """Perform principal component analysis and plot first two axes.

        Parameters
        ----------
        rank : {'auto', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'}, optional
            Analysis will be restricted to abundances of taxa at the specified level.
        normalize : 'auto' or `bool`, optional
            Convert read counts to relative abundances such that each sample sums to 1.0. Setting
            'auto' will choose automatically based on the data.
        org_vectors : `int`, optional
            Plot this many of the top-contributing eigenvectors from the PCA results.
        org_vectors_scale : `float`, optional
            Multiply the length of the lines representing the eigenvectors by this constant.
        title : `string`, optional
            Text label at the top of the plot.
        xlabel : `string`, optional
            Text label along the horizontal axis.
        ylabel : `string`, optional
            Text label along the vertical axis.
        size : `string` or `tuple`, optional
            A string or a tuple containing strings representing metadata fields. The size of points
            in the resulting plot will change based on the metadata associated with each sample.
        color : `string` or `tuple`, optional
            A string or a tuple containing strings representing metadata fields. The color of points
            in the resulting plot will change based on the metadata associated with each sample.
        tooltip : `string` or `list`, optional
            A string or list containing strings representing metadata fields. When a point in the
            plot is hovered over, the value of the metadata associated with that sample will be
            displayed in a modal.
        label : `string` or `callable`, optional
            A metadata field (or function) used to label each analysis. If passing a function, a
            dict containing the metadata for each analysis is passed as the first and only
            positional argument. The callable function must return a string.
        mark_size: `int`, optional
            The size of the points in the scatter plot.

        Examples
        --------
        Perform PCA on relative abundances at the species-level and color the resulting points by
        'geo_loc_name', a metadata field representing the geographical origin of each sample.

        >>> plot_pca(rank='species', normalize=True, color='geo_loc_name')

        Change the size of each point in the plot based on the abundance of Bacteroides.

        >>> plot_pca(size='Bacteroides')

        Display the abundances of Bacteroides, Prevotella, and Bifidobacterium in each sample when
        hovering over points in the plot.

        >>> plot_pca(tooltip=['Bacteroides', 'Prevotella', 'Bifidobacterium'])
        """
        # Deferred imports
        import altair as alt
        import numpy as np
        import pandas as pd
        from sklearn.decomposition import PCA

        if rank is None:
            raise OneCodexException("Please specify a rank or 'auto' to choose automatically")

        if len(self._results) < 3:
            raise PlottingException(
                "There are too few samples for PCA after filtering. Please select 3 or more "
                "samples to plot."
            )

        df = self.to_df(rank=rank, normalize=normalize)

        if len(df.columns) < 2:
            raise PlottingException(
                "There are too few taxa for PCA after filtering. Please select a rank that "
                "includes at least 2 taxa."
            )

        if tooltip:
            if not isinstance(tooltip, list):
                tooltip = [tooltip]
        else:
            tooltip = []

        tooltip.insert(0, "Label")

        if color and color not in tooltip:
            tooltip.insert(1, color)

        if size and size not in tooltip:
            tooltip.insert(2, size)

        magic_metadata, magic_fields = self._metadata_fetch(tooltip, label=label)

        pca = PCA()
        pca_vals = pca.fit(df.values).transform(df.values)
        pca_vals = pd.DataFrame(pca_vals, index=df.index)
        pca_vals.rename(columns=lambda x: "PC{}".format(x + 1), inplace=True)

        # label the axes
        if xlabel is None:
            xlabel = "PC1 ({}%)".format(round(pca.explained_variance_ratio_[0] * 100, 2))
        if ylabel is None:
            ylabel = "PC2 ({}%)".format(round(pca.explained_variance_ratio_[1] * 100, 2))

        # don't send all the data to vega, just what we're plotting
        plot_data = pd.concat(
            [pca_vals.loc[:, ("PC1", "PC2")], magic_metadata], axis=1
        ).reset_index()

        alt_kwargs = dict(
            x=alt.X("PC1", axis=alt.Axis(title=xlabel)),
            y=alt.Y("PC2", axis=alt.Axis(title=ylabel)),
            tooltip=[magic_fields[t] for t in tooltip],
            href="url:N",
            url=get_base_classification_url() + alt.datum.classification_id,
        )

        # only add these parameters if they are in use
        if color:
            color_kwargs = {
                "legend": alt.Legend(title=magic_fields[color]),
            }
            if not is_continuous(plot_data[color]) or has_missing_values(plot_data[color]):
                plot_data[color] = plot_data[color].fillna("N/A").astype(str)
                domain = plot_data[color].values
                color_range = interleave_palette(domain)
                color_kwargs["scale"] = alt.Scale(domain=domain, range=color_range)

            alt_kwargs["color"] = alt.Color(magic_fields[color], **color_kwargs)
        if size:
            alt_kwargs["size"] = magic_fields[size]

        chart = (
            alt.Chart(plot_data)
            .transform_calculate(url=alt_kwargs.pop("url"))
            .mark_circle(size=mark_size)
        )

        vector_chart = None
        # plot the organism eigenvectors that contribute the most
        if org_vectors > 0:
            plot_data = {
                "x": [],
                "y": [],
                "o": [],  # order these points should be connected in
                "Eigenvectors": [],
            }

            magnitudes = np.sqrt(pca.components_[0] ** 2 + pca.components_[1] ** 2)
            magnitudes.sort()
            cutoff = magnitudes[-1 * org_vectors]

            if org_vectors_scale is None:
                org_vectors_scale = 0.8 * np.max(pca_vals.abs().values)

            for tax_id, var1, var2 in zip(
                df.columns.values, pca.components_[0, :], pca.components_[1, :]
            ):
                if np.sqrt(var1**2 + var2**2) >= cutoff:
                    plot_data["x"].extend([0, var1 * float(org_vectors_scale)])
                    plot_data["y"].extend([0, var2 * float(org_vectors_scale)])
                    plot_data["o"].extend([0, 1])
                    plot_data["Eigenvectors"].extend([self.taxonomy["name"][tax_id]] * 2)

                    org_vectors -= 1

                    if org_vectors == 0:
                        break

            plot_data = pd.DataFrame(plot_data)

            vector_chart = (
                alt.Chart(plot_data)
                .mark_line(point=False)
                .encode(
                    x=alt.X("x", axis=None),
                    y=alt.Y("y", axis=None),
                    order="o",
                    color="Eigenvectors",
                )
            )

        chart = chart.encode(**alt_kwargs)

        if vector_chart:
            chart = alt.layer(chart, vector_chart).resolve_scale(color="independent")

        chart = chart.properties(**prepare_props(title=title, height=height, width=width))

        if return_chart:
            return chart
        else:
            chart.interactive().display()
