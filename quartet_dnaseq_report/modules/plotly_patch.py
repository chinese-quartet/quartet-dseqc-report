# -*- coding: utf-8 -*-

from __future__ import absolute_import

from collections import OrderedDict

from plotly import exceptions, optional_imports
from plotly.graph_objs import graph_objs
import plotly.figure_factory._dendrogram as dd

# Optional imports, may be None for users that only use our core functionality.
np = optional_imports.get_module("numpy")
scp = optional_imports.get_module("scipy")
sch = optional_imports.get_module("scipy.cluster.hierarchy")
scs = optional_imports.get_module("scipy.spatial")


def create_dendrogram(
    X,
    orientation="bottom",
    labels=None,
    colorscale=None,
    distfun=None,
    linkagefun=lambda x: sch.linkage(x, "complete"),
    hovertext=None,
    color_threshold=None,
):
    """
    BETA function that returns a dendrogram Plotly figure object.

    :param (ndarray) X: Matrix of observations as array of arrays
    :param (str) orientation: 'top', 'right', 'bottom', or 'left'
    :param (list) labels: List of axis category labels(observation labels)
    :param (list) colorscale: Optional colorscale for dendrogram tree
    :param (function) distfun: Function to compute the pairwise distance from
                               the observations
    :param (function) linkagefun: Function to compute the linkage matrix from
                               the pairwise distances
    :param (list[list]) hovertext: List of hovertext for constituent traces of dendrogram
                               clusters
    :param (double) color_threshold: Value at which the separation of clusters will be made

    Example 1: Simple bottom oriented dendrogram

    >>> from plotly.figure_factory import create_dendrogram

    >>> import numpy as np

    >>> X = np.random.rand(10,10)
    >>> fig = create_dendrogram(X)
    >>> fig.show()

    Example 2: Dendrogram to put on the left of the heatmap

    >>> from plotly.figure_factory import create_dendrogram

    >>> import numpy as np

    >>> X = np.random.rand(5,5)
    >>> names = ['Jack', 'Oxana', 'John', 'Chelsea', 'Mark']
    >>> dendro = create_dendrogram(X, orientation='right', labels=names)
    >>> dendro['layout'].update({'width':700, 'height':500})
    >>> dendro.show()

    Example 3: Dendrogram with Pandas

    >>> from plotly.figure_factory import create_dendrogram

    >>> import numpy as np
    >>> import pandas as pd

    >>> Index= ['A','B','C','D','E','F','G','H','I','J']
    >>> df = pd.DataFrame(abs(np.random.randn(10, 10)), index=Index)
    >>> fig = create_dendrogram(df, labels=Index)
    >>> fig.show()
    """
    if not scp or not scs or not sch:
        raise ImportError(
            "FigureFactory.create_dendrogram requires scipy, \
                            scipy.spatial and scipy.hierarchy"
        )

    s = X.shape
    if len(s) != 2:
        exceptions.PlotlyError("X should be 2-dimensional array.")

    if distfun is None:
        distfun = scs.distance.pdist

    dendrogram = Dendrogram(
        X,
        orientation,
        labels,
        colorscale,
        distfun=distfun,
        linkagefun=linkagefun,
        hovertext=hovertext,
        color_threshold=color_threshold,
    )

    return graph_objs.Figure(data=dendrogram.data, layout=dendrogram.layout)


class Dendrogram(dd._Dendrogram):
    """Refer to FigureFactory.create_dendrogram() for docstring."""

    def __init__(
        self,
        X,
        orientation="bottom",
        labels=None,
        colorscale=None,
        width=np.inf,
        height=np.inf,
        xaxis="xaxis",
        yaxis="yaxis",
        distfun=None,
        linkagefun=lambda x: sch.linkage(x, "complete"),
        hovertext=None,
        color_threshold=None,
    ):
        super(Dendrogram, self).__init__(
            X,
            orientation=orientation,
            labels=labels,
            colorscale=colorscale,
            width=width,
            height=height,
            xaxis=xaxis,
            yaxis=yaxis,
            distfun=distfun,
            linkagefun=linkagefun,
            hovertext=hovertext,
            color_threshold=color_threshold,
        )

    def get_color_dict(self, colorscale):
        """
        Returns colorscale used for dendrogram tree clusters.

        :param (list) colorscale: Colors to use for the plot in rgb format.
        :rtype (dict): A dict of default colors mapped to the user colorscale.

        """

        # These are the color codes returned for dendrograms
        # We're replacing them with nicer colors
        d = {
            "C0": "red",
            "C1": "green",
            "C2": "blue",
            "C3": "cyan",
            "C4": "magenta",
            "C5": "yellow",
            "C6": "black",
            "C7": "white",
        }

        default_colors = OrderedDict(sorted(d.items(), key=lambda t: t[0]))

        if colorscale is None:
            colorscale = [
                "rgb(0,116,217)",  # blue
                "rgb(35,205,205)",  # cyan
                "rgb(61,153,112)",  # green
                "rgb(40,35,35)",  # black
                "rgb(133,20,75)",  # magenta
                "rgb(255,65,54)",  # red
                "rgb(255,255,255)",  # white
                "rgb(255,220,0)",
            ]  # yellow

        for i in range(len(default_colors.keys())):
            k = list(default_colors.keys())[i]  # PY3 won't index keys
            if i < len(colorscale):
                default_colors[k] = colorscale[i]

        return default_colors
