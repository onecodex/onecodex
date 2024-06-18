*******************************
One Codex Python Client Library
*******************************

Introduction
============

This package provides 3 major pieces of functionality: (1) a core Python client
library; (2) a simple CLI for interacting with the `One Codex platform
<https://onecodex.com>`_ that uses that core library; and (3) optional
extensions to the client library, which offers many features aimed at advanced
users and provides functionality for use in interactive notebook environments
(e.g., Jupyter notebooks).

For installation instructions and instructions on using the command-line
client, see the :doc:`readme`. You can find the source code on `GitHub
<https://github.com/onecodex/onecodex>`_.

.. toctree::
   :maxdepth: 5
   :caption: Contents
   :hidden:

   sample_collection
   Visualization <visualization>
   Statistics <statistics>
   Taxonomy <taxonomy>
   Models <models>
   notebooks

   GitHub <https://github.com/onecodex/onecodex>


Installation
============

.. code-block:: bash

   # The CLI (and core Python library)
   pip install onecodex

   # Optional dependencies (visualization, analysis and reports)
   pip install onecodex[all]


Quickstart
==========

.. altair-plot::
    :strict:

    import onecodex

    # Instantiate the API (run onecodex login first)
    ocx = onecodex.Api()

    # Fetch some samples
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=10)

    # Generate a Pandas DataFrame from classification results
    results = samples.to_df()

    # Plot a bar graph
    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_bargraph(return_chart=True)


Tutorials
=========

For documentation and examples of common uses of the client library including
interactive data analysis and plotting within Jupyter notebooks see
:doc:`notebooks`.

CLI
===

The Python package also comes with a command-line tool which can be used to
upload samples to One Codex as well as results from the API. See
:doc:`readme` for more information.


See Also
========

* `Source on GitHub <https://github.com/onecodex/onecodex>`_
* `Package on PyPI <https://pypi.org/project/onecodex/>`_
* `One Codex API Documentation <https://developer.onecodex.com>`_
* `One Codex <https://onecodex.com>`_
