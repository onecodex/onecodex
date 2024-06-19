*************
Visualization
*************

.. admonition:: See Also
   :class: note

   Visualization functions are implemented as part of ``SampleCollection``. For
   more information, see :doc:`sample_collection`.

``plot_bargraph``
=================

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=10)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_bargraph(return_chart=True)

.. automethod:: onecodex.models.collection.SampleCollection.plot_bargraph


``plot_distance``
=================

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=10)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_distance(return_chart=True)

.. automethod:: onecodex.models.collection.SampleCollection.plot_distance

``plot_functional_heatmap``
===========================

.. automethod:: onecodex.models.collection.SampleCollection.plot_functional_heatmap

``plot_heatmap``
================

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=10)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_heatmap(return_chart=True)

.. automethod:: onecodex.models.collection.SampleCollection.plot_heatmap

``plot_mds``
============

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=10)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_mds(return_chart=True, color="country")

.. automethod:: onecodex.models.collection.SampleCollection.plot_mds

``plot_metadata``
=================

A general plotting tool which can be used to plot boxplots and scatter plots of
individual abundances or alpha-diversity metrics.

Alpha Diversity
---------------

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=20)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_metadata(return_chart=True, haxis="country")


2D Abundance Scatterplot
------------------------

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=20)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_metadata(return_chart=True, haxis="Bacteroides", vaxis="Firmicutes")

Boxplot
-------

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=20)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_metadata(return_chart=True, vaxis="Bacteroides", haxis="country")

.. automethod:: onecodex.models.collection.SampleCollection.plot_metadata

``plot_pca``
============

.. altair-plot::
    :strict:

    import onecodex

    ocx = onecodex.Api()
    project = ocx.Projects.get("d53ad03b010542e3")
    samples = ocx.Samples.where(project=project, public=True, limit=10)

    # note: return_chart is not needed if using a Jupyter notebook
    samples.plot_pca(return_chart=True, color="country")

.. automethod:: onecodex.models.collection.SampleCollection.plot_pca
