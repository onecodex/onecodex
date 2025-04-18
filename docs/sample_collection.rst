****************
SampleCollection
****************

Data export, analysis and visualization functions are all contained within the
``SampleCollection`` class.

A ``SampleCollection`` is returned whenever multiple Samples are returned via
the One Codex API using a model.

Usage
=====

``SampleCollection`` contains useful tools for data export, analysis,
visualization and statistics. See the following sections for more information:

- :doc:`Visualization <visualization>`
- :doc:`Statistics <statistics>`
- :doc:`Taxonomy <taxonomy>`

.. code-block:: python

   import onecodex

   ocx = onecodex.Api()

   project = ocx.Project.get("d53ad03b010542e3")
   samples = ocx.Samples.where(project=project)

   type(samples) # SampleCollection


A ``SampleCollection`` can also be created manually from a list of samples:

.. code-block:: python

   import onecodex.models.collection.SampleCollection

   sample_list = [
       ocx.Samples.get("cee3b512605a43c6"),
       ocx.Samples.get("01f703ac505e4a30")
   ]

   samples = SampleCollection(sample_list)

   # convert classification results to a Pandas DataFrame
   samples.to_df()


``filter``
==========

.. automethod:: onecodex.models.collection.SampleCollection.filter

``to_otu``
==========

.. automethod:: onecodex.models.collection.SampleCollection.to_otu

``to_df``
=========

.. automethod:: onecodex.models.collection.SampleCollection.to_df

``to_classification_df``
=========================

.. automethod:: onecodex.models.collection.SampleCollection.to_classification_df

``to_functional_df``
=====================

.. automethod:: onecodex.models.collection.SampleCollection.to_functional_df
