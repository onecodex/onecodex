******
Models
******

Each model maps to a resource in the One Codex API:

- A :class:`Samples <onecodex.models.sample.Samples>` record is a sequencing
  file you have uploaded, along with its
  :class:`Metadata <onecodex.models.sample.Metadata>`.
- A :class:`Jobs <onecodex.models.misc.Jobs>` record is something that can be
  run against a sample -- either a built-in One Codex analysis or a custom
  workflow you have defined yourself.
- Running a job on a sample produces an
  :class:`Analyses <onecodex.models.analysis.Analyses>` record, which links back
  to both.
- Each type of job produces a specific type of analysis, with its own model and
  its own fields: a classification job produces a
  :class:`Classifications <onecodex.models.analysis.Classifications>`, an
  alignment job an :class:`Alignments <onecodex.models.analysis.Alignments>`,
  and likewise for
  :class:`FunctionalProfiles <onecodex.models.analysis.FunctionalProfiles>`,
  :class:`Panels <onecodex.models.analysis.Panels>` and
  :class:`Mlsts <onecodex.models.analysis.Mlsts>`. Query
  :class:`Analyses <onecodex.models.analysis.Analyses>` directly, filtering on
  ``analysis_type``, when you want a sample's analyses regardless of type.
- A custom, user-defined job produces a
  :class:`Workflows <onecodex.models.analysis.Workflows>` record. Reference data
  those jobs depend on is uploaded as
  :class:`Assets <onecodex.models.misc.Assets>`.

Every model inherits from
:class:`OneCodexBase <onecodex.models.base.OneCodexBase>`, which provides the
shared ``.get()``, ``.all()`` and ``.where()`` methods. See :doc:`/querying` for
the operators ``.where()`` accepts.

Core Models
===========

``OneCodexBase``
----------------

.. autoclass:: onecodex.models.base.OneCodexBase
   :members:
   :show-inheritance:

Note: this is a base class that is not meant to be instantiated directly. See
this class's documentation for details about common methods such as `.where`,
`.all` and `.get` which are used by the models listed below.

``Samples``
-----------

.. autoclass:: onecodex.models.sample.Samples
   :members:
   :show-inheritance:

``Metadata``
------------

.. autoclass:: onecodex.models.sample.Metadata
   :members:
   :show-inheritance:

``Projects``
------------

.. autoclass:: onecodex.models.misc.Projects
   :members:
   :show-inheritance:

``Tags``
--------

.. autoclass:: onecodex.models.misc.Tags
   :members:
   :show-inheritance:

``Users``
---------

.. autoclass:: onecodex.models.misc.Users
   :members:
   :show-inheritance:

``Documents``
-------------

.. autoclass:: onecodex.models.misc.Documents
   :members:
   :show-inheritance:

``Jobs``
--------

.. autoclass:: onecodex.models.misc.Jobs
   :members:
   :show-inheritance:

Analysis Results
================

``Analyses``
------------

.. autoclass:: onecodex.models.analysis.Analyses
   :members:
   :show-inheritance:

``Alignments``
--------------

.. autoclass:: onecodex.models.analysis.Alignments
   :members:
   :show-inheritance:

``Classifications``
-------------------

.. autoclass:: onecodex.models.analysis.Classifications
   :members:
   :show-inheritance:

``FunctionalProfiles``
----------------------

.. autoclass:: onecodex.models.analysis.FunctionalProfiles
   :members:
   :show-inheritance:

``Panels``
----------

.. autoclass:: onecodex.models.analysis.Panels
   :members:
   :show-inheritance:

``Mlsts``
---------

.. autoclass:: onecodex.models.analysis.Mlsts
   :members:
   :show-inheritance:

Custom Workflows
================

``Workflows``
-------------

.. autoclass:: onecodex.models.analysis.Workflows
   :members:
   :show-inheritance:

``Assets``
----------

.. autoclass:: onecodex.models.misc.Assets
   :members:
   :show-inheritance:
