********
Querying
********

Every model (see :doc:`/models`) exposes a ``.where()`` classmethod for fetching
records from the One Codex server. Each keyword argument is a field to filter on,
and its value is either a literal to match exactly or a Mongo-style operator
``dict``::

    ocx.Samples.where(filename="SRR1234.fastq.gz")
    ocx.Samples.where(filename={"$icontains": "patient"})

Filters combine with a logical AND, and filtering happens server-side, so a
single ``.where()`` call is much faster than fetching records and looping::

    ocx.Samples.where(visibility="private", size={"$gte": 1_000_000})

Which operators a field accepts depends on the *type* of that field, not on the
model it belongs to — ``$icontains`` behaves the same on ``Samples.filename`` as
it does on ``Documents.filename``. Each model's ``.where()`` signature declares
the type of every field it accepts, so ``filename: str | StrFilter`` means
``filename`` takes either a plain string or a ``StrFilter`` operator dict. The
sections below describe each of those operator shapes; they are defined in
``onecodex.models.filters`` and may be imported if you want to annotate a
filter you are building up::

    from onecodex.models.filters import NumFilter

    big: NumFilter = {"$gte": 1_000_000}
    ocx.Samples.where(size=big)

All operator keys are optional.


Strings
=======

``StrFilter`` — e.g. ``Samples.filename``, ``Samples.status``,
``Projects.name``, ``Documents.filename``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Operator
     - Matches
   * - ``$eq``, ``$ne``
     - Exact (in)equality. ``None`` becomes an ``IS NULL`` / ``IS NOT NULL`` test.
   * - ``$in``
     - Any value in the given sequence.
   * - ``$contains``, ``$icontains``
     - Substring, case-sensitive and case-insensitive.
   * - ``$startswith``, ``$istartswith``
     - Prefix, case-sensitive and case-insensitive.
   * - ``$endswith``, ``$iendswith``
     - Suffix, case-sensitive and case-insensitive.

::

    ocx.Samples.where(filename={"$iendswith": ".fastq.gz"})
    ocx.Samples.where(status={"$in": ["ready", "processing"]})
    ocx.Samples.where(error_msg=None)


Enum-backed strings
===================

``EnumStrFilter`` — e.g. ``Metadata.platform``, ``Metadata.library_type``,
``Metadata.sample_type``.

These are string fields backed by an enum column on the server. They support
``$eq``, ``$ne``, ``$in``, ``$contains`` and ``$icontains``.

.. warning::

   The prefix and suffix matchers (``$startswith``, ``$endswith``, and their
   case-insensitive variants) are **not** available on these fields — they fail
   server-side rather than returning an empty result.

::

    ocx.Samples.where(platform="Illumina NovaSeq 6000")
    ocx.Samples.where(platform={"$icontains": "novaseq"})


Equality-only strings
=====================

``EqStrFilter`` — currently ``Users.email``.

Some fields are restricted server-side to equality tests only, accepting just
``$eq`` and ``$ne``. Substring and membership matching are unavailable::

    ocx.Users.where(email="you@example.com")


Numbers
=======

``NumFilter`` — e.g. ``Samples.size``, ``Documents.size``,
``Metadata.location_lat``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Operator
     - Matches
   * - ``$eq``, ``$ne``
     - Exact (in)equality. ``None`` tests for ``NULL``.
   * - ``$lt``, ``$lte``, ``$gt``, ``$gte``
     - Ordered comparison.
   * - ``$between``
     - A ``[low, high]`` range.
   * - ``$in``
     - Any value in the given sequence.

::

    ocx.Samples.where(size={"$gte": 1_000_000})
    ocx.Samples.where(size={"$between": [1_000_000, 5_000_000]})


Datetimes
=========

``DatetimeFilter`` — e.g. ``Samples.created_at``, ``Samples.updated_at``,
``Metadata.date_collected``.

The same operators as ``NumFilter``: ``$eq``, ``$ne``, ``$lt``, ``$lte``,
``$gt``, ``$gte``, ``$between`` and ``$in``.

.. warning::

   Datetime values must be **RFC3339 strings with a timezone offset**, not
   ``datetime.datetime`` objects — the latter are not JSON-serializable and will
   raise when the query is encoded. Call ``.isoformat()`` at the call site.

::

    from datetime import datetime, timedelta, timezone

    since = (datetime.now(timezone.utc) - timedelta(days=30)).isoformat()
    ocx.Samples.where(created_at={"$gte": since})


Booleans
========

``BoolFilter`` — e.g. ``Classifications.complete``, ``Classifications.success``,
``Metadata.starred``.

Accepts ``$eq`` and ``$ne`` only. Passing a bare boolean is usually clearer::

    ocx.Classifications.where(complete=True, success=True)


References to other models
==========================

``RefFilter`` — e.g. ``Samples.project``, ``Samples.owner``,
``Classifications.sample``.

Reference fields accept ``$eq``, ``$ne`` and ``$in``. Values may be a model
instance or an id string, used interchangeably::

    project = ocx.Projects.get("d53ad03b010542e3")

    ocx.Samples.where(project=project)
    ocx.Samples.where(project="d53ad03b010542e3")
    ocx.Samples.where(project={"$in": [project, "another_project_id"]})
    ocx.Samples.where(project=None)  # samples not in any project


Lists of references
===================

``ListRefFilter`` — fields holding a list of references, such as
``Samples.tags``.

Accepts ``$containsall`` (every given reference must be present) and
``$containsany`` (at least one must be present). For tags specifically,
``Samples.where()`` takes a plain list and builds the ``$containsall`` query for
you, resolving each entry by tag name, id, or :class:`Tags <onecodex.models.misc.Tags>`
instance::

    ocx.Samples.where(tags=["trimmed", "human-depleted"])


Fetching by ID
==============

Positional arguments to ``.where()`` are treated as record IDs, which is a
convenient way to fetch several known records in a single request::

    ocx.Samples.where("761bc54b97f64980", "0a1b2c3d4e5f6789")

Passing a single ``dict`` positionally sends it as a raw query, for cases the
keyword arguments do not cover::

    ocx.Samples.where({"visibility": "private", "filename": "tmp.fa"})

To fetch one record, use ``.get()``, which returns ``None`` when the record does
not exist or is not visible to you. To fetch every record, use ``.all()``.


Filtering client-side
=====================

Not every question can be answered by the server. The ``filter`` keyword takes a
callable applied to each fetched record, keeping those for which it returns a
truthy value::

    ocx.Samples.where(public=True, limit=100, filter=lambda s: not s.tags)

Because this runs after fetching, it narrows results but does not reduce the
amount of data transferred — prefer a server-side filter where one exists.
:doc:`SampleCollection </sample_collection>` results also offer a
:meth:`.filter() <onecodex.models.collection.SampleCollection.filter>` method for
the same purpose.

.. note::

   ``Metadata.custom`` is a free-form dictionary and is not server-filterable.
   Fetch a broader set of samples and filter on it client-side.
