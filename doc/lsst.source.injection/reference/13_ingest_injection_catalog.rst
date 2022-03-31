.. _lsst.source.injection-ref-ingest:

=============================
 Ingest an Injection Catalog
=============================

-------------------------------------------------------------
 Ingesting a Synthetic Source Catalog Into a Data Repository
-------------------------------------------------------------

Once a synthetic source catalog has been constructed, it can be ingested into a data repository for subsequent use.
This is the final required step before running source injection tasks.

The sections below detail how to ingest an injection catalog, on the command line using the :doc:`ingest_injection_catalog <../scripts/ingest_injection_catalog>` tool, or in Python using the :py:func:`~lsst.source.injection.ingest_injection_catalog` Python function.

.. _lsst.source.injection-ref-ingest-cli:

Ingest an Injection Catalog on the Command Line
===============================================

The :doc:`ingest_injection_catalog <../scripts/ingest_injection_catalog>` tool is used to ingest an injection catalog into a data repository from the command line.
For example, to ingest the catalog ``my_injection_catalog.csv`` into the ``u/$USER/my_injection_inputs`` RUN collection in the ``$REPO`` repository, associating this catalog with the ``g`` band, run:

.. code-block:: bash

    ingest_injection_catalog \
    -b $REPO \
    -i my_injection_catalog.csv g \
    -o u/$USER/my_injection_inputs

*where*

    `$REPO`
        The path to the butler repository.

    `$USER`
        The user's username.

If successful,  the tool will print a message similar to:

.. code-block:: bash

    Ingested 3 g band injection_catalog DatasetRefs into the butler.

In this example, the right ascension and declination values in the input catalog fall across three HTM7 trixels, resulting in three datasets being ingested.

If a single injection catalog is to be associated with multiple bands (i.e., no variation in injected source flux as a function of band), then multiple space-separated bands can be specified at the same time above for convenience, e.g.:

.. code-block:: bash

    -i my_injection_catalog.csv g r i z y

Any catalog format supported by the `Astropy Table <http://docs.astropy.org/en/stable/table/>`_ class can be ingested.
An attempt to auto-detect the catalog format will be made, but this can be overridden using the ``--format`` argument.

By default (and by convention), the output dataset type name for ingested data is ``injection_catalog``.
This can be overridden using the ``-t`` argument.

.. _lsst.source.injection-ref-ingest-python:

Ingest an Injection Catalog in Python
=====================================

The :py:func:`~lsst.source.injection.ingest_injection_catalog` Python function is used to generate a synthetic source injection catalog in Python:

.. code-block:: python

    from lsst.source.injection import ingest_injection_catalog

More information on the operation of this function may be obtained by calling ``ingest_injection_catalog?`` in a Python interpreter.

For example, the snippet below ingests the ``my_injection_catalog`` object into a writeable data butler, associating this catalog with the ``g`` band, and storing the resulting dataset in the ``u/$USER/my_injection_inputs`` RUN collection:

.. code-block:: python

    import os
    from lsst.daf.butler import Butler

    # Get username.
    user = os.getenv("USER")

    # Instantiate a writeable Butler.
    writeable_butler = Butler(REPO, writeable=True)

    # Ingest the injection catalog.
    my_injected_datasetRefs = ingest_injection_catalog(
        writeable_butler=writeable_butler,
        table=my_injection_catalog,
        band="g",
        output_collection=f"u/{user}/my_injection_inputs",
    )

*where*

    `REPO`
        The path to the butler repository.

.. caution::

    Be careful when utilizing a writeable Butler, as edits to the data repository can inadvertantly be made.

If successful, a list of dataset reference IDs will be returned, one per HTM7 trixel that the input catalog spans.
The output dataset type name will be ``injection_catalog`` by default (and convention), but this can be overridden by setting the ``dataset_type_name`` argument if so desired.

.. _lsst.source.injection-ref-ingest-wrap:

Wrap Up
=======

This page has described how to ingest a synthetic source injection catalog for use with the LSST Science Pipelines, both on the command line and in Python.
For source injections into multiple bands, the above commands may be called multiple times to associate different injection catalogs with different bands.

Move on to :ref:`another quick reference guide <lsst.source.injection-ref>`, consult the :ref:`FAQs <lsst.source.injection-faqs>`, or head back to the `main page <..>`_.
