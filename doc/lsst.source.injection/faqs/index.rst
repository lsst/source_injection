.. _lsst.source.injection-faqs:

============================
 Frequently Asked Questions
============================

---------------------------------------------------------------------------
 Frequently Asked Questions (and Answers) about Synthetic Source Injection
---------------------------------------------------------------------------

This page contains a list of frequently asked questions (and answers) about the LSST Science Pipelines' synthetic source injection framework.

.. seealso::

    Full documentation for this repository can be found on the :ref:`main page <lsst.source.injection>`.

.. figure:: _assets/t9813p42i_zoom_stamp_prepost_injection.gif
    :name: t9813p42i_zoom_stamp_prepost_injection
    :alt: An HSC i-band cutout from tract 9813, patch 42, showcasing the injection of the Rubin Observatory logo.
    :align: center
    :width: 100%

    ..

    An HSC i-band cutout from tract 9813, patch 42, showcasing the injection of the Vera C. Rubin Observatory logo.
    Images are ~100 arcseconds on the short axis, log scaled across the central 99.5% flux range, and smoothed with a Gaussian kernel of FWHM 3 pixels.

    .. list-table::
        :widths: 1 1

        * - .. figure:: _assets/t9813p42i_zoom_stamp_pre_injection.png
                :name: t9813p42i_zoom_stamp_pre_injection
                :alt: Tract 9813, patch 42, HSC i-band cutout, before postage stamp injection.
                :align: center
                :width: 100%

                ..

                Before injection.
          - .. figure:: _assets/t9813p42i_zoom_stamp_post_injection.png
                :name: t9813p42i_zoom_stamp_post_injection
                :alt: Tract 9813, patch 42, HSC i-band cutout, after postage stamp injection.
                :align: center
                :width: 100%

                ..

                After injection.

.. _lsst.source.injection-faqs-source-injection-first:

Do I need to run source injection on its own first, or as part of a subset?
===========================================================================

That depends on what you're trying to do.
The instructions in :ref:`lsst.source.injection-ref-make` will add your injection task into any associated subsets.
For example, the ``inject_exposure`` task will be added to the ``step1`` subset for the HSC DRP-RC2 pipeline (as of ``w_2023_39``).
Therefore, instead of running ``inject_exposure`` directly, you can instead run ``step1`` to reduce raw data, inject sources, and generate calibrated outputs all in one go.

To check what your pipeline looks like and see which tasks get run and in what order, see the instructions in :ref:`lsst.source.injection-ref-make-visualize`.

.. _lsst.source.injection-faqs-missing-injection-task:

Why don't I see the injection task in my output pipeline YAML?
==============================================================

It may be difficult to see which tasks are being run in your pipeline by looking at the output YAML directly.
Instead, use the instructions in :ref:`lsst.source.injection-ref-make-visualize` to visualize your pipeline and see which tasks are being run and in what order.

.. _lsst.source.injection-faqs-no-injection-overlap:

Why do some of my detectors raise the error "No injection sources overlap the data query"?
=============================================================================================

This error is raised when your input injection catalog does not cover either the spatial region or the band (or both) specified by the input data query.
For example, at the visit-level, your input injection catalog may only cover a subset of the CCDs in the visit.
Alternatively, your input injection catalog may have been registered in a band which has not been specified in the data query.

If the above error is raised due to a lack of spatial overlap, to override this error you may set the config option ``process_all_data_ids`` to ``True``.
This will cause the injection task to instead copy any non-overlapping data IDs into your output collection, prefixing the dataset type name with ``injected_`` and appending the ``INJECTED`` mask plane.
This can be useful if, for example, you require an entire visit with a common dataset type name for downstream data processing tasks (such as full focal plane sky correction).

All configuration options for a pipeline, a subset or a given task can be shown on the command line using ``pipetask build``, e.g.:

.. code-block:: shell

    pipetask build \
    -p DRP-RC2+injection.yaml#inject_exposure \
    --show config

.. _lsst.source.injection-faqs-not-enough-datasets:

Why am I seeing FileNotFoundError: Not enough datasets (0) found for non-optional connection?
=================================================================================================

This is usually an indicator that you're not providing enough input data to the pipeline you're trying to run.
Check that all required inputs are available in your input collection, e.g., by using the instructions in :ref:`lsst.source.injection-ref-make-visualize`.

.. _lsst.source.injection-faqs-spatial-overlap-1:

How can I determine which visits/tracts/patches overlap my source injection catalogs?
=====================================================================================

There are a number of ways to do this with the LSST Science Pipelines.
Much of this can be achieved on the command line, however, the examples below will focus on working in a Python shell to leverage the greater flexibility in that environment.

.. note::

    This method looks for visit/tract/patch overlap with the HTM7 trixel IDs that injection catalogs are registered against.
    An HTM7 trixel may span a larger area than your input RA and Dec coordinates, so this method may return more results than you expect.

    For an alternative approach that uses injection RA/Dec coordinates directly, see below.

First, lets extract the HTM7 trixel IDs for all `injection_catalog` datasets in your input collection:

.. code-block:: python

    import os
    from lsst.daf.butler import Butler

    # Instantiate a butler and a user.
    butler = Butler("/path/to/my/butler/repo")
    user = os.getenv("USER")

    # Get all injection catalog dataset references.
    injection_refs = butler.registry.queryDatasets(
        "injection_catalog",
        collections=f"u/{user}/my_injection_inputs",
    )

    # Extract the HTM7 trixel IDs.
    injection_trixels = tuple(x.dataId["htm7"] for x in injection_refs)

This will return a tuple of all HTM7 trixels that your injection catalog spans.

.. tip::

    You can add a ``where`` argument to the ``queryDatasets`` call above to filter the results to a specific data query only, e.g.:

    .. code-block:: python

        where="band='i'",

Now you have a tuple of injection trixel IDs, you can again query the butler registry for other data IDs that overlap these trixels.

For example, to find the visit numbers for all ``raw`` dataset types registered in the ``HSC/runs/RC2/w_2023_32/DM-40356`` collection that overlap our list of injection trixels in the i-band:

.. code-block:: python

    # Get the data IDs for all overlapping raw datasets.
    injection_raw_ids = set(
        butler.registry.queryDataIds(
            ["visit"],
            datasets="raw",
            where=f"htm7 IN {injection_trixels} AND band='i'",
            collections="HSC/runs/RC2/w_2023_32/DM-40356",
        )
    )

    # Extract the visit number for all overlapping raw datasets.
    injection_visits = sorted({x["visit"] for x in injection_raw_ids})

The process for finding overlapping tracts and patches is similar, instead querying the ``deepCoadd`` dataset type:

.. code-block:: python

    # Get the data IDs for all overlapping deepCoadd datasets.
    injection_deepCoadd_ids = set(
        butler.registry.queryDataIds(
            ["tract", "patch"],
            datasets="deepCoadd",
            where=f"htm7 IN {injection_trixels} AND band='i'",
            collections="HSC/runs/RC2/w_2023_32/DM-40356",
        )
    )

    # Extract the tract IDs for all overlapping deepCoadd datasets.
    injection_tracts = sorted({x["tract"] for x in injection_deepCoadd_ids})

    # Format the results into injection_tract_patch_dict.
    injection_tract_patch_dict = {}
    for injection_deepCoadd_id in injection_deepCoadd_ids:
        tract_id = injection_deepCoadd_id["tract"]
        patch_id = injection_deepCoadd_id["patch"]
        if tract_id in injection_tract_patch_dict:
            injection_tract_patch_dict[tract_id].append(patch_id)
        else:
            injection_tract_patch_dict[tract_id] = [patch_id]
    # Sort the patch ID list for each tract.
    for patch_list in injection_tract_patch_dict.values():
        patch_list.sort()

The output ``injection_tract_patch_dict`` will be a dictionary of overlapping tracts and patches for all considered HTM7 trixels.
Tract IDs are used as keys, with a list of patch IDs as values.

.. _lsst.source.injection-faqs-spatial-overlap-2:

How can I determine which tracts/patches overlap my source injection coordinates using a sky map?
=================================================================================================

Any registered sky map can be used to determine which tracts overlap a set of injection coordinates.

First, lets query the butler registry to get all i-band injection catalog dataset references:

.. code-block:: python

    import os
    from lsst.daf.butler import Butler

    # Instantiate a butler and a user.
    butler = Butler("/path/to/my/butler/repo")
    user = os.getenv("USER")

    # Get all injection catalog dataset references.
    injection_refs = butler.registry.queryDatasets(
        "injection_catalog",
        collections=f"u/{user}/my_injection_inputs",
        where="band='i'",
    )

Now we can use the butler to get all of our injection catalog data, using the astropy `astropy.table.vstack` function to concatenate them all together and converting the RA and Dec data for each source to `lsst.geom.SpherePoint` objects:

.. code-block:: python

    from astropy.table import vstack
    from lsst.geom import SpherePoint, degrees

    # Recombine all injection catalogs into a single table.
    injection_catalogs = [butler.get(x) for x in injection_refs]
    injection_catalog = vstack(injection_catalogs)
    injection_catalog.sort("injection_id")

    # Convert RA and Dec to SpherePoint objects.
    injection_sphere_points = [
        SpherePoint(ra, dec, units=degrees)
        for ra, dec in zip(injection_catalog["ra"], injection_catalog["dec"])
    ]

Finally, we can query a sky map such as "``hsc_rings_v1``" for all tract and patch overlaps:

.. code-block:: python

    # Get the sky map.
    skymap = butler.get(
        "skyMap",
        collections="skymaps",
        skymap="hsc_rings_v1",
    )

    # Find all tract and patch overlaps.
    injection_tract_patch_info = skymap.findTractPatchList(injection_sphere_points)

    # Format the results into injection_tract_patch_dict.
    injection_tract_patch_dict = {}
    for tract_info, patch_info in injection_tract_patch_info:
        tract_id = tract_info.tract_id
        patch_ids = [patch.sequential_index for patch in patch_info]
        injection_tract_patch_dict[tract_id] = sorted(patch_ids)

The output ``injection_tract_patch_dict`` will be a dictionary of overlapping tracts and patches for all considered injection coordinates.
Tract IDs are used as keys, with a list of patch IDs as values.

.. _lsst.source.injection-faqs-final-visit-summary:

Why am I seeing 'CRITICAL: No datasets of type finalVisitSummary in collection' when I try to run source injection?
=====================================================================================================================

This error is raised when you try to run visit-level source injection and your input collection does not contain a ``finalVisitSummary`` dataset.
The ``finalVisitSummary`` table is a visit-level data product from a prior data reduction run.
It contains our best-estimate final measurements for various visit-level data products, such as the PSF, photometric calibration and WCS solution.

Source injection into visit-level data typically requires a ``finalVisitSummary`` dataset to be present in the input collection, unless use of these external data have been explicitly disabled by the user (e.g., by setting ``external_psf=False`` in the case of the PSF).

.. caution::

    Where possible, it is strongly recommended to use a ``finalVisitSummary`` table for visit-level source injection.

    If use of external data from the ``finalVisitSummary`` table *has* been disabled by the user at runtime, the injection task will instead attempt to use internal data products attached to the visit-level dataset.
    This may result in degraded performance, particularly for the WCS solution, or may prevent data reduction from progressing at all for injection into some early dataset types.

To resolve the CRITICAL error, check your input collections to ensure that they contain a ``finalVisitSummary`` dataset from a prior data-reduction run.
You can query the butler registry to find all available ``finalVisitSummary`` datasets in your data repository.
For example, to query for all ``finalVisitSummary`` tables for HSC visit 1228:

.. code-block:: shell

    butler query-datasets $REPO finalVisitSummary \
    --where "instrument='HSC' AND visit=1228"

*where*

    `$REPO`
        The path to the butler repository.

.. _lsst.source.injection-faqs-injection-task-only:

How can I run an injection task without making a larger pipeline YAML file?
===========================================================================

Although it is strongly recommended to make use of the instructions in :ref:`lsst.source.injection-ref-make` to add your injection task into an existing pipeline, it is of course possible to run an injection task on its own.
The benefit of using the :doc:`make_injection_pipeline <scripts/make_injection_pipeline>` command line script or the associated :py:func:`~lsst.source.injection.make_injection_pipeline` Python function are that a number of config overrides and task exclusions for downstream pipeline tasks are automatically set for you.

To just run a source injection task on its own, you may either directly call the source injection YAML file in the ``$SOURCE_INJECTION_DIR`` directory when calling ``pipetask run`` on the command line, or call the injection tasks directly from within Python (see :ref:`lsst.source.injection-ref-inject-python` for more information on this latter approach).

For example, to run the ``inject_exposure`` task on its own from the command line, replace the ``-p`` argument in :ref:`lsst.source.injection-ref-inject-cli` with the path to the ``inject_exposure.yaml`` file:

.. code-block:: shell

    ...
    -p $SOURCE_INJECTION_DIR/pipelines/inject_exposure.yaml
    ...

*where*

    `$SOURCE_INJECTION_DIR`
        The path to the source injection package directory.

.. tip::

    Even if you only want to run a source injection task on its own now, there's still nothing preventing you from first constructing a unified pipeline YAML.
    Afterwards, the source injection task can be selected from this YAML in isolation, allowing you to just run that task alone.

    For example, supposing you've constructed the pipeline YAML '``DRP-RC2+injection.yaml``' as in the :ref:`lsst.source.injection-ref-make` instructions, the ``inject_exposure`` task can be isolated by using:

    .. code-block:: shell

        ...
        -p DRP-RC2+injection.yaml#inject_exposure
        ...

.. _lsst.source.injection-faqs-catalog-ingest-requirement:

Do I need to ingest my input injection catalog into the data butler before running source injection?
====================================================================================================

No, you do not need to ingest injection catalogs into the data butler prior to source injection.
It is possible to feed a list of `astropy.table.Table` objects directly into your injection task from within a Python environment.
For more instructions, see the notes in :ref:`lsst.source.injection-ref-inject-python` (skipping the sections which construct an astropy table).

Note that it is not usually desirable to do operate in this manner however, as the output injected catalog may not necessarily contain a full copy of all potential injection sources.
For example, if you construct a source injection catalog that spans an entire tract, but then only perform source injection for a single patch, only injected sources spatially close to that patch will have their input catalog data persisted.
If you then lose your original astropy table, that data will be lost forever.
For this reason, it is usually recommended that your input catalog is ingested first, to facilitate future data analysis and provenance checking.

Please also note that it is not currently possible to perform source injection on the command line without first ingesting an input catalog into the data butler.

.. _lsst.source.injection-faqs-write-into-butler:

Is it possible to write injected exposures back into a data butler?
===================================================================

When working in a command line environment and using ``pipetask run``, injected outputs are automatically written back into the data butler.
However, when working in a Python/Jupyter notebook environment, calling injection tasks directly and running them will not ingest any outputs back into the data butler.

Whilst it is *possible* to write bespoke injected outputs into the data butler, e.g., to facilitate onward processing downstream, this is not normally recommended.
Rather, injecting sources into data within a notebook envionment is designed for quick-look analyses alone.
Datasets ingested into the butler have strict format requirements, without which it would be impossible to guarantee the integrity of a given pipeline.

If you do wish to perform source injection inside a Python/Jupyter notebook environment and store the outputs in the data butler for subsequent processing, it is strongly recommended to make use of the `SimplePipelineExecutor`.
The `SimplePipelineExecutor` is a lightweight high-level executor for running pipelines in Python.
The example below demonstrates how to use the `SimplePipelineExecutor` to run a source injection task in a pipeline already constructed (see :ref:`lsst.source.injection-ref-make`) and write the outputs back into the data butler:

    .. code-block:: python

        import os
        from lsst.ctrl.mpexec import SimplePipelineExecutor
        from lsst.pipe.base import Pipeline

        user = os.getenv("USER")

        # Create a Butler instance with collections appropriate for processing.
        butler = SimplePipelineExecutor.prep_butler(
            REPO,
            inputs=["HSC/runs/RC2/w_2024_10/DM-43178", "injection/defaults"],
            output=f"u/{user}/my_injected_outputs",
        )

        # Load the pipeline from a YAML file.
        pipeline = Pipeline.fromFile(os.path.join("PATH", "TO", "DRP-RC2+injection.yaml#inject_coadd"))

        # Optionally configure a task.
        pipeline.addConfigOverride(
            label="inject_coadd",
            key="process_all_data_ids",
            value=True,
        )

        # Check that the config has been applied.
        inject_coadd_task = pipeline.to_graph().tasks.get("inject_coadd")
        print(inject_coadd_task.config.process_all_data_ids)  # returns: True

        # Create an executor; build a QuantumGraph from an in-memory pipeline.
        executor = SimplePipelineExecutor.from_pipeline(
            pipeline=pipeline,
            where="instrument='HSC' AND skymap='hsc_rings_v1' AND tract=9813 AND patch=42 AND band='i'",
            butler=butler,
        )

        # Run all quanta in the QuantumGraph.
        quanta = executor.run(register_dataset_types=True)
        print(f"number of quanta executed: {len(quanta)}")

        # Get the run outputs.
        dataset_refs = {dtype.name:ref[0] for dtype, ref in quanta[0].outputs.items()}
        print(f"available dataset types: {list(dataset_refs.keys())}")
        injected_deepCoadd = butler.get(dataset_refs["injected_deepCoadd"])

.. _lsst.source.injection-faqs-restrictions-guidelines-conventions:

Are there restrictions or guidelines for ingesting source injection images and/or catalogs into the butler?
===========================================================================================================

Any imaging data product produced by the LSST Science Pipelines can, in theory, have synthetic sources injected into it.
There are currently no restrictions on this, but this may change as the user base performing source injection grows.

With regards guidelines, we recommend that users maintain the prefix naming conventions of '``injection_``' for input source injection data, and '``injected_``' for output source injection data.
This level of consistency should help to avoid confusion when working with multiple data collections and with other collaborators.
However, these naming conventions are not strictly enforced, and you are free to use whatever naming convention you prefer.

.. _lsst.source.injection-faqs-find-more-information:

Where can I find more information on source injection?
======================================================

For more information, consult a :ref:`quick reference guide <lsst.source.injection-ref>` or head back to the :ref:`main page <lsst.source.injection>`.
