.. _lsst.source.injection-ref-make:

============================
 Make an Injection Pipeline
============================

----------------------------------------------------
 Making a Fully Qualified Source Injection Pipeline
----------------------------------------------------

Fully qualified source injection pipeline definition YAML files are dynamically generated prior to use by combining a source injection pipeline definition stub with a reference pipeline file.
This allows the user to specify source injection configuration parameters in a simple, human-readable format and provides flexibility in choosing which dataset type synthetic sources are injected into.
The reference pipeline file must be a complete pipeline definition YAML file, typically one used to reduce data as part of a data reduction campaign.

Either the :doc:`make_injection_pipeline <../scripts/make_injection_pipeline>` command line script or the associated :py:func:`~lsst.source.injection.make_injection_pipeline` Python function may be used to generate a fully qualified injection pipeline.
Examples on this page illustrate the use of both methods.

.. note::

    Two legacy dynamic source injection pipelines are automatically generated inside the ``drp_pipe`` repository.
    These pipelines are located in the ``$DRP_PIPE_DIR/pipelines/HSC`` directory, facilitating source injection data reductions for the Hyper Suprime-Cam RC2 and RC2 subset datasets: ``DRP-RC2+injected_deepCoadd.yaml`` and ``DRP-RC2_subset+injected_deepCoadd.yaml``, respectively.
    As indicated by the appended name, synthetic sources are injected into the ``deepCoadd`` dataset type.

.. _lsst.source.injection-ref-make-stubs:

Injection Pipeline Stubs
=========================

A number of different source injection pipeline stubs have been constructed in the ``$SOURCE_INJECTION_DIR/pipelines`` directory.
Each of these pipeline stubs contain a single task that is used to inject sources into a particular dataset type.

Although these injection pipeline YAML stubs can be used directly, it is recommended that the :doc:`make_injection_pipeline <../scripts/make_injection_pipeline>` command line script or the associated :py:func:`~lsst.source.injection.make_injection_pipeline` Python function be used to generate a complete source injection pipeline definition YAML file for subsequent use.
A complete injection pipeline definition file will contain the pipeline stub as a subtask alongside any additional tasks required to complete the source injection process.
Tasks from the reference pipeline may either be removed or have specific configuration overrides applied as necessary to support subsequent injected source image data reduction.

.. note::

    When using the above utilities to construct a fully qualified injection pipeline, any existing subsets will also be updated to include the injection task where appropriate.
    Furthermore, a series of ``injected_*`` subsets will be constructed.
    These ``injected_*`` subsets are copies of existent subsets, but with any tasks not directly impacted by source injection removed.

    For example, if the ``inject_exposure.yaml`` pipeline stub is used to inject sources into a ``post_isr_image`` dataset type, the subset of the reference pipeline containing the ``isr`` task will be updated to also include the ``injectExposure`` task.

    This behavior can be disabled by passing the ``-e`` argument on the command line, or setting ``exclude_subsets`` to ``True`` in Python.
    Additionally, a new subset, ``injected_[MY_SUBSET]``, will also be created containing all tasks from the ``[MY_SUBSET]`` subset but with the ``isr`` task removed (as sources will be injected after this task has run).

.. note::

    After a fully qualified injection pipeline has been generated, a check is performed to ensure that all reference :ref:`pipeline contracts <pipeline_creating_contracts>` (if any) are satisfied.
    Pipeline contracts are a means by which to ensure that certain configuration values are set in a predictable manner.
    When generating an injection pipeline, it's possible that some of these contracts will become invalid.
    For example, if a contract specifies that the dataset type produced by a task prior to source injection matches the dataset type consumed by a task after source injection, this contract may become invalid if the tasks downstream of source injection have been modified to instead consume the new source injected input.
    The :doc:`make_injection_pipeline <../scripts/make_injection_pipeline>` command line script and the :py:func:`~lsst.source.injection.make_injection_pipeline` Python function will check for this and warn if any contracts are invalid.
    Invalid contracts will be removed from the final output pipeline YAML.

The table below lists the available pipeline YAML stubs inside the ``$SOURCE_INJECTION_DIR/pipelines`` directory and the dataset types they are designed to inject sources into:

.. list-table::
    :widths: 1 1 1 1
    :stub-columns: 1

    * - Injection Dataset Type
      - Exposure-level image (e.g. ``post_isr_image``)
      - Visit-level image (e.g. ``visit_image``)
      - Coadd-level image (e.g. ``deep_coadd_predetection``)
    * - Injection Pipeline Stub
      - inject_exposure.yaml_
      - inject_visit.yaml_
      - inject_coadd.yaml_
    * - Injection Task
      - :lsst-task:`~lsst.source.injection.ExposureInjectTask`
      - :lsst-task:`~lsst.source.injection.VisitInjectTask`
      - :lsst-task:`~lsst.source.injection.CoaddInjectTask`
    * - Injection Task Graph
      - .. image:: ../_assets/inject_exposure.png
            :width: 100%
      - .. image:: ../_assets/inject_visit.png
            :width: 100%
      - .. image:: ../_assets/inject_coadd.png
            :width: 100%
    * -
      - :download:`PDF <../_assets/inject_exposure.pdf>`
      - :download:`PDF <../_assets/inject_visit.pdf>`
      - :download:`PDF <../_assets/inject_coadd.pdf>`

.. _inject_exposure.yaml: https://github.com/lsst/source_injection/blob/main/pipelines/inject_exposure.yaml
.. _inject_visit.yaml: https://github.com/lsst/source_injection/blob/main/pipelines/inject_visit.yaml
.. _inject_coadd.yaml: https://github.com/lsst/source_injection/blob/main/pipelines/inject_coadd.yaml

A source injection pipeline stub may always be specified directly, however, both the :doc:`make_injection_pipeline <../scripts/make_injection_pipeline>` command line script and the :py:func:`~lsst.source.injection.make_injection_pipeline` Python function will attempt to infer the correct pipeline stub to use based on the injected dataset type specified.
This inference is based on a match of the injected dataset type to a predefined list of common types and their associated pipeline stubs.

.. _lsst.source.injection-ref-make-cli:

Make an Injection Pipeline on the Command Line
==============================================

The :doc:`make_injection_pipeline <../scripts/make_injection_pipeline>` command line script is used to generate a complete source injection pipeline definition YAML file.
More information on the operation of this script may be obtained by running ``make_injection_pipeline --help``.

As an example on the command line, to create a pipeline YAML which will inject a synthetic source into a `post_isr_image` exposure-type dataset type using the LSSTCam DRP pipeline as a reference:

.. code-block:: shell

    make_injection_pipeline \
    -t post_isr_image \
    -r $DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml \
    -f DRP-injection.yaml

*where*

    `$DRP_PIPE_DIR`
        The path to the `drp_pipe` package directory.

The above command will save a complete and fully expanded pipeline definition file into the file ``DRP-injection.yaml``.
In this example, synthetic sources are to be injected into the ``post_isr_image`` dataset type, using the ``LSSTCam/DRP.yaml`` pipeline definition file as a reference.
As the ``post_isr_image`` dataset type has dimensions of ``exposure``, the ``inject_exposure.yaml`` source injection pipeline definition file stub has been automatically inferred.
That particular injection pipeline YAML stub contains the :lsst-task:`~lsst.source.injection.ExposureInjectTask` task.

.. tip::

    To print the fully qualified output pipeline to the terminal window instead of saving it to a file, omit the ``-f`` option in the above example.

To specify an injection pipeline definition file stub explicitly rather than allowing the function to attempt to infer it from the injected dataset type, the ``-i`` option may be appended to the above command:

.. code-block:: shell

    ...
    -i $SOURCE_INJECTION_DIR/pipelines/inject_exposure.yaml

*where*

    `$SOURCE_INJECTION_DIR`
        The path to the source injection package directory.

.. _lsst.source.injection-ref-make-python:

Make an Injection Pipeline in Python
====================================

The :py:func:`~lsst.source.injection.make_injection_pipeline` Python function is used to generate a complete source injection pipeline definition YAML file in Python:

.. code-block:: python

    from lsst.source.injection import make_injection_pipeline

More information on the operation of this function may be obtained by calling ``make_injection_pipeline?`` in a Python interpreter.

As an example in Python, to create a pipeline which will inject a synthetic source into a `post_isr_image` exposure-type dataset type using the LSSTCam DRP pipeline as a reference:

.. code-block:: python

    # Construct the Pipeline object.
    pipeline = make_injection_pipeline(
        dataset_type_name="post_isr_image",
        reference_pipeline="$DRP_PIPE_DIR/pipelines/LSSTCam/DRP.yaml",
    )

    # Print the pipeline.
    print(pipeline)

To specify an injection pipeline definition file stub explicitly rather than attempting to infer it from the injected dataset type, the ``injection_pipeline`` argument may also be used, e.g.:

.. code-block:: python

    pipeline = make_injection_pipeline(
        ...
        injection_pipeline="$SOURCE_INJECTION_DIR/pipelines/inject_exposure.yaml",
    )

Once a pipeline object has been constructed, it may be written to disk using the ``write_to_uri`` method:

.. code-block:: python

    pipeline.write_to_uri("DRP-injection.yaml")

.. _lsst.source.injection-ref-make-visualize:

Visualize an Injection Pipeline
===============================

Any pipeline YAML, including an injection pipeline, can be visualized to clarify exactly what the pipeline does.
In this section we provide instructions for visualizing the ``DRP-injection.yaml`` pipeline generated in the above examples.
Options for text-based outputs on the command line and rich rendered outputs are presented.
The tasks and dataset types printed below are accurate as of ``w_2025_37`` of the LSST Science Pipelines.

.. tip::

    Only the ``isr``, ``injectExposure`` and ``calibrateImage`` tasks of the
    fully qualified injection pipeline are selected in the snippets below by
    appending the ``#`` symbol followed by a comma-separated list of the task
    label names to the YAML pipeline filename.
    Any subset or task within a pipeline YAML can be selected in this way.

.. _lsst.source.injection-ref-make-visualize-tasks:

Visualize pipeline tasks
------------------------

The snippet below will generate a text-based representation of only the tasks in the ``isr``, ``injectExposure`` and ``calibrateImage`` tasks from the pipeline.

.. code-block:: shell

    pipetask build \
    -p DRP-injection.yaml#isr,injectExposure,calibrateImage \
    --show task-graph

returning:

.. code-block:: shell

    ■  isr: {detector, exposure}
    │
    ■  injectExposure: {detector, exposure}
    │
    ■  calibrateImage: {detector, visit}

.. _lsst.source.injection-ref-make-visualize-pipeline:

Visualize pipeline tasks and datasets
-------------------------------------

The snippet below will generate a text-based representation of both the tasks and the input/output dataset types associated with the ``isr``, ``injectExposure`` and ``calibrateImage`` tasks from the pipeline.

.. code-block:: shell

    pipetask build \
    -p DRP-injection.yaml#isr,injectExposure,calibrateImage \
    --show pipeline-graph

returning:

.. code-block:: shell

                      ○  flat: {detector, physical_filter} ExposureF
                      │
                    ○ │  bfk: {detector} BrighterFatterKernel
                    │ │
                  ○ │ │  camera: {instrument} Camera
                  │ │ │
                ○ │ │ │  crosstalk: {detector} CrosstalkCalib
                │ │ │ │
              ○ │ │ │ │  cti: {detector} IsrCalib
              │ │ │ │ │
            ◍ │ │ │ │ │  dark, bias: {detector} ExposureF
            │ │ │ │ │ │
          ○ │ │ │ │ │ │  defects: {detector} Defects
          │ │ │ │ │ │ │
        ○ │ │ │ │ │ │ │  linearizer: {detector} Linearizer
        │ │ │ │ │ │ │ │
      ○ │ │ │ │ │ │ │ │  ptc: {detector} PhotonTransferCurveDataset
      │ │ │ │ │ │ │ │ │
    ○ │ │ │ │ │ │ │ │ │  raw: {detector, exposure} Exposure
    ╰─┴─┴─┴─┴─┴─┴─┴─┴─┤
                      ■  isr: {detector, exposure}
                    ╭─┤
                    ○ │  isrStatistics: {detector, exposure} StructuredDataDict
                      │
                      ○  post_isr_image: {detector, exposure} Exposure
                      │
                    ○ │  injection_catalog: {band, htm7} ArrowAstropy
                    │ │
                  ○ │ │  visit_summary: {visit} ExposureCatalog
                  ╰─┴─┤
                      ■  injectExposure: {detector, exposure}
                    ╭─┤
                    ○ │  injected_post_isr_image_catalog: {detector, exposure}...[1]
                      │
                      ○  injected_post_isr_image: {detector, exposure} Exposure
                      │
                    ○ │  the_monster_20250219: {htm7} SimpleCatalog
                    ╰─┤
                      ■  calibrateImage: {detector, visit}
                    ╭─┤
                    ○ │  injected_preliminary_visit_image: {detector, visit} E...[2]
                    ╭─┤
                    ○ │  injected_preliminary_visit_image_background: {detecto...[3]
                    ╭─┤
                    ◍ │  injected_single_visit_star_footprints, injected_singl...[4]
                    ╭─┤
                    ◍ │  injected_single_visit_star_unstandardized, injected_s...[5]
                      │
                      ◍  injected_initial_photometry_match_detector, injected_...[6]
    [1]
      injected_post_isr_image_catalog: {detector, exposure} ArrowAstropy
    [2]
      injected_preliminary_visit_image: {detector, visit} ExposureF
    [3]
      injected_preliminary_visit_image_background: {detector, visit} Background
    [4]
      injected_single_visit_star_footprints,
      injected_single_visit_psf_star_footprints: {detector, visit} SourceCatalog
    [5]
      injected_single_visit_star_unstandardized, injected_single_visit_psf_star:
      {detector, visit} ArrowAstropy
    [6]
      injected_initial_photometry_match_detector,
      injected_initial_astrometry_match_detector: {detector, visit} Catalog

.. _lsst.source.injection-ref-make-visualize-render:

Render a pipeline in graphical format
-------------------------------------

The ``pipetask build`` command can also output a pipeline in the GraphViz DOT graph description language format.
This format can be rendered into multiple visual formats such as PDF or PNG types using the ``dot`` command line tool.

The snippet below converts the ``isr``, ``injectExposure`` and ``calibrateImage`` tasks from the pipeline produced in the above example into a PNG file.
To help improve the layout of the graph, the ``unflatten`` preprocessing filter is also used.

.. code-block:: shell

    INPUT_PIPELINE=DRP-injection.yaml#isr,injectExposure,calibrateImage
    OUTPUT_FILE=DRP_with_injected_exposure.png
    OUTPUT_EXT=${OUTPUT_FILE##*.}  # Resolves to: pdf/svg/png/jpg/...

    # Create the directed graph from an input pipeline.
    pipetask build -p $INPUT_PIPELINE --pipeline-dot graph_pre.dot

    # Post-process the directed graph to improve layout.
    unflatten -l 3 -f -o graph_post.dot graph_pre.dot

    # Draw the directed graph.
    dot graph_post.dot -T$OUTPUT_EXT > $OUTPUT_FILE
    # NB: Also add an optional -Gdpi=[DPI] argument to change the resolution

The output PNG from the above example injection into a ``post_isr_image`` type is shown below (left panel).
Equivalent graphs for injections into ``preliminary_visit_image`` (central panel) and ``deep_coadd_predetection`` (right panel) types are also shown, for reference.

.. list-table::
    :widths: 1 1 1

    * - .. image:: ../_assets/DRP_with_injected_exposure.png
            :width: 100%
      - .. image:: ../_assets/DRP_with_injected_visit.png
            :width: 100%
      - .. image:: ../_assets/DRP_with_injected_coadd.png
            :width: 100%
    * - :download:`PDF <../_assets/DRP_with_injected_exposure.pdf>`
      - :download:`PDF <../_assets/DRP_with_injected_visit.pdf>`
      - :download:`PDF <../_assets/DRP_with_injected_coadd.pdf>`
    * - The ``injectExposure`` task merged into the LSSTCam DRP pipeline.
      - The ``injectVisit`` task merged into the LSSTCam DRP pipeline.
      - The ``injectCoadd`` task merged into the LSSTCam DRP pipeline.

.. _lsst.source.injection-ref-make-wrap:

Wrap Up
=======

This reference page has described how to make a fully qualified source injection pipeline definition YAML file, either on the command line or in Python.
Options for visualizing the resultant pipeline have also been presented.

Move on to :ref:`another quick reference guide <lsst.source.injection-ref>`, consult the :ref:`FAQs <lsst.source.injection-faqs>`, or head back to the `main page <..>`_.
