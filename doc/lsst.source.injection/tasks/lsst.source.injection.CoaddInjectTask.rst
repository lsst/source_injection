.. lsst-task-topic:: lsst.source.injection.CoaddInjectTask

=================
 CoaddInjectTask
=================

-----------------------------------------------
 Inject Synthetic Sources Into Coadded Imaging
-----------------------------------------------

``CoaddInjectTask`` injects synthetic sources into coadd-level science imaging.
A user-supplied input catalog defines the positions and characteristics for
synthetic sources to be generated.
If the user opts instead to inject pre-built FITS postage stamp images, the
path to these must be supplied.
The `GalSim`_ software package is used to generate synthetic postage stamps
(if required) before performing source injection into an image.

.. _GalSim: https://galsim-developers.github.io/GalSim/

As an output, an `injected_` variant of the input dataset type will be created.
For coadd-level imaging, the default input is a `deepCoadd` dataset type,
producing an `injected_deepCoadd` output.

.. _lsst.source.injection.CoaddInjectTask-summary:

Processing Summary
==================

``CoaddInjectTask`` runs this sequence of operations:

1. the input injection catalog is loaded and column names standardized;
2. injection sources are cleaned according to certain criteria, including:

  * missing or erroneous magnitude information;
  * centroids that do not fall close to the injection area;
  * sources not selected by virtue of their visit or selection flag;
  * a source type which is not a supported `GSObject`_ class;
  * extreme Sérsic indices (i.e., outside the range :math:`0.3 \le n \le 6.2`);

3. the input source is converted into a `GalSim`_ object:

  * no attempt at source injection will be made if the fully realized bounding
    box does not overlap with the injection image;

4. sources are injected into the coadd injection image using `GalSim`_;
5. provenance metadata is added to the resultant source injected dataset;
6. the injected image and associated injected catalog are output.

.. _GSObject: https://galsim-developers.github.io/GalSim/_build/html/sb.html

.. _lsst.source.injection.CoaddInjectTask-api:

Python API Summary
==================

.. lsst-task-api-summary:: lsst.source.injection.CoaddInjectTask

.. _lsst.source.injection.CoaddInjectTask-configs:

Configuration Fields
====================

.. lsst-task-config-fields:: lsst.source.injection.CoaddInjectTask

.. _lsst.source.injection.CoaddInjectTask-examples:


Examples
========

An example calling this task from the command line using `pipetask run`:

.. code-block:: shell

    SKYINST="skymap='hsc_rings_v1' AND instrument='HSC'"
    DATAQUERY="tract=9813 AND patch=42 AND band='i'"

    pipetask --long-log --log-file $LOGFILE \
    run --register-dataset-types \
    -b $REPO \
    -i $INPUT_DATA_COLL,$INJECTION_CATALOG_COLL \
    -o $OUTPUT_COLL \
    -p $SOURCE_INJECTION_DIR/pipelines/inject_coadd.yaml \
    -d "$SKYINST AND $DATAQUERY"

*where*

    `$LOGFILE`
        The full path to a user-defined output log file.

    `$REPO`
        The path to the butler repository.

    `$INPUT_DATA_COLL`
        The name of the input data collection.

    `$INJECTION_CATALOG_COLL`
        The name of the input injection catalog collection.

    `$OUTPUT_COLL`
        The name of the injected output collection.

    `$SOURCE_INJECTION_DIR`
        The path to the source injection module directory.

.. _lsst.source.injection.CoaddInjectTask-debug:

Debugging
=========

Additional debug-level log information may be shown by setting the `log-level` option to `DEBUG` at runtime:

.. code-block:: shell

    pipetask --log-level DEBUG ...
