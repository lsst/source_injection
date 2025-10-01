.. _lsst.source.injection-ref-generate:

===============================
 Generate an Injection Catalog
===============================

----------------------------------------
 Generating Synthetic Source Parameters
----------------------------------------

The LSST Science Pipelines provides a set of tools in the `source_injection` package to assist in generating synthetic source injection catalogs.
Such synthetic catalogs may be used as part of quality assurance or algorithmic performance testing processes.

To generate a synthetic source injection catalog, we first need to define the parameters of the sources we want to inject.
These parameters include the position, flux, size, and shape of the sources.
Both the :doc:`generate_injection_catalog <../scripts/generate_injection_catalog>` command line script or the associated :py:func:`~lsst.source.injection.generate_injection_catalog` Python function may be used to generate such a catalog.
Examples on this page illustrate the use of both methods.

.. note::

    These source generation scripts are meant as a convenience, rather than a requirement.
    If you prefer, you may generate your own catalog using any method you choose.

.. _lsst.source.injection-ref-generate-catalog:

Catalog Parameters
==================

The tools on this page assist in generating a synthetic source injection catalog with sources quasi-randomly distributed across the sky.
Sources will be scattered between the right ascension and declination limits specified by the user (using, for example, the `-a` and `-d` arguments respectively on the command line).
Optional magnitudes may also be generated using the same sequence (using the `-m` command line argument).
The number of repeats for each unique combination of profile parameters is also specified by the user (using, for example, the `-n` argument on the command line).

.. tip::

    See :ref:`lsst.source.injection-faqs-spatial-overlap-1` and :ref:`lsst.source.injection-faqs-spatial-overlap-2` for instructions on calculating spatial overlap between your input injection catalog and available on-disk data.

Other than this, all additional parameters which describe the surface brightness profile of sources to be injected are specified using the `-p` argument on the command line, or as additional keyword arguments in Python.
For example, to generate a catalog of Sérsic sources that have three different magnitudes (15, 17, and 19), three different values of the Sérsic index (1, 2, and 4) and two different half light radii (5 and 10), you would use the following on the command line:

.. code-block:: shell

    ...
    -p source_type Sersic \
    -p mag 15 17 19 \
    -p n 1 2 4 \
    -p half_light_radius 5 10

Magnitudes correspond only to one band.
If source injections are to be performed in multiple bands, multiple catalogs (one per band) must be generated.
If the magnitude of injected sources does not vary as a function of band however, a single catalog can be associated with multiple bands at the time of ingestion into the butler repository (see :ref:`lsst.source.injection-ref-ingest` for more information).

.. _lsst.source.injection-ref-generate-catalog-units:

Units
-----

By default, units are assumed to be those expected by `GalSim <https://galsim-developers.github.io/GalSim/_build/html/sb.html>`_.
Typically, this will be units of pixels for simple profile length measures, arcseconds for astrophysical length measures (such as the Sérsic profile) and degrees for angle measures.
Exceptions to this are special parameters reserved for use by this package (see :ref:`lsst.source.injection-ref-generate-catalog-parameters` below for more details) such as ``mag``, which is in units of AB magnitude.

To specify a parameter in units other than those expected by GalSim, values may be multiplied by an astropy unit object when working in Python (this cannot be achieved from the command line).
For example, to specify a half light radius of 5 arcseconds, you may use something similar to ``half_light_radius=5*u.arcsec``.

.. _lsst.source.injection-ref-generate-catalog-parameters:

Source Injection Package Parameters
-----------------------------------

There are two types of input catalog parameters: special column names reserved for use by this package, and those defined by GalSim.
The table below lists the parameter names reserved for use by this package.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Parameter Name
      - Description
    * - ``ra``
      - The right ascension of the source in degrees.
    * - ``dec``
      - The declination of the source in degrees.
    * - ``source_type``
      - The `GalSim Class of surface brightness profile type <https://galsim-developers.github.io/GalSim/_build/html/sb.html>`_ to use (see **note** below).
    * - ``mag``
      - The magnitude of the source in a single band (converted to `flux` for use by GalSim`).
    * - ``stamp``
      - The filename for a local FITS postage stamp image to inject instead of using GalSim to draw a source.
    * - ``draw_size``
      - The size of the bounding box within which a synthetic source is generated. If not provided, GalSim will determine an optimal draw size based on the profile parameters.
    * - ``trail_length``
      - The length of the a satellite trail in pixels.

.. note::

    Injection of FITS-file postage stamps only requires the ``ra``, ``dec``, ``source_type``, ``stamp`` and ``mag`` columns to be specified in the injection catalog.
    The ``source_type`` values should all be set to ``Stamp`` (see below).
    For more information on injection of postage stamps, see :ref:`lsst.source.injection-ref-inject-stamps`.

.. _lsst.source.injection-ref-generate-catalog-types:

Source Types
------------

With regards the ``source_type`` parameter, most of the `surface brightness profile parameters defined by GalSim <https://galsim-developers.github.io/GalSim/_build/html/sb.html>`_ are natively supported.
This includes standard classes such as `Sersic`, `Exponential`, `DeVaucouleurs`, `DeltaFunction` and `Gaussian`.
The `source_injection` package additionally provides these extra custom classes:

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Custom Source Type
      - Description
    * - ``Stamp``
      - A FITS postage stamp image of a source.
    * - ``Trail``
      - A satellite trail.
    * - ``Star``
      - A convenient alias for the GalSim `DeltaFunction <https://galsim-developers.github.io/GalSim/_build/html/simple.html#delta-function>`_ class.

To see a full list of supported GalSim surface brightness profile types and their allowed parameters, call the :doc:`show_source_types <../scripts/show_source_types>` command line script.

.. _lsst.source.injection-ref-generate-catalog-galsim:

Common GalSim Parameters
------------------------

The table below is a *non-exhaustive* list of some of the most commonly used GalSim surface brightness profile parameters.
See the `full GalSim surface brightness profile documentation <https://galsim-developers.github.io/GalSim/_build/html/sb.html>`_ for more details.
Note that these parameter names must be written exactly as defined in the GalSim documentation in order for them to be correctly identified.

.. list-table::
    :widths: auto
    :header-rows: 1

    * - Parameter Name
      - Description
    * - ``n``
      - Sérsic index
    * - ``half_light_radius``
      - The half-light radius of the source.
    * - ``q``
      - The minor-to-major axis ratio.
    * - ``beta``
      - The position angle of the object (in degrees).

.. _lsst.source.injection-ref-generate-catalog-ids:

Injection IDs
-------------

When using the tools below to construct a synthetic source injection catalog, a unique ID is assigned to each source under the ``injection_id`` column.
This ID may be used to uniquely identify an injected source and groups of associated injected sources in the output catalog.
*Injection IDs are for reference only* and are not used by the source injection process itself.
If so desired, this column may be replaced with user-defined values instead.

If the number of repeats for each unique combination of injection parameters is 1 (default), then the automatically generated ``injection_id`` values will start at zero and increase by one for each source (i.e., ``0, 1, 2, 3, ...``).

If the number of repeats for each unique combination of injection parameters is greater than 1 however, then ``injection_id`` values are given by base (N) + i, where i is the index of the repeat for a given unique combination of injection parameters and N is the base required for injection ID values to not clash.
For example, if 3 repeats for each unique parameter combination are requested, injection IDs will be: ``0, 1, 2, 10, 11, 12, 20, 21, 22, ...``.
If however 20 repeats for each unique parameter combination are requested, injection IDs will be: ``0, 1, 2, ..., 17, 18, 19, 100, 101, 102, ...``.

.. _lsst.source.injection-ref-generate-cli:

Generate an Injection Catalog on the Command Line
=================================================

The :doc:`generate_injection_catalog <../scripts/generate_injection_catalog>` command line script is used to generate a synthetic source injection catalog.
This script takes a number of arguments, including the right ascension and declination limits of the quasi-randomly generated positions and the number of sources to inject.
Optional magnitudes may also be generated using the same random sequence.
More information on the arguments accepted by this script may be found by running:

.. code-block:: shell

    generate_injection_catalog --help

The example below generates a synthetic source injection catalog with sources randomly scattered in the range 149.7 < RA < 150.1 and 2.0 < Dec < 2.4, with 3 repeats of each unique combination of profile parameters.
Additional parameters describing a series of Sérsic sources are also specified (see above for more details).

.. code-block:: shell

    generate_injection_catalog \
    -a 149.7 150.1 \
    -d 2.0 2.4 \
    -n 3 \
    -p source_type Sersic \
    -p mag 15 17 19 \
    -p n 1 2 4 \
    -p half_light_radius 5 10

.. _skylimits:

.. note::

    The RA and Dec limits above were chosen to fully overlap HSC tract 9813, patch 42; a tract in the COSMOS field.
    These limits were also designed to overlap with HSC i-band visit 1228, detectors 42, 43, 50, 51, 58 and 59.

Running the above will generate a catalog containing 54 sources: 18 combinations repeated 3 times, of which the first several lines will look something like this:

.. _catalogsnippet:

.. code-block:: shell

    injection_id         ra                dec         source_type mag   n  half_light_radius
    ------------ ------------------ ------------------ ----------- ---- --- -----------------
               0  149.8402947814415  2.210198210586508      Sersic 15.0 1.0               5.0
               1  150.0277947814415 2.2250130254013225      Sersic 15.0 1.0               5.0
               2 149.72154478144148 2.3830377167593473      Sersic 15.0 1.0               5.0
              10 149.72779478144147  2.091679692067989      Sersic 15.0 1.0              10.0
              11  149.9965447814415 2.2645191982408286      Sersic 15.0 1.0              10.0
              12 150.09654478144148 2.0571117908334213      Sersic 15.0 1.0              10.0
              20 149.84654478144148 2.1311858649074953      Sersic 15.0 2.0               5.0
              21 149.81529478144148  2.284272284660582      Sersic 15.0 2.0               5.0
              22 149.97779478144147 2.0719266056482364      Sersic 15.0 2.0               5.0
    ...

**To generate source positions using WCS information** (recommended), you may supply the `-b` (butler data repository), `-w` (WCS dataset type), `-c` (collection to query) and, optionally, the `--where` arguments to the script.
With these arguments, a lookup to the butler data repository is made to identify a dataset with WCS appropriate for this catalog.
If these arguments are *not* supplied, source positions will be generated using Cartesian geometry instead.

For example, to use WCS information from a ``deepCoadd_calexp`` dataset from HSC tract 9813, patch 42 in the i-band within the ``HSC/runs/RC2/w_2023_35/DM-40588`` collection, you would add something akin to the following to the above query:

.. code-block:: shell

    ...
    -b $REPO \
    -w deepCoadd_calexp \
    -c HSC/runs/RC2/w_2023_35/DM-40588 \
    --where "instrument='HSC' AND skymap='hsc_rings_v1' AND tract=9813 AND patch=42 AND band='i'"

*where*

    `$REPO`
        The path to the butler repository.

**To save this catalog to a file on disk**, you may add `-f` argument.
Any format type recognized by the astropy Table API may be used.
An attempt to recognize the format will be made based on the file extension (or explicitly specified using the `--format` argument):

.. code-block:: shell

    ...
    -f my_injection_catalog.csv

**To register this catalog directly into the butler**, you may add the `-b` (butler data repository), `-i` (input bands) and `-o` (output collection) arguments when calling the script.
For example, to write the catalog to a butler repository at ``$REPO`` under the collection ``u/$USER/my_injection_inputs`` and to register this catalog to the ``g``, ``r`` and ``i`` bands, you would add these arguments:

.. code-block:: shell

    ...
    -b $REPO \
    -i g r i \
    -o u/$USER/my_injection_inputs

*where*

    `$REPO`
        The path to the butler repository.

    `$USER`
        The users username.

.. _lsst.source.injection-ref-generate-python:

Generate an Injection Catalog in Python
=======================================

The :py:func:`~lsst.source.injection.generate_injection_catalog` Python function is used to generate a synthetic source injection catalog in Python:

.. code-block:: python

    from lsst.source.injection import generate_injection_catalog

More information on the operation of this function may be obtained by calling ``generate_injection_catalog?`` in a Python interpreter.

As an example in Python, the snippet below creates a source injection catalog with synthetic Sérsic sources quasi-randomly scattered in the range 149.7 < RA < 150.1 and 2.0 < Dec < 2.4 (see :ref:`this note for more information on this choice of limits<skylimits>`).
Source combinations consist of three distinct magnitudes, three distinct Sérsic indices and two distinct half light radii.
Three repeats of each unique combination of profile parameters are generated.

.. code-block:: python

    my_injection_catalog = generate_injection_catalog(
        ra_lim=[149.7, 150.1],
        dec_lim=[2.0, 2.4],
        number=3,
        source_type="Sersic",
        mag=[15, 17, 19],
        n=[1, 2, 4],
        half_light_radius=[5, 10],
    )

The resulting catalog is an `astropy.table.Table` object, which may be manipulated as desired (see above for :ref:`an example snippet of the resultant catalog<catalogsnippet>`).
Further arguments may also be supplied to this function.
For example, WCS information may be supplied via the ``wcs`` argument, a source density may be requested via ``density``, and a seed for the random number generator may be supplied via ``seed``.

The resultant catalog may be saved to disk and/or, alternatively, registered into a butler data repository (see :ref:`lsst.source.injection-ref-ingest` for more information).

.. _lsst.source.injection-ref-generate-wrap:

Wrap Up
=======

This page has shown how to generate a synthetic source injection catalog for use with the LSST Science Pipelines.
This catalog may be generated either on the command line or in Python.
Information about the specific format of this catalog are also discussed above.
Further modification of the catalog may also occur prior to its ingestion into the butler repository.

Move on to :ref:`another quick reference guide <lsst.source.injection-ref>`, consult the :ref:`FAQs <lsst.source.injection-faqs>`, or head back to the `main page <..>`_.
