.. py:currentmodule:: lsst.source.injection

.. _lsst.source.injection:

=======================
 lsst.source.injection
=======================

----------------------------
 Synthetic Source Injection
----------------------------

The ``lsst.source.injection`` package contains tools designed to assist in the injection of synthetic sources into scientific imaging.
Source generation and subsequent source injection is powered by the `GalSim`_ software package.

.. _GalSim: https://galsim-developers.github.io/GalSim/

.. figure:: _assets/t9813p42i_zoom_sersic_prepost_injection.gif
    :name: t9813p42i_zoom_sersic_prepost_injection
    :alt: An HSC i-band cutout from tract 9813, patch 42, showing the injection of a series of synthetic Sérsic sources.
    :align: center
    :width: 100%

    ..

    An HSC i-band cutout from tract 9813, patch 42, showcasing the injection of a series of synthetic Sérsic sources.
    Images are ~100 arcseconds on the short axis, log scaled across the central 99.5% flux range, and smoothed with a Gaussian kernel of FWHM 3 pixels.

    .. list-table::
        :widths: 1 1

        * - .. figure:: _assets/t9813p42i_zoom_sersic_pre_injection.png
                :name: t9813p42i_zoom_sersic_pre_injection
                :alt: Tract 9813, patch 42, HSC i-band cutout, before source injection.
                :align: center
                :width: 100%

                ..

                Before injection.
          - .. figure:: _assets/t9813p42i_zoom_sersic_post_injection.png
                :name: t9813p42i_zoom_sersic_post_injection
                :alt: Tract 9813, patch 42, HSC i-band cutout, after source injection.
                :align: center
                :width: 100%

                ..

                After injection.

.. _lsst.source.injection-ref:

Quick Reference Guide
=====================

References for each aspect of the synthetic source injection process.

.. toctree::
    :maxdepth: 1
    :glob:

    reference/*

.. _lsst.source.injection-questions:

Frequently Asked Questions
==========================

A collection of :ref:`lsst.source.injection-faqs`, and answers.

.. _lsst.source.injection-scripts:

Command Line Scripts
====================

Utility scripts which may be called from the command line.

.. toctree::
    :maxdepth: 1
    :glob:

    scripts/*

.. _lsst.source.injection-tasks:

Pipeline Tasks
==============

Documentation for ``lsst.source.injection`` pipeline tasks within this package.

.. lsst-pipelinetasks::
    :root: lsst.source.injection
    :toctree: tasks

.. _lsst.source.injection-contributing:

Contributing
============

The ``lsst.source.injection`` package is developed at
`github.com/lsst/source_injection <https://github.com/lsst/source_injection>`_.

Jira issues relating to this package can be found using the
`source_injection <https://ls.st/sourceinjectionjira>`_ component.

.. _lsst.source.injection-pyapi:

Python API Reference
====================

.. automodapi:: lsst.source.injection
