.. py:currentmodule:: lsst.source.injection

.. _lsst.source.injection:

=======================
 lsst.source.injection
=======================

----------------------------
 Synthetic Source Injection
----------------------------

The ``lsst.source.injection`` package contains tools designed to assist in the
injection of synthetic sources into scientific imaging. Source generation and
subsequent source injection is powered by the
`GalSim`_ software package.

.. _GalSim: https://galsim-developers.github.io/GalSim/

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

Documentation for ``lsst.pipe.base`` pipeline tasks within this package.

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
