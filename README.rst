################
source_injection
################

The ``source_injection`` package contains tools designed to assist in the injection of synthetic sources into scientific imaging.
Source generation and subsequent source injection is powered by the `GalSim`_ software package.

.. _GalSim: https://galsim-developers.github.io/GalSim/

For more information, see the ``source_injection`` documentation on `pipelines.lsst.io`_.

.. _pipelines.lsst.io: https://pipelines.lsst.io/v/daily/modules/lsst.source.injection

.. figure:: doc/lsst.source.injection/_assets/t9813p42i_zoom_sersic_prepost_injection.gif
    :name: t9813p42i_zoom_sersic_prepost_injection
    :alt: An HSC i-band cutout from tract 9813, patch 42, showing the injection of a series of synthetic Sérsic sources.
    :align: center
    :width: 100%

    ..

    An HSC i-band cutout from tract 9813, patch 42, showcasing the injection of a series of synthetic Sérsic sources.
    Images are ~100 arcseconds on the short axis, log scaled, and smoothed with a Gaussian kernel.

    .. list-table::
        :widths: 1 1

        * - .. figure:: doc/lsst.source.injection/_assets/t9813p42i_zoom_sersic_pre_injection.png
                :name: t9813p42i_zoom_sersic_pre_injection
                :alt: Tract 9813, patch 42, HSC i-band cutout, before source injection.
                :align: center
                :width: 100%

                ..

                Before injection.
          - .. figure:: doc/lsst.source.injection/_assets/t9813p42i_zoom_sersic_post_injection.png
                :name: t9813p42i_zoom_sersic_post_injection
                :alt: Tract 9813, patch 42, HSC i-band cutout, after source injection.
                :align: center
                :width: 100%

                ..

                After injection.
