"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
For more information, see:
https://developer.lsst.io/stack/building-single-package-docs.html
"""

from documenteer.conf.pipelinespkg import *  # noqa: F403

project = "source_injection"
html_theme_options["logotext"] = project  # noqa: F405
html_title = project
html_short_title = project

# Update intersphinx configuration to allow finding butler/pipe_base entries.
intersphinx_mapping["lsst"] = ("https://pipelines.lsst.io/v/daily/", None)  # noqa
