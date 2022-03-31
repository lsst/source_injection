# This file is part of source_injection.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ["make_injection_pipeline"]

from lsst.pipe.base import Pipeline


def _get_dataset_type_names(conns, fields):
    """Return the name of a connection's dataset type."""
    dataset_type_names = set()
    for field in fields:
        dataset_type_names.add(getattr(conns, field).name)
    return dataset_type_names


def make_injection_pipeline(
    dataset_type_name: str,
    reference_pipeline: str,
    injection_pipeline: str | None = None,
    exclude_subsets: bool = False,
    prefix: str = "injected_",
    instrument: str | None = None,
) -> Pipeline:
    """Make an expanded source injection pipeline.

    This function takes a reference pipeline definition file in YAML format and
    prefixes all post-injection dataset type names with the injected prefix. If
    an optional injection pipeline definition YAML file is also provided, the
    injection task will be merged into the pipeline.

    Unless explicitly excluded, all subsets from the reference pipeline which
    contain the task which generates the injection dataset type will also be
    updated to include the injection task.

    Parameters
    ----------
    dataset_type_name : `str`
        Name of the dataset type being injected into.
    reference_pipeline : `str`
        Location of a reference pipeline definition YAML file.
    injection_pipeline : `str`, optional
        Location of an injection pipeline definition YAML file.
    exclude_subsets : `bool`, optional
        If True, do not update pipeline subsets to include the injection task.
    prefix : `str`, optional
        Prefix to prepend to each affected post-injection dataset type name.
    instrument : `str`, optional
        Add instrument overrides. Must be a fully qualified class name.

    Returns
    -------
    pipeline : `lsst.pipe.base.Pipeline`
        An expanded source injection pipeline.
    """
    pipeline = Pipeline.fromFile(reference_pipeline)

    # Add an instrument override, if provided.
    if instrument:
        pipeline.addInstrument(instrument)

    # Determine the set of dataset type names affected by source injection
    injected_types = {dataset_type_name}
    precursor_injection_task_labels = set()
    # Loop over all tasks in the pipeline.
    for taskDef in pipeline.toExpandedPipeline():
        conns = taskDef.connections
        input_types = _get_dataset_type_names(conns, conns.inputs)
        output_types = _get_dataset_type_names(conns, conns.outputs)
        if dataset_type_name in output_types:
            precursor_injection_task_labels.add(taskDef.label)
        # If the task has any injected dataset type names as inputs, add all of
        # its outputs to the set of injected types.
        if len(input_types & injected_types) > 0:
            injected_types |= output_types
            # Add the injection prefix to all affected dataset type names.
            for field in conns.inputs | conns.outputs:
                if (conn_type := getattr(conns, field).name) in injected_types:
                    pipeline.addConfigOverride(taskDef.label, "connections." + field, prefix + conn_type)

    # Merge the injection pipeline to the modified pipeline, if provided.
    if injection_pipeline:
        pipeline2 = Pipeline.fromFile(injection_pipeline)
        if len(pipeline2) != 1:
            raise RuntimeError(
                f"The injection pipeline contains {len(pipeline2)} tasks; only one task is allowed."
            )
        pipeline.mergePipeline(pipeline2)
        # Loop over all injection tasks and modify the connection names.
        for injection_taskDef in pipeline2.toExpandedPipeline():
            conns = injection_taskDef.connections
            pipeline.addConfigOverride(
                injection_taskDef.label, "connections.input_exposure", dataset_type_name
            )
            pipeline.addConfigOverride(
                injection_taskDef.label, "connections.output_exposure", prefix + dataset_type_name
            )
            # Optionally update subsets to include the injection task.
            if not exclude_subsets:
                for label in precursor_injection_task_labels:
                    precursor_subsets = pipeline.findSubsetsWithLabel(label)
                    for subset in precursor_subsets:
                        pipeline.addLabelToSubset(subset, injection_taskDef.label)
    return pipeline
