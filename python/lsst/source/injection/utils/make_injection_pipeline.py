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

import logging

from lsst.analysis.tools.interfaces import AnalysisPipelineTask
from lsst.pipe.base import LabelSpecifier, Pipeline
from lsst.pipe.tasks.calibrate import CalibrateTask
from lsst.pipe.tasks.calibrateImage import CalibrateImageTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask
from lsst.pipe.tasks.multiBand import MeasureMergedCoaddSourcesTask


def _get_dataset_type_names(conns, fields):
    """Return the name of a connection's dataset type."""
    dataset_type_names = set()
    for field in fields:
        dataset_type_names.add(getattr(conns, field).name)
    return dataset_type_names


def _parse_config_override(config_override: str) -> tuple[str, str, str]:
    """Parse a config override string into a label, a key and a value.

    Parameters
    ----------
    config_override : `str`
        Config override string to parse.

    Returns
    -------
    label : `str`
        Label to override.
    key : `str`
        Key to override.
    value : `str`
        Value to override.

    Raises
    ------
    TypeError
        If the config override string cannot be parsed.
    """
    try:
        label, keyvalue = config_override.split(":", 1)
    except ValueError:
        raise TypeError(
            f"Unrecognized syntax for option 'config': '{config_override}' (does not match pattern "
            "(?P<label>.+):(?P<value>.+=.+))"
        ) from None
    try:
        key, value = keyvalue.split("=", 1)
    except ValueError as e:
        raise TypeError(
            f"Could not parse key-value pair '{config_override}' using separator '=', with multiple values "
            f"not allowed: {e}"
        ) from None
    return label, key, value


def make_injection_pipeline(
    dataset_type_name: str,
    reference_pipeline: Pipeline | str,
    injection_pipeline: Pipeline | str | None = None,
    exclude_subsets: bool = False,
    excluded_tasks: set[str]
    | str = {
        "jointcal",
        "gbdesAstrometricFit",
        "fgcmBuildFromIsolatedStars",
        "fgcmFitCycle",
        "fgcmOutputProducts",
    },
    prefix: str = "injected_",
    instrument: str | None = None,
    config: str | list[str] | None = None,
    log_level: int = logging.INFO,
) -> Pipeline:
    """Make an expanded source injection pipeline.

    This function takes a reference pipeline definition file in YAML format and
    prefixes all post-injection dataset type names with the injected prefix. If
    an optional injection pipeline definition YAML file is also provided, the
    injection task will be merged into the pipeline.

    Unless explicitly excluded, all subsets from the reference pipeline
    containing the task which generates the injection dataset type will also be
    updated to include the injection task. A series of new injected subsets
    will also be created. These new subsets are copies of existent subsets, but
    containing only the tasks which are affected by source injection. New
    injected subsets will be the original subset name with the prefix
    'injected_' prepended.

    Parameters
    ----------
    dataset_type_name : `str`
        Name of the dataset type being injected into.
    reference_pipeline : Pipeline | `str`
        Location of a reference pipeline definition YAML file.
    injection_pipeline : Pipeline | `str`, optional
        Location of an injection pipeline definition YAML file stub. If not
        provided, an attempt to infer the injection pipeline will be made based
        on the injected dataset type name.
    exclude_subsets : `bool`, optional
        If True, do not update pipeline subsets to include the injection task.
    excluded_tasks : `set` [`str`] | `str`
        Set or comma-separated string of task labels to exclude from the
        injection pipeline.
    prefix : `str`, optional
        Prefix to prepend to each affected post-injection dataset type name.
    instrument : `str`, optional
        Add instrument overrides. Must be a fully qualified class name.
    config : `str` | `list` [`str`], optional
        Config override for a task, in the format 'label:key=value'.
    log_level : `int`, optional
        The log level to use for logging.

    Returns
    -------
    pipeline : `lsst.pipe.base.Pipeline`
        An expanded source injection pipeline.
    """
    # Instantiate logger.
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # Load the pipeline and apply config overrides, if supplied.
    if isinstance(reference_pipeline, str):
        pipeline = Pipeline.fromFile(reference_pipeline)
    else:
        pipeline = reference_pipeline
    if config:
        if isinstance(config, str):
            config = [config]
        for conf in config:
            config_label, config_key, config_value = _parse_config_override(conf)
            pipeline.addConfigOverride(config_label, config_key, config_value)

    # Add an instrument override, if provided.
    if instrument:
        pipeline.addInstrument(instrument)

    # Remove all tasks which are not to be included in the injection pipeline.
    if isinstance(excluded_tasks, str):
        excluded_tasks = set(excluded_tasks.split(","))
    all_tasks = {taskDef.label for taskDef in pipeline.toExpandedPipeline()}
    preserved_tasks = all_tasks - excluded_tasks
    label_specifier = LabelSpecifier(labels=preserved_tasks)
    # EDIT mode removes tasks from parent subsets but keeps the subset itself.
    pipeline = pipeline.subsetFromLabels(label_specifier, pipeline.PipelineSubsetCtrl.EDIT)
    if len(not_found_tasks := excluded_tasks - all_tasks) > 0:
        grammar = "Task" if len(not_found_tasks) == 1 else "Tasks"
        logger.warning(
            "%s marked for exclusion not found in the reference pipeline: %s.",
            grammar,
            ", ".join(sorted(not_found_tasks)),
        )

    # Determine the set of dataset type names affected by source injection.
    injected_tasks = set()
    all_connection_type_names = set()
    injected_types = {dataset_type_name}
    precursor_injection_task_labels = set()
    # Loop over all tasks in the pipeline.
    for taskDef in pipeline.toExpandedPipeline():
        # Add override for Analysis Tools taskDefs. Connections in Analysis
        # Tools are dynamically assigned, and so are not able to be modified in
        # the same way as a static connection. Instead, we add a config
        # override here to the connections.outputName field. This field is
        # prepended to all Analysis Tools connections, and so will prepend the
        # injection prefix to all plot/metric outputs. Further processing of
        # this taskDef will be skipped thereafter.
        if issubclass(taskDef.taskClass, AnalysisPipelineTask):
            pipeline.addConfigOverride(
                taskDef.label, "connections.outputName", prefix + taskDef.config.connections.outputName
            )
            continue

        def configure_measurement(name):
            """Configure this task's SingleFrameMeasurement subtask (called
            ``name`` in the task config) to include INJECTED and INJECTED_CORE
            in the catalog pixel flags.
            """
            pipeline.addConfigPython(
                taskDef.label,
                f"config.{name}.plugins['base_PixelFlags'].masksFpAnywhere.extend(['INJECTED'])",
            )
            pipeline.addConfigPython(
                taskDef.label,
                f"config.{name}.plugins['base_PixelFlags'].masksFpCenter.extend(['INJECTED_CORE'])",
            )

        def configure_star_selector(name):
            """Configure a star selector in this task (where ``name`` is the
            name of the bad flag list in a star selector task or subtask) to
            add the injected_coreCenter pixel flag to the list of bad flags.
            """
            pipeline.addConfigPython(
                taskDef.label,
                f"config.{name}.extend(['base_PixelFlags_flag_injected_coreCenter'])",
            )

        # Add injection flag configs to relevant tasks.
        # TODO: do this for coadds
        injected_flag_tasks = (CharacterizeImageTask, CalibrateTask, MeasureMergedCoaddSourcesTask)
        if issubclass(taskDef.taskClass, injected_flag_tasks):
            configure_measurement("measurement")
        if issubclass(taskDef.taskClass, CharacterizeImageTask):
            configure_star_selector("measurePsf.starSelector['objectSize'].badFlags")
            configure_star_selector("measureApCorr.sourceSelector['science'].flags.bad")
        if issubclass(taskDef.taskClass, CalibrateTask):
            configure_star_selector("astrometry.sourceSelector['science'].flags.bad")
        if issubclass(taskDef.taskClass, CalibrateImageTask):
            configure_measurement("psf_source_measurement")
            configure_measurement("star_measurement")
            configure_star_selector("psf_measure_psf.starSelector['objectSize'].badFlags")
            configure_star_selector("measure_aperture_correction.sourceSelector['science'].flags.bad")
            configure_star_selector("star_selector['science'].flags.bad")
            injected_core_flag = "base_PixelFlags_flag_injected_coreCenter"
            pipeline.addConfigPython(
                taskDef.label, f"config.star_selector['science'].flags.bad.extend([{injected_core_flag}])"
            )

        conns = taskDef.connections
        input_types = _get_dataset_type_names(conns, conns.initInputs | conns.inputs)
        output_types = _get_dataset_type_names(conns, conns.initOutputs | conns.outputs)
        all_connection_type_names |= input_types | output_types
        # Identify the precursor task: allows appending inject task to subset.
        if dataset_type_name in output_types:
            precursor_injection_task_labels.add(taskDef.label)
        # If the task has any injected dataset type names as inputs, add the
        # task to a set of tasks touched by injection, and add all of the
        # outputs of this task to the set of injected types.
        if len(input_types & injected_types) > 0:
            injected_tasks |= {taskDef.label}
            injected_types |= output_types
            # Add the injection prefix to all affected dataset type names.
            for field in conns.initInputs | conns.inputs | conns.initOutputs | conns.outputs:
                if hasattr(taskDef.config.connections.ConnectionsClass, field):
                    # If the connection type is not dynamic, modify as usual.
                    if (conn_type := getattr(conns, field).name) in injected_types:
                        pipeline.addConfigOverride(taskDef.label, "connections." + field, prefix + conn_type)
                else:
                    # Add log warning if the connection type is dynamic.
                    logger.warning(
                        "Dynamic connection %s in task %s is not supported here. This connection will "
                        "neither be modified nor merged into the output injection pipeline.",
                        field,
                        taskDef.label,
                    )
    # Raise if the injected dataset type does not exist in the pipeline.
    if dataset_type_name not in all_connection_type_names:
        raise RuntimeError(
            f"Dataset type '{dataset_type_name}' not found in the reference pipeline; "
            "no connection type edits to be made."
        )

    # Attempt to infer the injection pipeline from the dataset type name.
    if not injection_pipeline:
        match dataset_type_name:
            case "postISRCCD":
                injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_exposure.yaml"
            case "icExp" | "calexp":
                injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_visit.yaml"
            case "deepCoadd" | "deepCoadd_calexp" | "goodSeeingCoadd":
                injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_coadd.yaml"
            case _:
                # Print a warning rather than a raise, as the user may wish to
                # edit connection names without merging an injection pipeline.
                logger.warning(
                    "Unable to infer injection pipeline stub from dataset type name '%s' and none was "
                    "provided. No injection pipeline will be merged into the output pipeline.",
                    dataset_type_name,
                )
        if injection_pipeline:
            logger.info(
                "Injected dataset type '%s' used to infer injection pipeline: %s",
                dataset_type_name,
                injection_pipeline,
            )

    # Merge the injection pipeline to the modified pipeline, if provided.
    if injection_pipeline:
        if isinstance(injection_pipeline, str):
            injection_pipeline = Pipeline.fromFile(injection_pipeline)
        if len(injection_pipeline) != 1:
            raise RuntimeError(
                f"The injection pipeline contains {len(injection_pipeline)} tasks; only 1 task is allowed."
            )
        pipeline.mergePipeline(injection_pipeline)
        # Loop over all injection tasks and modify the connection names.
        for injection_taskDef in injection_pipeline.toExpandedPipeline():
            injected_tasks |= {injection_taskDef.label}
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

    # Create injected subsets.
    injected_label_specifier = LabelSpecifier(labels=injected_tasks)
    injected_pipeline = pipeline.subsetFromLabels(injected_label_specifier, pipeline.PipelineSubsetCtrl.EDIT)
    injected_subset_labels = set()
    for injected_subset in injected_pipeline.subsets.keys():
        injected_subset_label = "injected_" + injected_subset
        injected_subset_description = (
            "All tasks from the '" + injected_subset + "' subset impacted by source injection."
        )
        if len(injected_subset_tasks := injected_pipeline.subsets[injected_subset]) > 0:
            injected_subset_labels |= {injected_subset_label}
            pipeline.addLabeledSubset(
                injected_subset_label, injected_subset_description, injected_subset_tasks
            )

    grammar1 = "task" if len(pipeline) == 1 else "tasks"
    grammar2 = "subset" if len(injected_subset_labels) == 1 else "subsets"
    logger.info(
        "Made an injection pipeline containing %d %s and %d new injected %s.",
        len(pipeline),
        grammar1,
        len(injected_subset_labels),
        grammar2,
    )
    return pipeline
