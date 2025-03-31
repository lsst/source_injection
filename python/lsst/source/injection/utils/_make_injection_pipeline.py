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

import itertools
import logging

from lsst.analysis.tools.interfaces import AnalysisPipelineTask
from lsst.pipe.base import LabelSpecifier, Pipeline
from lsst.pipe.base.pipelineIR import ContractError


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
    excluded_tasks: set[str] | str = {
        "jointcal",
        "gbdesAstrometricFit",
        "fgcmBuildFromIsolatedStars",
        "fgcmFitCycle",
        "fgcmOutputProducts",
    },
    prefix: str = "injected_",
    instrument: str | None = None,
    config: str | list[str] | None = None,
    additional_pipelines: list[Pipeline] | list[str] | None = None,
    subset_name: str | None = None,
    subset_description: str = "",
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

    When the injection pipeline is constructed, a check on all existing
    pipeline contracts is performed. If any contracts are violated, they are
    removed from the pipeline. A warning is logged for each contract that is
    removed.

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
    additional_pipelines: `list`[Pipeline] | `list`[`str`]
        Location(s) of additional input pipeline definition YAML file(s).
        Tasks from these additional pipelines will be added to the output
        injection pipeline.
    subset_name: `str`, optional
        All tasks from any additional pipelines will be added to this subset.
        The subset will be created if it does not already exist.
    subset_description: `str`, optional
        The description given to a new subset which holds tasks from additional
        pipelines provided. Note: this argument is ignored if the subset
        already exists.
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
    all_tasks = set(pipeline.task_labels)
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

    # Check for any empty subsets and remove them.
    removed_subsets = set()
    for subset_label, subset_tasks in pipeline.subsets.items():
        if not subset_tasks:
            removed_subsets.add(subset_label)
            pipeline.removeLabeledSubset(subset_label)
    if (removed_subsets_count := len(removed_subsets)) > 0:
        grammar = "subset" if removed_subsets_count == 1 else "subsets"
        logger.warning(
            "Removed %d empty %s from the pipeline: %s.",
            removed_subsets_count,
            grammar,
            ", ".join(sorted(removed_subsets)),
        )

    # Determine the set of dataset type names affected by source injection.
    injected_tasks = set()
    all_connection_type_names = set()
    injected_types = {dataset_type_name}
    precursor_injection_task_labels = set()
    # Loop over all tasks in the pipeline.
    for task_node in pipeline.to_graph().tasks.values():
        # Add override for Analysis Tools task outputs (but not inputs).
        # Connections in Analysis Tools are dynamically assigned, and so are
        # not able to be modified in the same way as a static connection.
        # Instead, we add an override to the connections.outputName field.
        # This field is prepended to all Analysis Tools connections, and so
        # will prepend the injection prefix to all plot/metric outputs.
        if isAnalysisPipelineTask := issubclass(task_node.task_class, AnalysisPipelineTask):
            pipeline.addConfigOverride(
                task_node.label,
                "connections.outputName",
                prefix + task_node.config.connections.outputName,
            )

        input_types = {
            read_edge.parent_dataset_type_name
            for read_edge in itertools.chain(task_node.inputs.values(), task_node.init.inputs.values())
        }
        output_types = {
            write_edge.parent_dataset_type_name
            for write_edge in itertools.chain(task_node.outputs.values(), task_node.init.outputs.values())
        }

        all_connection_type_names |= input_types | output_types
        # Identify the precursor task: allows appending inject task to subset.
        if dataset_type_name in output_types:
            precursor_injection_task_labels.add(task_node.label)
        # If the task has any injected dataset type names as inputs, add the
        # task to a set of tasks touched by injection, and add all of the
        # outputs of this task to the set of injected types.
        if len(input_types & injected_types) > 0:
            injected_tasks |= {task_node.label}
            injected_types |= output_types
            # Add the injection prefix to all affected dataset type names.
            for edge in itertools.chain(
                task_node.init.inputs.values(),
                task_node.inputs.values(),
                task_node.init.outputs.values(),
                task_node.outputs.values(),
            ):
                # Continue if this is an analysis task and edge is an output.
                if isAnalysisPipelineTask and (
                    edge in set(task_node.init.outputs.values()) | set(task_node.outputs.values())
                ):
                    continue
                if hasattr(task_node.config.connections.ConnectionsClass, edge.connection_name):
                    # If the connection type is not dynamic, modify as usual.
                    if edge.parent_dataset_type_name in injected_types:
                        pipeline.addConfigOverride(
                            task_node.label,
                            "connections." + edge.connection_name,
                            prefix + edge.dataset_type_name,
                        )
                else:
                    # Add log warning if the connection type is dynamic.
                    logger.warning(
                        "Dynamic connection %s in task %s is not supported here. This connection will "
                        "neither be modified nor merged into the output injection pipeline.",
                        edge.connection_name,
                        task_node.label,
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
            case "postISRCCD" | "post_isr_image":
                injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_exposure.yaml"
            case "icExp" | "calexp" | "initial_pvi" | "pvi | preliminary_visit_image | visit_image":
                injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_visit.yaml"
            case (
                "deepCoadd"
                | "deepCoadd_calexp"
                | "goodSeeingCoadd"
                | "deep_coadd_predetection"
                | "deep_coadd"
                | "template_coadd"
            ):
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
        for injection_task_label in injection_pipeline.task_labels:
            injected_tasks.add(injection_task_label)
            pipeline.addConfigOverride(injection_task_label, "connections.input_exposure", dataset_type_name)
            pipeline.addConfigOverride(
                injection_task_label, "connections.output_exposure", prefix + dataset_type_name
            )
            pipeline.addConfigOverride(
                injection_task_label, "connections.output_catalog", prefix + dataset_type_name + "_catalog"
            )
            # Optionally update subsets to include the injection task.
            if not exclude_subsets:
                for label in precursor_injection_task_labels:
                    precursor_subsets = pipeline.findSubsetsWithLabel(label)
                    for subset in precursor_subsets:
                        pipeline.addLabelToSubset(subset, injection_task_label)

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

    # Optionally include additional tasks in the injection pipeline.
    if additional_pipelines:
        additional_tasks: set[str] = set()
        # Record all input task labels and merge all input pipelines into the
        # injection pipeline.
        for additional_pipeline in additional_pipelines:
            if isinstance(additional_pipeline, str):
                additional_pipeline = Pipeline.fromFile(additional_pipeline)
            additional_tasks.update(additional_pipeline.task_labels)
            pipeline.mergePipeline(additional_pipeline)

        # Add all tasks to subset_name. If the subset does not exist create it.
        if isinstance(subset_name, str):
            if subset_name in pipeline.subsets.keys():
                for additional_task in additional_tasks:
                    pipeline.addLabelToSubset(subset_name, additional_task)
                    subset_grammar = f"the existing subset {subset_name}"
            else:
                pipeline.addLabeledSubset(subset_name, subset_description, additional_tasks)
                subset_grammar = f"a new subset {subset_name}"

        # Logging info.
        task_grammar = "task" if len(additional_tasks) == 1 else "tasks"
        logger.info(
            "Added %d %s to %s",
            len(additional_tasks),
            task_grammar,
            subset_grammar,
        )

    # Validate contracts, and remove any that are violated
    try:
        _ = pipeline.to_graph()
    except ContractError:
        contracts_initial = pipeline._pipelineIR.contracts
        pipeline._pipelineIR.contracts = []
        contracts_passed = []
        contracts_failed = []
        for contract in contracts_initial:
            pipeline._pipelineIR.contracts = [contract]
            try:
                _ = pipeline.to_graph()
            except ContractError:
                contracts_failed.append(contract)
                continue
            contracts_passed.append(contract)
        pipeline._pipelineIR.contracts = contracts_passed
        if contracts_failed:
            logger.warning(
                "The following contracts were violated and have been removed: \n%s",
                "\n".join([str(contract) for contract in contracts_failed]),
            )

    return pipeline
