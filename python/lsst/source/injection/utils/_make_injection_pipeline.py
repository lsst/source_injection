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
import warnings

from lsst.pipe.base import LabelSpecifier, Pipeline, PipelineGraph
from lsst.pipe.base.pipelineIR import ContractError


def _infer_injection_pipeline(dataset_type_name: str, logger: logging.Logger) -> str | None:
    """Infer the injection pipeline from the dataset type name.

    Parameters
    ----------
    dataset_type_name : `str`
        Name of the dataset type being injected into.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Returns
    -------
    injection_pipeline : `str` | `None`
        Location of an injection pipeline definition YAML file stub, or None if
        no suitable injection pipeline could be inferred.
    """
    injection_pipeline = None
    match dataset_type_name:
        case "postISRCCD" | "post_isr_image":
            injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_exposure.yaml"
        case "icExp" | "calexp" | "initial_pvi" | "pvi" | "preliminary_visit_image" | "visit_image":
            injection_pipeline = "$SOURCE_INJECTION_DIR/pipelines/inject_visit.yaml"
        case (
            "deepCoadd"
            | "deepCoadd_calexp"
            | "goodSeeingCoadd"
            | "deep_coadd_predetection"
            | "deep_coadd"
            | "deep_coadd_cell_predetection"
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
    return injection_pipeline


def _merge_injection_pipeline(
    pipeline: Pipeline,
    injection_pipeline: Pipeline | str | None,
    dataset_type_name: str,
    prefix: str,
) -> None | str:
    """Merge an injection pipeline into an existing pipeline.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to merge the injection pipeline into.
    injection_pipeline : `~lsst.pipe.base.Pipeline` | `str` | `None`
        Injection pipeline to merge, or location of an injection pipeline
        definition YAML file stub. If None, no injection pipeline is merged.
    dataset_type_name : `str`
        Name of the dataset type being injected into.
    prefix : `str`
        Prefix to prepend to each affected post-injection dataset type name.

    Returns
    -------
    injection_task_label : `str` | `None`
        Label of the injection task, if an injection pipeline was merged, or
        None if no injection pipeline was merged.

    Notes
    -----
    This function modifies the input pipeline in place.
    """
    if injection_pipeline is None:
        return None
    if isinstance(injection_pipeline, str):
        injection_pipeline = Pipeline.fromFile(injection_pipeline)
    if len(injection_pipeline) != 1:
        raise RuntimeError(
            f"The injection pipeline contains {len(injection_pipeline)} tasks; only 1 task is allowed."
        )
    pipeline.mergePipeline(injection_pipeline)

    injection_task_label = next(iter(injection_pipeline.task_labels))
    pipeline.addConfigOverride(injection_task_label, "connections.input_exposure", dataset_type_name)
    pipeline.addConfigOverride(
        injection_task_label, "connections.output_exposure", prefix + dataset_type_name
    )
    pipeline.addConfigOverride(
        injection_task_label, "connections.output_catalog", prefix + dataset_type_name + "_catalog"
    )
    return injection_task_label


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


def _configure_injection_pipeline(
    pipeline: Pipeline,
    config: str | list[str],
    logger: logging.Logger,
) -> None:
    """Apply user-supplied config overrides to the pipeline.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to apply config overrides to. Pipeline is modified in place.
    config : `str` | `list` [`str`]
        Config override(s) to apply, in the format 'label:key=value'.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Notes
    -----
    This function modifies the input pipeline in place.
    """
    if isinstance(config, str):
        config = [config]
    for conf in config:
        config_label, config_key, config_value = _parse_config_override(conf)
        try:
            pipeline.addConfigOverride(config_label, config_key, config_value)
        except LookupError:
            logger.debug(
                "Config override '%s' for label '%s' not found in the reference "
                "pipeline, either due to a typo or the label not existing in "
                "the reference pipeline.",
                conf,
                config_label,
            )


def _remove_excluded_tasks(
    pipeline: Pipeline,
    excluded_tasks: set[str] | str,
    logger: logging.Logger,
) -> Pipeline:
    """Remove excluded tasks from the pipeline and any subsets,
    and remove any empty subsets.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to remove tasks from. This pipeline is modified in place.
    excluded_tasks : `set` [`str`] | `str`
        Task labels to exclude from the injection pipeline.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Returns
    -------
    pipeline : `~lsst.pipe.base.Pipeline`
        The input pipeline with excluded tasks and empty subsets removed.
    """
    if isinstance(excluded_tasks, str):
        excluded_tasks = set(excluded_tasks.split(","))
    all_tasks = set(pipeline.task_labels)
    preserved_tasks = all_tasks - excluded_tasks

    preserved_task_labels = LabelSpecifier(labels=preserved_tasks)
    # EDIT mode removes tasks from parent subsets but keeps the subset itself.
    pipeline = pipeline.subsetFromLabels(preserved_task_labels, pipeline.PipelineSubsetCtrl.EDIT)

    if len(found_tasks := excluded_tasks & all_tasks) > 0:
        grammar = "task" if len(found_tasks) == 1 else "tasks"
        logger.info(
            "%d %s excluded from the output pipeline: %s",
            len(found_tasks),
            grammar,
            ", ".join(sorted(found_tasks)),
        )

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

    return pipeline


def _get_pipeline_graph(pipeline: Pipeline, logger: logging.Logger) -> PipelineGraph:
    """Get the pipeline graph, handling any contract errors.

    Pipeline contracts that are violated by any modifications made to the
    pipeline will be removed, with a warning logged for each contract that's
    removed.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to validate contracts for.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Returns
    -------
    pipeline_graph : `~lsst.pipe.base.PipelineGraph`
        The pipeline graph for the input pipeline, with any violated contracts
        removed from the input pipeline.

    Notes
    -----
    This function modifies the input pipeline in place, removing any violated
    contracts.
    """
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r".*formatted like a Pipeline parameter but was not found within the Pipeline.*",
                category=UserWarning,
            )
            pipeline_graph = pipeline.to_graph()
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
        pipeline_graph = pipeline.to_graph()

        if contracts_failed:
            logger.warning(
                "The following contracts were violated and have been removed: \n%s",
                "\n".join([str(contract) for contract in contracts_failed]),
            )
    return pipeline_graph


def _collect_injected_task_labels(
    pipeline_graph: PipelineGraph,
    dataset_type_name: str,
) -> set[str]:
    """Collect tasks downstream of the injection point.

    Parameters
    ----------
    pipeline_graph : `~lsst.pipe.base.PipelineGraph`
        Pipeline graph to inspect.
    dataset_type_name : `str`
        Name of the dataset type being injected into.

    Returns
    -------
    injected_task_labels : `set` [`str`]
        Labels of all tasks that consume the injected dataset type directly or
        indirectly, including the injection task itself if present.
    """
    injected_task_labels = set()

    dataset_type_frontier = {dataset_type_name}
    seen_dataset_types = set(dataset_type_frontier)

    # Note: here we opt to walk the pipeline graph instead of using
    # `pipeline_graph._xgraph.successors`. The `_xgraph` attribute is a private
    # implementation detail and therefore not a guaranteed interface.

    while dataset_type_frontier:
        next_frontier = set()
        for current_dataset_type in dataset_type_frontier:
            for task_node in pipeline_graph.consumers_of(current_dataset_type):
                if task_node.label in injected_task_labels:
                    continue

                injected_task_labels.add(task_node.label)

                output_edges = task_node.iter_all_outputs()

                for edge in output_edges:
                    output_dataset_type = edge.parent_dataset_type_name
                    if output_dataset_type not in seen_dataset_types:
                        seen_dataset_types.add(output_dataset_type)
                        next_frontier.add(output_dataset_type)
        dataset_type_frontier = next_frontier

    return injected_task_labels


def _add_injected_subsets(
    pipeline: Pipeline,
    injected_task_labels: set[str],
    prefix: str,
    logger: logging.Logger,
) -> int:
    """Create injected variants of existing subsets.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to modify in place.
    injected_task_labels : `set` [`str`]
        Labels of tasks downstream of the injection point.
    prefix : `str`
        Prefix to prepend to the subset names.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Returns
    -------
    subset_count : `int`
        Number of injected subsets created.
    """
    if not injected_task_labels:
        return 0

    injected_label_specifier = LabelSpecifier(labels=injected_task_labels)
    injected_pipeline = pipeline.subsetFromLabels(injected_label_specifier, pipeline.PipelineSubsetCtrl.EDIT)

    injected_subset_labels = set()
    for subset_label, subset_tasks in injected_pipeline.subsets.items():
        if not subset_tasks:
            continue
        injected_subset_label = prefix + subset_label
        injected_subset_description = (
            f"All tasks from the '{subset_label}' subset impacted by source injection."
        )

        # If the subset already exists, add any new tasks to the existing
        # subset rather than creating a new subset with the same name.
        # Used in pipelines with more than one injection tasks.
        if injected_subset_label in pipeline.subsets.keys():
            existing_subset_tasks = set(pipeline.subsets[injected_subset_label])
            for injected_subset_task in set(subset_tasks) - existing_subset_tasks:
                pipeline.addLabelToSubset(injected_subset_label, injected_subset_task)
        else:
            pipeline.addLabeledSubset(injected_subset_label, injected_subset_description, subset_tasks)
            injected_subset_labels.add(injected_subset_label)

    return len(injected_subset_labels)


def _reconfigure_injection_pipeline(
    pipeline: Pipeline,
    dataset_type_name: str,
    prefix: str,
    injection_task_label: str | None,
    update_subsets: bool,
    logger: logging.Logger,
) -> None:
    """Reconfigure the injection pipeline by prefixing post-injection dataset
    type names and updating subsets.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to configure. This pipeline is modified in place.
    dataset_type_name : `str`
        Name of the dataset type being injected into.
    prefix : `str`
        Prefix to prepend to each affected post-injection dataset type name.
    injection_task_label : `str` | `None`
        Label of the injection task.
    update_subsets : `bool`
        If True, update pipeline subsets to include the injection task.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Notes
    -----
    This function modifies the input pipeline in place.
    """
    # Use pipeline graph to determine tasks with connections to be modified
    pipeline_graph = _get_pipeline_graph(pipeline, logger)
    injected_task_labels = _collect_injected_task_labels(pipeline_graph, dataset_type_name)
    post_injection_tasks = pipeline_graph.consumers_of(dataset_type_name)
    if len(post_injection_tasks) == 0:
        logger.warning(
            "Dataset type '%s' not found in the reference pipeline; no input connection edits to be made.",
            dataset_type_name,
        )
    if post_injection_tasks:
        post_injection_tasks = [task for task in post_injection_tasks if task.label != injection_task_label]
    else:
        post_injection_tasks = []

    # Loop over each post injection task; prefix input connections only
    for task_node in post_injection_tasks:
        input_edges = task_node.iter_all_inputs()

        for edge in input_edges:
            if hasattr(task_node.config.connections.ConnectionsClass, edge.connection_name):
                if edge.parent_dataset_type_name == dataset_type_name:
                    pipeline.addConfigOverride(
                        task_node.label,
                        "connections." + edge.connection_name,
                        prefix + edge.dataset_type_name,
                    )

    # Update subsets to include the injection task
    if (
        update_subsets
        and injection_task_label is not None
        and (pre_injection_task := pipeline_graph.producer_of(dataset_type_name)) is not None
    ):
        precursor_subsets = pipeline.findSubsetsWithLabel(pre_injection_task.label)
        for subset in precursor_subsets:
            pipeline.addLabelToSubset(subset, injection_task_label)

    injected_subset_count = 0
    if update_subsets:
        injected_subset_count = _add_injected_subsets(pipeline, injected_task_labels, prefix, logger)

    logger.info(
        "Made an injection pipeline containing %d task%s and %d injected subset%s.",
        len(pipeline),
        "" if len(pipeline) == 1 else "s",
        injected_subset_count,
        "" if injected_subset_count == 1 else "s",
    )


def _add_additional_pipelines(
    pipeline: Pipeline,
    additional_pipelines: list[Pipeline] | list[str],
    additional_subset: list[str] | str | None,
    logger: logging.Logger,
) -> None:
    """Add additional pipelines to the injection pipeline, and optionally add
    all additional tasks to existing or new subsets.

    Parameters
    ----------
    pipeline : `~lsst.pipe.base.Pipeline`
        Pipeline to add additional pipelines to. Pipeline is modified in place.
    additional_pipelines : `list` [`~lsst.pipe.base.Pipeline`] | `list` [`str`]
        Additional pipelines to merge, or locations of additional pipeline
        definition YAML file stubs.
    additional_subset : `list` [`str`] | `str` | `None`
        A list of subset definitions in the form
        "subset_name[:subset_description]".
        These subsets will be created if they don't already exist. All tasks
        from the additional pipelines will be added to these subsets.
        If None, additional tasks will not be added to any subsets.
    logger : `~logging.Logger`
        Logger for warning and info messages.

    Notes
    -----
    This function modifies the input pipeline in place.
    """
    # Merge all additional pipelines into the main pipeline
    additional_tasks: set[str] = set()
    for additional_pipeline in additional_pipelines:
        if isinstance(additional_pipeline, str):
            additional_pipeline = Pipeline.fromFile(additional_pipeline)
        additional_tasks.update(additional_pipeline.task_labels)
        pipeline.mergePipeline(additional_pipeline)

    # Add all tasks to subset_name; create the subset if it does not exist
    subset_text = ""
    if additional_subset is not None:
        if not isinstance(additional_subset, list):
            additional_subset = [additional_subset]
        subset_names_old = []
        subset_names_new = []
        for subset in additional_subset:
            # Parse the subset definition
            if ":" in subset:
                subset_name, subset_description = subset.split(":", 1)
            else:
                subset_name = subset
                subset_description = ""
            # Add or create the subset with all additional tasks
            if subset_name in pipeline.subsets:
                subset_names_old.append(subset_name)
                for additional_task in additional_tasks:
                    pipeline.addLabelToSubset(subset_name, additional_task)
            else:
                subset_names_new.append(subset_name)
                pipeline.addLabeledSubset(subset_name, subset_description, additional_tasks)
        if subset_names_old:
            subset_text += f", and to existing subset{'s' if len(subset_names_old) > 1 else ''} "
            subset_text += f"{', '.join(sorted(subset_names_old))}"
        if subset_names_new:
            subset_text += f", and to new subset{'s' if len(subset_names_new) > 1 else ''} "
            subset_text += f"{', '.join(sorted(subset_names_new))}"

    # Revalidate the pipeline graph
    _ = _get_pipeline_graph(pipeline, logger)

    grammar = "task" if len(additional_tasks) == 1 else "tasks"
    logger.info(
        "Added %d %s to the pipeline%s: %s",
        len(additional_tasks),
        grammar,
        subset_text,
        ", ".join(sorted(additional_tasks)),
    )


def make_injection_pipeline(
    dataset_type_name: str,
    reference_pipeline: Pipeline | str,
    injection_pipeline: Pipeline | str | None = None,
    update_subsets: bool = True,
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
    additional_subset: list[str] | str | None = None,
    log_level: int = logging.INFO,
) -> Pipeline:
    """Make an expanded source injection pipeline.

    This function takes a reference pipeline definition file and prefixes all
    immediately post-injection dataset type names with the injected prefix. If
    an optional injection pipeline definition YAML file is also provided, the
    injection task will be merged into the pipeline.

    Unless subset updates are explicitly disabled, all subsets from the
    reference pipeline containing the task which generates the injection
    dataset type will also be updated to include the injection task.

    When the injection pipeline is constructed, a check on all existing
    pipeline contracts is performed. If any contracts are violated, they're
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
    update_subsets : `bool`, optional
        If True, update pipeline subsets to include the injection task.
    excluded_tasks : `set` [`str`] | `str`
        Set of task labels to exclude, or a comma-separated string of labels.
    prefix : `str`, optional
        Prefix to prepend to each affected post-injection dataset type name.
    instrument : `str`, optional
        Add instrument overrides. Must be a fully qualified class name.
    config : `str` | `list` [`str`], optional
        Config override for a task, in the format 'label:key=value'.
    additional_pipelines: `list`[Pipeline] | `list`[`str`], optional
        Additional pipelines to merge into the output pipeline, or their YAML
        file locations. Tasks from these additional pipelines will be added to
        the output injection pipeline.
    additional_subset: `list`[`str`] | `str`, optional
        A list of subset definitions in the form
        "subset_name[:subset_description]".
        These subsets will be created if they don't already exist.
        All tasks from additional_pipelines will be added to these subsets.
    log_level : `int`, optional
        The log level to use for logging.

    Returns
    -------
    pipeline : `lsst.pipe.base.Pipeline`
        An expanded source injection pipeline.
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # Get the main reference pipeline
    if isinstance(reference_pipeline, str):
        pipeline = Pipeline.fromFile(reference_pipeline)
    else:
        pipeline = reference_pipeline

    # Add an instrument override
    if instrument:
        pipeline.addInstrument(instrument)

    # Infer the injection pipeline if not provided, and where possible
    if not injection_pipeline:
        injection_pipeline = _infer_injection_pipeline(
            dataset_type_name,
            logger,
        )

    # Merge the injection pipeline into the main pipeline
    injection_task_label = _merge_injection_pipeline(pipeline, injection_pipeline, dataset_type_name, prefix)

    # Apply all user-supplied config overrides
    if config is not None:
        _configure_injection_pipeline(pipeline, config, logger)

    # Remove excluded tasks from the pipeline, and remove any empty subsets
    pipeline = _remove_excluded_tasks(pipeline, excluded_tasks, logger)

    # Prefix post-injection dataset type name connections and update subsets
    _reconfigure_injection_pipeline(
        pipeline, dataset_type_name, prefix, injection_task_label, update_subsets, logger
    )

    # Optionally include additional tasks in the injection pipeline.
    if additional_pipelines:
        _add_additional_pipelines(pipeline, additional_pipelines, additional_subset, logger)

    return pipeline
