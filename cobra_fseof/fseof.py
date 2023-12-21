import collections
from typing import Any, Callable, List, Literal, Optional, Union

import cobra
import numpy as np
import pandas as pd
import tqdm
from cobra.flux_analysis.variability import flux_variability_analysis

ReactionLike = Union[cobra.Reaction, str]
FSEOFResults = collections.namedtuple("FSEOFResults", ["step", "scan", "fva"])


def as_id(obj: Any) -> str:
    return obj if isinstance(obj, str) else obj.id


def default_reactions(model: cobra.Model, exclude: List[str]) -> List[cobra.Reaction]:
    return [
        rxn.id for rxn in model.reactions
        if (rxn.id not in exclude) and (not rxn.boundary) and rxn.functional
    ]


def default_optimizer(model: cobra.Model) -> cobra.Solution:
    return model.optimize()


def fseof(
    model: cobra.Model,
    num_steps: int,
    enforced: ReactionLike,
    objective: ReactionLike,
    enforced_direction: Literal["min", "max"] = "max",
    enforced_frac_opt: float = 0.9,
    reactions: Optional[List[ReactionLike]] = None,
    optimizer: Optional[Callable[[cobra.Model], cobra.Solution]] = None,
) -> FSEOFResults:
    """Runs the FSEOF algorithm [1].

    Args:
        model: the metabolic model on which to run FSEOF.
        num_steps: the number of steps to interpolate between the initial
            and highest enforced fluxes.
        enforced: the enforced (e.g., product) reaction.
        objective: the base objective (e.g., biomass) reaction.
        enforced_direction: the direction in which FSEOF
            optimizes the enforced reaction.
        enforced_frac_opt: the fraction of the theoretical optimal enforced flux
            to use as the upper interpolation boundary.
        reactions: the space of reactions that FSEOF considers when searching for
            targets. If None, this is set to all reactions in the input model.
        optimizer: the FBA optimizer that is used. If None, this is set to
            the regular COBRApy optimizer.

    Returns:
        An named tuple, FSEOFResults, of pandas DataFrames whose fields are:
        - step: Each row summarizes an interpolation step during FSEOF.
            The index is the step ID and the enforced and objective fluxes are
            stored as columns.
        - scan: Each row corresponds to the flux data of a reaction during FSOEF.
            The index is the reaction ID and the columns are: "target"
            (a boolean for whether the reaction is considered a target),
            "vmin" and "vmax" from the FSEOF algorithm, and "flux<i>" (the
            flux of the reaction in the i-th interpolation step).
        - fva: After FSEOF, Choi et al. [1] use FVA to further refine the
            targets. Each row summarizes the FVA data of a target reaction
            identified in the scan DataFrame. The index is the reaction ID
            and the columns are: "rank" (the score of the target, higher meaning
            better), and "min<i>" and "max<i>" (the flux bounds of the reaction
            identified by FVA in the i-th interpolation step).

    References:
        [1] https://doi.org/10.1128/AEM.00115-10
    """

    model = model.copy()

    # Deal with default kwargs
    if reactions is None:
        reactions = default_reactions(model, [enforced, objective])
    if optimizer is None:
        optimizer = default_optimizer

    # Standardize all arguments into a single data type
    enforced = as_id(enforced)
    objective = as_id(objective)
    reactions = list(map(as_id, reactions))

    # Find scanning bounds
    with model:
        model.objective = objective
        enforced_base = optimizer(model).fluxes[enforced]
    with model:
        model.objective = enforced
        model.objective_direction = enforced_direction
        enforced_opt = optimizer(model).fluxes[enforced]

    # FSEOF fails here since there is nothing to scan over
    if enforced_opt == enforced_base:
        raise ValueError()

    # Scan across enforced fluxes
    step_data = []
    scan_solutions = []
    scan_fluxes = np.linspace(enforced_base, enforced_frac_opt * enforced_opt, num=num_steps)

    for step, v in enumerate(tqdm.tqdm(scan_fluxes, desc="FSEOF; Scanning")):
        with model:
            model.reactions.get_by_id(enforced).bounds = (v, v)
            model.objective = objective
            solution = optimizer(model)

            vobj = solution.fluxes[objective]
            step_data.append({"enforced_flux": v, "objective_flux": vobj})
            scan_solutions.append(solution)

    # Compile into tables
    step_data = pd.DataFrame(step_data, index=range(len(step_data)))

    # From paper:
    #   "The targets were selected by identifying fluxes that increased upon the
    #   application of the enforced objective flux without changing the reaction’s
    #   direction."
    scan_data = []
    for rxn_id in reactions:
        fluxes = {f"flux{i}": S[rxn_id] for i, S in enumerate(scan_solutions)}
        vmin = min(fluxes.values())
        vmax = max(fluxes.values())
        target = (vmax * vmin >= 0) and (abs(vmax) > abs(fluxes["flux0"]))
        scan_data.append({"target": target, "vmin": vmin, "vmax": vmax, **fluxes})
    scan_data = pd.DataFrame(scan_data, index=reactions)

    # Use FVA to further rank targets
    candidates = scan_data[scan_data["target"]].index.tolist()

    scan_fva_results = []
    for step, v in enumerate(tqdm.tqdm(scan_fluxes, desc="FSEOF; Running FVA")):
        with model:
            model.reactions.get_by_id(enforced).bounds = (v, v)
            model.objective = objective
            fva_results = flux_variability_analysis(model, candidates)
            scan_fva_results.append(fva_results.to_dict("index"))

    fva_data = []
    for rxn_id in candidates:
        min_fluxes = [d[rxn_id]["minimum"] for d in scan_fva_results]
        max_fluxes = [d[rxn_id]["maximum"] for d in scan_fva_results]

        rank = 0
        if sorted(min_fluxes) == min_fluxes:
            rank += 1
            if min(max_fluxes) < max(min_fluxes):
                rank += 1
                if max(M - m for m, M in zip(min_fluxes, max_fluxes)) < 1e-3:
                    rank += 1

        fva_data.append({
            "rank": rank,
            **{f"min{i}": f for i, f in enumerate(min_fluxes)},
            **{f"max{i}": f for i, f in enumerate(max_fluxes)},
        })
    fva_data = pd.DataFrame(fva_data, index=candidates)
    fva_data = fva_data.sort_values("rank", ascending=False)

    return FSEOFResults(step=step_data, scan=scan_data, fva=fva_data)
