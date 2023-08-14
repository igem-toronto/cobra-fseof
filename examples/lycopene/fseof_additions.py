from typing import Dict, Optional

import cobra
import cobra.manipulation
from cobra import Reaction


def add_single_gene_reaction_pair_lyc(
    model: cobra.Model,
    gene_id: str,
    reaction_id: str,
    reaction_name: str,
    metabolites: Dict[str, float],
    gene_name: Optional[str] = None,
) -> None:
    assert not model.genes.query(lambda k: k == gene_id, attribute="id")
    assert not model.reactions.query(lambda k: k == reaction_id, attribute="id")

    Reaction = cobra.Reaction
    rxn = Reaction(id=reaction_id)

    if gene_name is None:
        gene_name = gene_id
    gene = cobra.Gene(gene_id, name=gene_name)

    model.add_reactions([rxn])
    model.genes.add(gene)

    rxn.name = reaction_name
    rxn.bounds = (-1000, 1000)  # TODO: set to something? Assumes that reaction is reversible.
    rxn.add_metabolites(metabolites)
    rxn.gene_reaction_rule = gene_id


def add_crtE(model: cobra.Model) -> None:
    # References: Choi et. al 2010 supplementary file 3

    ggpp = cobra.Metabolite(
        id="ggpp",
        formula="C20H33O7P2",
        name="geranylgeranyl diphosphate",
        charge=-3,
        compartment="c",
    )
    model.add_metabolites([ggpp])

    phyto = cobra.Metabolite(
        id="phyto",
        formula="C40H64",
        name="phytoene",
        charge=0,
        compartment="c",
    )
    model.add_metabolites([phyto])

    lyco = cobra.Metabolite(
        id="lyco",
        formula="C40H56",
        name="lycopene",
        charge=0,
        compartment="c",
    )
    model.add_metabolites([lyco])

    add_single_gene_reaction_pair_lyc(
        model=model,
        gene_id="crtE",
        reaction_id="ZCRTE",
        reaction_name="Synthesis of geranylgeranyl pyrophosphate",
        metabolites={
            "ipdp_c": -1.0,
            "frdp_c": -1.0,
            "ggpp": 1.0,
            "ppi_c": 1.0
        },
    )


def add_crtB(model: cobra.Model) -> None:
    # References: Choi et. al 2010 supplementary file 3

    add_single_gene_reaction_pair_lyc(
        model=model,
        gene_id="crtB",
        reaction_id="ZCRTB",
        reaction_name="Synthesis of phytoene",
        metabolites={
            "ggpp": -2.0,
            "phyto": 1.0,
            "ppi_c": 1.0
        },
    )


def add_crtI(model: cobra.Model) -> None:
    # References: Choi et. al 2010 supplementary file 3

    add_single_gene_reaction_pair_lyc(
        model=model,
        gene_id="crtI",
        reaction_id="ZCRTI",
        reaction_name="synthesis of lycopene from phytoene (dehydrogenation rxn)",
        metabolites={
            "phyto": -1.0,
            "fad_c": -8.0,
            "lyco": 1.0,
            "fadh2_c": 8.0
        },
    )

def add_lyco_dem(model: cobra.Model) -> None:
    # References: Choi et. al 2010 supplementary file 3,
    # https://cnls.lanl.gov/external/qbio2018/Slides/FBA%202/qBio-FBA-lab-slides.pdf (slide 21)

    lyco_dem = Reaction('LYCOdem')
    model.add_reactions([lyco_dem])
    lyco_dem.name = 'Lycopene demand rxn'
    lyco_dem.lower_bound = 0
    lyco_dem.upper_bound = 1000
    lyco_dem.add_metabolites({"lyco": -1.0})










