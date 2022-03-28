import argparse
import itertools
import logging
import re
from pathlib import Path

import ihm
import ihm.dataset
import ihm.dumper
import ihm.location
import ihm.model
import ihm.protocol
import ihm.representation
import ihm.restraint

from haddock2mmcif.modules.docking import DockingModel
from haddock2mmcif.modules.pdb import PDB
from haddock2mmcif.modules.restraints import AmbigRestraint, UnambigRestraint

log = logging.getLogger("log")
log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(
    " %(asctime)s %(module)s:%(lineno)d %(levelname)s - %(message)s"
)
ch.setFormatter(formatter)
log.addHandler(ch)

PARAM_REGEX = r"{===>}\s(\w+)=(\d.*);"


def rank_clusters(cluster_out, file_list):
    """Rank the clusters based on their combined score."""
    score_dic = {}
    score_regex = r"{\s(-?\d*.?\d*)\s}"
    with open(file_list, "r") as f_fh:
        for model_idx, line in enumerate(f_fh.readlines()):
            match = re.search(score_regex, line)
            if match:
                score = float(match.group(1))
                score_dic[model_idx] = score

    clt_dic = {}
    cluster_line_regex = r"Cluster\s(\d)\s->\s\d+\s(.*)"
    with open(cluster_out, "r") as c_fh:
        for line in c_fh.readlines():
            match = re.search(cluster_line_regex, line)
            if match:
                cluster_name = int(match.group(1))
                cluster_elements = match.group(2).split()
                cluster_elements = map(int, cluster_elements)

                cluster_elements = list(cluster_elements)
                top4 = cluster_elements[:4]

                # avg_total = (sum(score_dic[e] for e in cluster_elements) /
                #              len(cluster_elements))

                avg_top4 = sum(score_dic[e] for e in top4) / len(top4)

                clt_dic[cluster_name] = avg_top4

    ranking_dic = {}
    for ranking, cluster_name in enumerate(sorted(clt_dic, key=clt_dic.get)):
        ranking_dic[ranking + 1] = cluster_name

    return ranking_dic


def get_final_models(path):
    """Get the final clusterN_N.pdb stuctures."""
    cluster_regex = r"cluster(\d)_\d.pdb"
    cluster_dic = {}
    for element in sorted(path.glob("*pdb")):
        match = re.search(cluster_regex, str(element))
        if match:
            cluster_name = int(match.group(1))
            if cluster_name not in cluster_dic:
                cluster_dic[cluster_name] = []
            cluster_dic[cluster_name].append(element)
    return cluster_dic


def list_to_range(lst):
    """Convert a list to ranges."""
    # thanks: https://stackoverflow.com/a/4629241
    for a, b in itertools.groupby(enumerate(lst), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def get_probability(run_cns):
    """Read the run.cns and look for noecv/ncvpart."""
    noecv = False
    ncvpart = 0.0
    with open(run_cns, "r") as fh:
        for line in fh.readlines():
            if "noecv" in line and "true" in line:
                noecv = True
            if "ncvpart" in line:
                match = re.search(PARAM_REGEX, line)
                if match:
                    value = float(match.group(2))
                    ncvpart = float(value)

    if not noecv:
        probability = 1
    elif ncvpart == 0.0:
        # noecv = true but nvcpart not defined
        #  handle this here
        pass
    else:
        probability = 1 / ncvpart

    return probability


def get_flcut(run_cns):
    """Retrieve the flcut parameter."""
    cutoff = 5.0
    with open(run_cns, "r") as fh:
        for line in fh.readlines():
            if "flcut" in line:
                match = re.search(PARAM_REGEX, line)
                if match:
                    cutoff = float(match.group(2))

    return cutoff


def main():

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("rundir", type=str, help="")
    parser.add_argument("--output", type=str, help="")

    args = parser.parse_args()

    rundir = Path(args.rundir)
    log.info(f"Input run directory: {rundir}")
    run_cns = Path(rundir, "run.cns")

    # ==============================================================
    # Initialize the system
    log.info("Initializing System")
    system = ihm.System()

    # ==============================================================
    # Get contents of the simulation and create entities + asym units

    # the complex_1.pdb has all the chains and it has been processed
    #  by haddock, use this one to define the entities
    log.info("Creating Asymetric Units")
    begin_regex = r"(complex_1.pdb)"
    begin_dir = Path(rundir, "begin")
    entity_list = []
    asym_dic = {}
    for element in sorted(begin_dir.glob("*")):
        match = re.search(begin_regex, str(element))
        if match:
            log.info(f"Reading {element}")
            pdb = PDB(element)
            pdb.load()

            for chainID in pdb.seq_dic:
                # create an entitity
                seq = pdb.seq_dic[chainID]
                log.info(f"Creating entity based on chain {chainID}")
                entity = ihm.Entity(seq, description=f"Chain {chainID}")
                entity_list.append(entity)

                # create the assymetric unit
                mapping = pdb.map_dic[chainID]
                asym = ihm.AsymUnit(
                    entity, auth_seq_id_map=mapping, details=f"Subunit {chainID}"
                )
                asym_dic[chainID] = asym

    # Add them to the system
    log.info("Adding Asymetric Units to the System")
    system.entities.extend(entity_list)
    system.asym_units.extend(list(asym_dic.values()))

    # ==============================================================
    # Organize into an assembly
    log.info("Organizing Asymetric Units into Modeled Assembly")
    modeled_assembly = ihm.Assembly(list(asym_dic.values()), name="Modeled assembly")

    # ==============================================================
    # Add the protocol
    log.info("Defining the protocol")
    protocol = ihm.protocol.Protocol(name="HADDOCK")

    # ==============================================================
    # Rank the clusters since the cluster_name is not its ranking
    water_dir = Path(rundir, "structures", "it1", "water")
    water_analysis_dir = Path(water_dir, "analysis")
    cluster_out = Path(water_analysis_dir, "cluster.out")
    file_list = Path(water_dir, "file.list")

    log.info(f"Ranking the clusters from {cluster_out} based on {file_list}")
    cluster_ranking = rank_clusters(cluster_out, file_list)

    # ==============================================================
    # Generate the models based on the clusters
    clustered_structures = get_final_models(rundir)

    log.info(f"Getting the interface cutoff from {run_cns}")
    interface_cutoff = get_flcut(run_cns)

    group_list = []
    for ranking in cluster_ranking:
        cluster_name = cluster_ranking[ranking]
        model_list = []
        for structure in clustered_structures[cluster_name]:
            log.info(f"Processing {structure.name}")
            cluster_pdb = PDB(structure)
            cluster_pdb.load()

            model_id = int(structure.stem.split("_")[1])

            # each model has 2 representations
            #  one is the whole structure as rigid
            #  second is its interface as flexible
            cluster_pdb.get_interface(cutoff=interface_cutoff)

            rep_list = []
            for chainID in cluster_pdb.interface_dic:
                asym = asym_dic[chainID]

                interface_reslist = cluster_pdb.interface_dic[chainID]

                rigid_rep = ihm.representation.AtomicSegment(asym, rigid=True)
                rep_list.append(rigid_rep)

                for elements in list_to_range(interface_reslist):
                    start, end = elements
                    rng = asym(start, end)
                    flex_rep = ihm.representation.AtomicSegment(rng, rigid=False)
                    rep_list.append(flex_rep)

            rep = ihm.representation.Representation(rep_list)

            model = DockingModel(
                assembly=modeled_assembly,
                protocol=protocol,
                representation=rep,
                name=f"model {model_id}",
                assymetric_dic=asym_dic,
                atom_list=cluster_pdb.atom_list,
            )

            model_list.append(model)

        # one group per cluster
        log.info(f"Finalizing group cluster Rank: {ranking} Number:{cluster_name}")
        _group_name = f"Cluster {ranking} (#{cluster_name})"
        model_group = ihm.model.ModelGroup(model_list, name=_group_name)
        group_list.append(model_group)

    # ==============================================================
    # Groups are then placed into states, which can in turn be grouped.
    log.info("Adding the groups to a state")
    state = ihm.model.State(group_list)
    system.state_groups.append(ihm.model.StateGroup([state]))

    # ==============================================================
    # Add the ambiguous restraints
    log.info("Adding Ambiguous Restraints")
    restraint_dir = Path(rundir, "data", "distances")

    ambig_tbl_f = Path(restraint_dir, "ambig.tbl")

    log.info(f"Reading {ambig_tbl_f}")
    ambig = AmbigRestraint(ambig_tbl_f)
    ambig.load()

    loc = ihm.location.InputFileLocation(str(ambig_tbl_f))
    amig_dataset = ihm.dataset.Dataset(loc)
    log.info(f"Getting probability from {run_cns}")
    prob = get_probability(run_cns)
    for i, active in enumerate(ambig.tbl_dic):
        active_res, active_segid = active
        active_asym = asym_dic[active_segid]

        active_rng = active_asym(active_res, active_res)

        passive_l = ambig.tbl_dic[active]
        passive_ranges = []
        for element in passive_l:
            passive_res, passive_segid = element
            passive_asym = asym_dic[passive_segid]

            passive_rng = passive_asym(passive_res, passive_res)
            passive_ranges.append(passive_rng)

        active = ihm.restraint.ResidueFeature(
            [active_rng], details=f"Ambig AIR {i+1} Active"
        )
        passive = ihm.restraint.ResidueFeature(
            passive_ranges, details=f"Ambig AIR {i+1} Passive"
        )

        dist = ihm.restraint.UpperBoundDistanceRestraint(2.0)

        restraint = ihm.restraint.DerivedDistanceRestraint(
            dataset=amig_dataset,
            feature1=active,
            feature2=passive,
            distance=dist,
            probability=prob,
        )

        log.info(f"Adding Restraint {i+1} to System")
        # system.orphan_features.append(active)
        # system.orphan_features.append(passive)
        system.restraints.append(restraint)

    # ==============================================================
    # Add the unambiguous restraints
    log.info("Adding unambiguous restraints")
    unambig_tbl_f = Path(restraint_dir, "unambig.tbl")
    unambig = UnambigRestraint(unambig_tbl_f)
    unambig.load()

    loc = ihm.location.InputFileLocation(str(unambig_tbl_f))
    unamig_dataset = ihm.dataset.Dataset(loc)

    for i, element in enumerate(unambig.tbl_list):
        res_i, segid_i, res_j, segid_j, distance, lower_bound, upper_bound = element
        asym_i = asym_dic[segid_i]
        asym_j = asym_dic[segid_j]

        rng_i = asym_i(res_i, res_i)
        rng_j = asym_j(res_j, res_j)

        rest_i = ihm.restraint.ResidueFeature([rng_i], details=f"Unambig AIR {i+1}_i")
        rest_j = ihm.restraint.ResidueFeature([rng_j], details=f"Unambig AIR {i+1}_j")

        lower = distance - lower_bound
        upper = distance + upper_bound

        distance = ihm.restraint.LowerUpperBoundDistanceRestraint(lower, upper)

        restraint = ihm.restraint.DerivedDistanceRestraint(
            dataset=unamig_dataset,
            feature1=rest_i,
            feature2=rest_j,
            distance=distance,
            probability=prob,
        )

        system.restraints.append(restraint)

    # ==============================================================
    # System is complete, write it to an mmCIF file:
    if args.output:
        output_fname = args.output
    else:
        output_fname = "output.cif"

    log.info(f"Dumping to {output_fname}")

    with open(f"{output_fname}.cif", "w") as fh:
        ihm.dumper.write(fh, [system])


if __name__ == "__main__":
    main()
