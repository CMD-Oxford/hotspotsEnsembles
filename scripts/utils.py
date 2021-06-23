from pathlib import Path
from ccdc.io import MoleculeWriter, MoleculeReader
from ccdc.protein import Protein
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.wrapper_protoss import *
from hotspots.hs_utilities import Helper
from hotspots.hs_ensembles import EnsembleResult
from hotspots.result import _Scorer
from hotspots.data import common_solvents
from hotspots.grid_extension import Grid, _GridEnsemble
import numpy as np
import pandas as pd
from collections import Counter

def get_clusters_centre_mass(cluster_array, hotspot_map):
    """
    May not be the best place for it. Check if it isn't already in _GridEnsemble
    :param cluster_array:
    :param hotspot_map: Note - this is a numpy array!
    :return:
    """
    coords = {}
    for c in set(cluster_array[cluster_array > 0]):
        arr = hotspot_map * (cluster_array == c)
        coords[c] = _GridEnsemble.get_center_of_mass(arr)
    return coords

def find_bs_residues(prot, list_of_ligands):
    """
    Finds the bindinding site based on residues within 5A of any bound ligand in the ensemble
    :return:
    """
    bs_resi_list = []
    for l in list_of_ligands:
        binding_site = Protein.BindingSiteFromMolecule(protein=prot, molecule=l, distance=5.0)
        for r in binding_site.residues:
            if r not in list(set(bs_resi_list)):
                bs_resi_list.append(r)
    return bs_resi_list

def find_bs_ligands(prot, resi_list, threshold=4.0):
    """
    Finds the ligands in the binding site from the list of binding site residues
    :param prot:
    :param resi_list:
    :return:
    """
    prot.detect_ligand_bonds()
    ccdc_resi = [r for r in prot.residues  for a in resi_list if r.identifier==a]
    binding_site = Protein.BindingSiteFromListOfResidues(prot, ccdc_resi)
    bs_ligands = []
    bs_coords = [a.coordinates for r in binding_site.residues for a in r.atoms]
    # If any atom of the ligand is within H-bonding distance (generous = 4 A) of any of the residues, it is from the binding site

    for lig in prot.ligands:
        # Check that the ligand is not just common solvent/dmso
        if not any([x in lig.identifier[2:5] for x in common_solvents() if len(x)>2]):
            lig_coords = [a.coordinates for a in lig.heavy_atoms]
            dists = np.array([Helper.get_distance(lig_coord, bs_coord) for lig_coord in lig_coords for bs_coord in bs_coords])
            if np.any(dists <= threshold):
                bs_ligands.append(lig)

    return bs_ligands

def get_pdb_resolution(pdb_file_path):
    """
    Given the path to a pdb file, finds the resolution field. TBD whether it works in all
    :param pdb_file_path:
    :return:
    """
    print(pdb_file_path)
    prot_str = Path(pdb_file_path).read_text()
    lines = prot_str.split('\n')
    res_lines = [l for l in lines if 'RESOLUTION.' in l]
    try:
        res_line = res_lines[0]
        fields = [x for x in res_line.strip().split(' ') if x != ""]
        res = float(fields[fields.index('RESOLUTION.') + 1])
    except:
        if len(res_lines) > 1:
            res_line = res_lines[1]
            fields = [x for x in res_line.strip().split(' ') if x != ""]
            res = float(fields[fields.index('RESOLUTION.') + 1])
        else:
            res = 'Cannot find resolution in pdb file'

    return res

def get_subset(path_to_ensemble_data, data_source, mw_threshold=300.0, unique_ligands=True):
    """
    Given the path to an ensemble dataframe, selects all the structures that correspond to the conditions
    :param path_to_ensemble_data: Path
    :param data_source: str, either 'KLIFS' or SIENA
    :param mw_threshold: float
    :param unique_ligands: bool
    :return:
    """
    df = pd.read_csv(str(path_to_ensemble_data))
    if data_source == 'SIENA':
        # Get all fields where the ligand has mw less than the threhold
        df = df[df["Ligand MW"] <= mw_threshold]

        # Get unique ligands with the highest resolution
        if unique_ligands:
            unique_ligands = set(df["Ligand SMILES"].values)
            high_res_idx = []
            for ul in unique_ligands:
                sub_df = df[df["Ligand SMILES"]==ul]
                high_res_idx.append(sub_df['Resolution'].idxmin())
            print(sorted(high_res_idx))
            df = df.loc[sorted(high_res_idx)]

    elif data_source == 'KLIFS':
        # the size of the ligands is already taken into account. 
        df = df[df["Ligand MW"] <= mw_threshold]
        df = df[df['alt'] != 'B']
        print(df.columns)

        # First, get the most common sequence (we assume that's the true one).
        seq_counter = Counter(df['pocket'])
        most_common_list = seq_counter.most_common()
        most_common_seq = most_common_list[0][0]

        # If there is more than one sequence, say what it is and truncate the dataframe
        if len(most_common_list) > 1:
            print("Chosen binding site sequence: \n {} \n Next most commonly occurring: \n {}".format(
                most_common_list[0],
                most_common_list[1]))
            df = df[df['pocket'] == most_common_seq]

        # Remove structures in complex with the same ligands. Keep only the highest resolution one.
        ligs = list(set(df['orthosteric_PDB']))
        idx_list = []
        for lig in ligs:
            # Check that the structure doesn't have allosteric ligands:
            sli = df[(df['orthosteric_PDB'] == lig )& (df['allosteric_PDB']== '-')]
            if len(sli) > 0:
                # get the entry with the highest resolution, or if there is a tie, take the first.
                idx = sli[sli['resolution'] == sli['resolution'].min()].index[0]
                idx_list.append(idx)
        print(idx_list)
        print(df)
        df = df.loc[idx_list]

    else:
        print(f"Unrecognised data source {data_source}. Accepted values: 'SIENA', 'KLIFS'")

    return df

def KLIFS_molecular_weight(path_to_ensemble_data):
    df = pd.read_csv(str(path_to_ensemble_data))
    lig_mws = []
    
    for idx, row in df.iterrows():
        if row['alt'] == ' ':
            pname = "{}_chain{}".format(row['pdb'].strip(), row['chain'].strip())
        else:
            print(row['alt'])
            pname = "{}_alt{}_chain{}".format(row['pdb'].strip(), row['alt'].strip(), row['chain'].strip())
        lig = MoleculeReader(str(Path(path_to_ensemble_data.parent, pname, 'ligand.mol2')))[0]
        lig_mws.append(lig.molecular_weight)
    df['Ligand MW'] = lig_mws
    df.to_csv(path_to_ensemble_data)

def process_KLIFS_pdbs(res_table_path, transformation=None):
    """

    :param res_table_path:
    :param transformation:
    :return:
    """
    df = pd.read_csv(res_table_path)
    out_dir = Path(res_table_path.parent, "hotspot_results")
    if not out_dir.exists(): out_dir.mkdir()
    out_paths = []

    for idx, row in df.iterrows():
        # Find the protein file:
        if row['alt'] == ' ':
            pname = "{}_chain{}".format(row['pdb'].strip(), row['chain'].strip())
        else:
            print(row['alt'])
            pname = "{}_alt{}_chain{}".format(row['pdb'].strip(), row['alt'].strip(), row['chain'].strip())

        complex_prot = Protein.from_file(str(Path(res_table_path.parent, pname, 'complex.mol2')))
        complex_prot.detect_ligand_bonds()
        ligands = complex_prot.ligands
        input_prot = Protein.from_file(str(Path(res_table_path.parent, pname, 'protein_aligned.mol2')))

        if transformation:
            complex_prot.transform(transformation)
            input_prot.transform(transformation)
            for l in ligands:
                l.transform(transformation)

        hs_dir = Path(out_dir, pname)
        if not hs_dir.exists(): hs_dir.mkdir()
        with MoleculeWriter(str(Path(out_dir, pname, 'complex.pdb'))) as wc:
            wc.write(complex_prot)

        out_path = Path(out_dir, pname, 'protein.mol2')
        with MoleculeWriter(str(out_path)) as wp:
            wp.write(input_prot)

        # with MoleculeWriter(str(Path(out_dir, pname, 'ligands.sdf'))) as wl:
        #     for l in ligands:
        #         wl.write(l)
        lig = Path(res_table_path.parent, pname, 'ligand_aligned.mol2').read_text()
        Path(out_dir, pname, 'ligand.mol2').write_text(lig)

        out_paths.append(out_path)

    return out_paths


def process_siena_pdbs(res_table_path, transformation=None):
    """

    :param res_table_path: Path to the results table outputted by SIENA
    :return:
    """
    df = pd.read_csv(res_table_path)
    out_paths = []
    out_dir = Path(res_table_path.parent, "hotspot_results")
    if not out_dir.exists(): out_dir.mkdir()
    chain_superposition = Protein.ChainSuperposition()
    sets = chain_superposition.Settings()
    sets.superposition_atoms = "CALPHA"

    protoss = Protoss() # queries the Protoss RESTful API

    for idx, row in df.iterrows():
        # Find the target protein

        prot_path = Path(row['Filename'])
        print(prot_path)

        # Load it up and process
        siena_prot = Protein.from_file(str(prot_path))


        protoss_res = protoss.add_hydrogens(pdb_code=row['PDB ID'])
        cur_prot = protoss_res.protein

        # Remove any ligands apart from that in the binding site
        #bs_ligands = find_bs_ligands(cur_prot, row['Binding site'][2:-2].strip().split("', '"))
        bs_ligands = find_bs_ligands(siena_prot, row['Binding site'][2:-2].strip().split("', '"))
        bs_ligand_ids = [x.identifier for x in bs_ligands]

        # Remove any chains that are not involved in binding/ are symmetry mates
        chains = list(row['Chains'])

        for ch in cur_prot.chains:
            if ch.identifier not in chains:
                cur_prot.remove_chain(ch.identifier)

        for l in cur_prot.ligands:
            if l.identifier not in bs_ligand_ids:
                cur_prot.remove_ligand(l.identifier)

        # back align to the siena protein
        siena_prot_chains = [c for c in siena_prot.chains if c.identifier in chains]
        (rmsd, s_transformation) = chain_superposition.superpose(siena_prot_chains[0], cur_prot.chains[0])
        report_string = f"Structure {row['ID']}, RMSD to siena original: {rmsd} \n"
        print(report_string)

        if transformation:
            cur_prot.transform(transformation)

        comple = cur_prot.copy()

        # cur_prot.add_hydrogens()  # don't need to do this as already protonated with protoss
        for l in cur_prot.ligands: cur_prot.remove_ligand(l.identifier)
        cur_prot.remove_all_metals()
        cur_prot.remove_all_waters()

        out_path = Path(out_dir, prot_path.stem)
        if not out_path.exists(): out_path.mkdir()

        comp_path = Path(out_path, f"{prot_path.stem}_complex.mol2")
        with MoleculeWriter(str(comp_path)) as wr:
            wr.write(comple)

        ppath = Path(out_path, f"{prot_path.stem}_prepped.mol2")
        with MoleculeWriter(str(ppath)) as wr:
            wr.write(cur_prot)

        lig_path = Path(out_path, f"{prot_path.stem}_ligands.sdf")
        with MoleculeWriter(str(lig_path)) as wr:
            for bsl in bs_ligands:
                wr.write(bsl)

        protoss_report = Path(protoss_res.files['log'][0])
        protoss_report_file = Path(out_path, f"{prot_path.stem}_protoss_log.txt")
        protoss_report_file.write_text(protoss_report.read_text())

        out_paths.append(ppath)

    return out_paths

def shrink_bs_maps(hotspot_paths, ligand_paths, padding=4.0):
    """

    :param hotspot_paths:
    :param ligand_paths:
    :param padding:
    :return:
    """
    loaded_ligs = [x for lp in ligand_paths for x in MoleculeReader(str(lp))]
    coords = np.array([a.coordinates for lig in loaded_ligs for a in lig.heavy_atoms])
    min_coords = np.min(coords, axis=0) - padding
    max_coords = np.max(coords, axis=0) + padding

    shrunk_hs_results = []
    shrunk_hs_result_paths = []

    for hs_path in hotspot_paths:
        with HotspotReader(str(hs_path)) as hs_reader:
            hs_res = hs_reader.read()
        probes = hs_res.super_grids.keys()

        # now to shrink the grids for each probe
        for p in probes:
            hs_res.super_grids[p] = EnsembleResult.shrink_to_binding_site(in_grid=hs_res.super_grids[p],
                                                                          new_origin=min_coords,
                                                                          new_far_corner=max_coords)
        shrunk_hs_results.append(hs_res)

        h_out_dir = Path(hs_path.parent, 'binding_site_maps')
        if not h_out_dir.exists(): h_out_dir.mkdir()
        with HotspotWriter(str(h_out_dir.resolve()), grid_extension=".ccp4", zip_results=True) as writer:
            writer.write(hs_res)

        shrunk_hs_result_paths.append(h_out_dir)

    return shrunk_hs_results, shrunk_hs_result_paths

def score_hotspot_protein(in_tup):
    ind = in_tup[0]
    hs_path = in_tup[1]
    out_dir = in_tup[2]

    with HotspotReader(str(hs_path)) as hs_reader:
        hs_result = hs_reader.read()
    sc = _Scorer(hs_result, hs_result.protein, tolerance=2)
    scored_prot = sc._score_protein_backup(hs_result.protein)
    print('scored protein')
    with MoleculeWriter(str(Path(out_dir, '{}_{}.mol2'.format(str(ind), scored_prot.identifier)))) as mwr:
        mwr.write(scored_prot)
    return scored_prot









