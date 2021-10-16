from pathlib import Path
import pandas as pd
from utils import get_subset, get_clusters_centre_mass, shrink_bs_maps
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from ccdc.io import MoleculeReader, MoleculeWriter
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.hs_ensembles import EnsembleResult
from hotspots.grid_extension import _GridEnsemble, Grid
import numpy as np
import matplotlib.pyplot as plt


def ensemble_cluster_summary(ensemble_res, save_dir):
    clust_ids = []
    clust_size = []
    clust_medians = []
    map_size = []
    clust_centroid = []
    clust_contributing = []
    clust_contributions = []

    for probe in ensemble_res.ensemble_hotspot_result.super_grids.keys():
        # Get the contributing grids for all the clusters
        probe_ens = ensemble_res.grid_ensembles[probe]

        # Convert into numpy array:
        probe_array = _GridEnsemble.array_from_grid(ensemble_res.ensemble_hotspot_result.super_grids[probe])

        # Find the hotspot features using HDBSCAN
        if probe == 'apolar':
            ensemble_probe_clusters = _GridEnsemble.HDBSCAN_cluster(probe_array,
                                                                    min_cluster_size=27,
                                                                    allow_single_cluster=True)
        else:
            ensemble_probe_clusters = _GridEnsemble.HDBSCAN_cluster(probe_array,
                                                                    min_cluster_size=7,
                                                                    allow_single_cluster=True)
        contribs = probe_ens.get_contributing_maps(ensemble_probe_clusters)

        probe_clust_grid = ensemble_res.grid_ensembles[probe].as_grid(ensemble_probe_clusters)
        probe_clust_grid.write(str(Path(save_dir, f"{probe}_ensemble_clusters_ranges.ccp4")))
        coords = get_clusters_centre_mass(ensemble_probe_clusters, probe_array)

        for c in coords.keys():
            clust_ids.append(f"{probe}_{c}")
            clust_size.append(len(ensemble_probe_clusters[ensemble_probe_clusters == c]))
            clust_centroid.append(coords[c])
            clust_medians.append(np.median(probe_array[ensemble_probe_clusters == c]))
            map_size.append(probe_clust_grid.count_grid())
            clust_contributing.append([ensemble_res.index_dict[ci[0]] for ci in contribs[c]])
            clust_contributions.append([ci[1] for ci in contribs[c]])

    probe_df = pd.DataFrame()
    probe_df['clust_id'] = clust_ids
    probe_df['clust_size'] = clust_size
    probe_df['num_points_map'] = map_size
    probe_df['centroid'] = clust_centroid
    probe_df['clust_median'] = clust_medians
    probe_df['contributing_structures'] = clust_contributing
    probe_df['structure_contributions'] = clust_contributions
    return probe_df

def make_brd1_summary_ensemble(brd1_selmaps, ensemble_dir):
    hs_results = []
    for p in brd1_selmaps:
        with HotspotReader(str(p)) as hrd:
            res = hrd.read()
            for probe, g in res.super_grids.items():
                #if probe == 'apolar':
                g = g* (g > 10.0)
                    
                res.super_grids[probe] = g
            res.protein.identifier = str(p).split('over_')[1].split('_')[0]
            hs_results.append(res)
            
    ensemble_settings = EnsembleResult.Settings()
    ensemble_settings.combine_mode = 'median'
    ensemble_settings.polar_frequency_threshold = 0.0
    ensemble = EnsembleResult(hs_results_list=hs_results,
                              ensemble_id=f"BRD1_selectivity",
                              settings=ensemble_settings)
    ensemble.make_ensemble_maps(save_grid_ensembles=True)
    ensemble_hs_result = ensemble.ensemble_hotspot_result
    ensemble_df = ensemble_cluster_summary(ensemble, save_dir=ens_dir)
    
    with HotspotWriter(str(Path(ensemble_dir, "ensemble_selectivities"))) as hwr:
        hwr.write(ensemble_hs_result)
    return ensemble_df

def make_struct_df(selectivity_df, aligned_compounds_dir):
    struct_list = []
    cid_list = []
    on_targets = []
    off_targets = []
    crystal_targets = []
    act_types = []
    act_ratios = []
    
    for i, row in selectivity_df.iterrows():
        compound_save_dir = Path(aligned_compounds_dir, row['compound_chembl'])
        structs_to_score = list(compound_save_dir.glob("*_aligned.sdf"))
        
        for struct in structs_to_score:
            crystal_target = struct.name.split('_')[0]
            struct_name = struct.stem
            
            struct_list.append(struct_name)
            cid_list.append(row['compound_chembl'])
            on_targets.append(row['on_target'])
            off_targets.append(row['off_target'])
            crystal_targets.append(crystal_target)
            act_types.append(row['activity_type'])
            act_ratios.append(row['activity_ratio'])
            
    all_df = pd.DataFrame()
    all_df['structure_name'] = struct_list
    all_df['crystal_target'] = crystal_targets
    all_df['compound_chembl_id'] = cid_list
    all_df['on_target'] = on_targets
    all_df['off_target'] = off_targets
    all_df['activity_type'] = act_types
    all_df['activity_ratio'] = act_ratios
    
    return all_df
    
    
if __name__ == "__main__":
    np.random.seed(1500)
    
    sel_dir = Path('../results/select_maps_all').resolve()
    ens_dir = Path('../results/BRD1_selectivity_ensemble')
    compounds_dir = Path('../results/aligned_chembl_compounds')
    
    if not ens_dir.exists():
        ens_dir.mkdir()
    
    lig_paths = list(compounds_dir.glob("*/*_aligned.sdf"))
    selmap_paths = list(sel_dir.glob('*/out/hotspot'))
    #brd1_select_maps, brd1_select_maps_paths = shrink_bs_maps(selmap_paths, lig_paths, padding=4.0)
    brd1_select_maps_paths = [x for x in sel_dir.glob('*/out/binding_site_maps/out.zip') if 'BRD1_over' in str(x)]
    
    ens_df = make_brd1_summary_ensemble(brd1_select_maps_paths, ens_dir)
    ens_df.to_csv(str(Path(ens_dir, "BRD1_sel_cluster_summary.csv")))
    #ens_df = pd.read_csv(str(Path(ens_dir, "BRD1_sel_cluster_summary.csv")), index_col=0)
        
    
    acc_grid = Grid.from_file(str(Path(ens_dir, "acceptor_ensemble_clusters_ranges.ccp4")))
    #cropped_acc_grid = (acc_grid > 0) * (acc_grid < 2)
    cropped_acc_grid = acc_grid > 1
    apolar_grid = Grid.from_file(str(Path(ens_dir, "apolar_ensemble_clusters_ranges.ccp4")))
    cropped_apolar_grid1 = (apolar_grid > 0) * (apolar_grid < 2)

    
    sele_df = pd.read_csv(str(Path("../results/chembl_compounds/selectivities_chembl_compounds.csv")), index_col=0)
    compound_df = make_struct_df(selectivity_df=sele_df, aligned_compounds_dir=compounds_dir)
    compound_df.to_csv(Path(compounds_dir, "all_structures_selratios.csv"))
    
    sub_df = compound_df[compound_df['on_target'] == 'BRD1'].copy()
    
    cluster_scores_acceptor = []
    cluster_scores_apolar1 = []
    
    
    ens_compounds_dir = Path(ens_dir, 'ensemble_compounds')
    if not ens_compounds_dir.exists():
        ens_compounds_dir.mkdir()
    
    for i, row in sub_df.iterrows():
        comp_name = row['compound_chembl_id']
        sname = row['structure_name']
        mol = MoleculeReader(str(Path(compounds_dir, comp_name, f"{sname}.sdf")))[0]
        if len(mol.components) > 1:
            mol = mol.components [0]
        mol.add_hydrogens()
        rd_mol = Chem.MolFromMol2Block(mol.to_string())
        AllChem.Compute2DCoords(rd_mol)
        Draw.MolToFile(rd_mol,str(Path(ens_compounds_dir, f"{sname}.png"))) 
        with MoleculeWriter(str(Path(ens_compounds_dir, f"{sname}.sdf"))) as mwr:
            mwr.write(mol)
        
        sum_scores = sum([cropped_acc_grid.value_at_coordinate(a.coordinates, tolerance=1, position=False) for a in mol.heavy_atoms if a.is_acceptor])
        if sum_scores> 0:
            cluster_scores_acceptor.append(1)
        else:
            cluster_scores_acceptor.append(0)
            
        sum_scores_apolar1 = sum([cropped_apolar_grid1.value_at_coordinate(a.coordinates, tolerance=1, position=False) for a in mol.heavy_atoms if (not (a.is_acceptor or a.is_donor)) and (not (a.is_donor and a.is_acceptor))]) 
        
        if sum_scores_apolar1> 0:
            cluster_scores_apolar1.append(1)
        else:
            cluster_scores_apolar1.append(0)
                
    sub_df['BRD1_acceptor_hit'] = cluster_scores_acceptor
    sub_df['BRD1_apolar_hit'] = cluster_scores_apolar1
    
    sub_df.to_csv(str(Path(ens_dir, 'BRD1_cluster_hits.csv')))
    
    all_df = pd.read_csv('../results/chembl_compounds/chembl_all_targets_curated.csv', index_col=0)
    
    chembl_comp_ids = sub_df['compound_chembl_id'].unique()
    PDB_comp_ids = []
    PDB_prot_structs = []
    smiles = []
    for cid in chembl_comp_ids:
        structs = list(Path(compounds_dir, cid).glob("*_aligned.sdf"))
        pdb_structs = list(set([s.name.split('_')[2] for s in structs]))
        pdb_comps = list(set([s.name.split('_')[3] for s in structs]))
        smile = all_df[all_df['molecule_chembl_id'] == cid]['canonical_smiles'].values
        
        PDB_comp_ids.append(pdb_comps)
        PDB_prot_structs.append(pdb_structs)
        smiles.append(list(set(smile)))
        
    sup_df = pd.DataFrame()
    sup_df['PDB identifier'] = PDB_comp_ids
    sup_df['Chembl identifier'] = chembl_comp_ids
    sup_df['PDB structures'] = PDB_prot_structs
    sup_df['Smiles'] = smiles
    
    sup_df.to_csv(Path(ens_dir, "supplementary_table_7.csv"))