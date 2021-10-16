from pathlib import Path
from ccdc.protein import Protein
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.hs_ensembles import EnsembleResult, SelectivityResult
from siena_ensemble import SienaQuery
from utils import process_siena_pdbs, get_subset, shrink_bs_maps
from run_hotspots_job import run_parallel_hotspot_jobs
import random
import json

# Set up the targets, query pdbs and query ligands
ensemble_dict = {'on_target':{'target_name': 'BRD7',
                              'pdb_code': '4UIU',
                              'ligand_name':'TVU_A_1123'},
                 'off_target': {'target_name': 'BRD9',
                                'pdb_code': '5MQ1',
                                'ligand_name': '5U6_A_700'}}

# Get the ensembles from SIENA. Remove redundant ligands (keep highest res) and anything > 300 Da
for ens in ensemble_dict.keys():
    ens_dir = Path(f"../ensembles/{ensemble_dict[ens]['target_name']}")
    if not ens_dir.exists(): ens_dir.mkdir()
    # Check that we haven't already queried SIENA:
    if not Path(ens_dir, 'subset_ensemble_data.csv').exists():
        print('Ensemble summary data .csv not found. Querying SIENA')
        ens_query = SienaQuery(target_name=ensemble_dict[ens]['target_name'],
                               pdb_code=ensemble_dict[ens]['pdb_code'],
                               ligand_name=ensemble_dict[ens]['ligand_name'], 
                               siena_dir=ens_dir, 
                              )

        ens_query.get_ensemble_protein_df()
    ens_sub_ens =get_subset(Path(ens_dir, 'ensemble_data.csv'), data_source='SIENA', mw_threshold=1000.0)
    ensemble_dict[ens]['ensemble_df_path'] = Path(ens_dir, 'subset_ensemble_data.csv')
    ens_sub_ens.to_csv(str(ensemble_dict[ens]['ensemble_df_path']))

# Align the two ensembles prior to recalculating the hotspots
# Works in this case because the reference structure only has one chain - fix later.
ref_prot_path = Path("../../case_studies/bromodomains_protoss_csd2021/BRD1/pdb_files/5AMF_7.pdb")
ref_prot = Protein.from_file(str(ref_prot_path))

on_target_prot_path = list(Path(ensemble_dict['on_target']['ensemble_df_path'].parent, 'pdb_files').glob(f'{ensemble_dict["on_target"]["pdb_code"]}_*.pdb'))[0]

on_target_ref_prot = Protein.from_file(str(on_target_prot_path))

off_target_prot_path = list(Path(ensemble_dict['off_target']['ensemble_df_path'].parent, 'pdb_files').glob(f'{ensemble_dict["off_target"]["pdb_code"]}_*.pdb'))[0]

off_target_ref_prot = Protein.from_file(str(off_target_prot_path))

# Superpose the off-target reference protein based on the binding site residues:
sup = Protein.ChainSuperposition()
sets = sup.Settings()
sets.superposition_atoms = "BACKBONE"
tar_binding_site = Protein.BindingSiteFromMolecule(ref_prot,ref_prot.ligands[0], 6.0)
(rmsd, on_transform) = sup.superpose(ref_prot.chains[0], on_target_ref_prot.chains[0], tar_binding_site)
(rmsd, off_transform) = sup.superpose(ref_prot.chains[0], off_target_ref_prot.chains[0], tar_binding_site)

on_tar_paths_chains = process_siena_pdbs(res_table_path=ensemble_dict['on_target']['ensemble_df_path'], transformation=on_transform)
print(on_tar_paths_chains)
ensemble_dict['on_target']['input_protein_paths'] = [x for x in on_tar_paths_chains]
print(ensemble_dict['on_target']['input_protein_paths'])

off_tar_paths_chains = process_siena_pdbs(res_table_path=ensemble_dict['off_target']['ensemble_df_path'],
                                                                           transformation=off_transform)

ensemble_dict['off_target']['input_protein_paths'] = [n for n in off_tar_paths_chains]
print(off_tar_paths_chains)


# Calculate the hotspot maps:
for ens in ensemble_dict.keys():
    random.seed(123)
    ensemble_dict[ens]['fullsize_hotspots'] = run_parallel_hotspot_jobs([str(p) for p in ensemble_dict[ens]['input_protein_paths']],
                                                   hs_charged=False,
                                                   hs_rotations=3000,
                                                   hs_protonated=False,
                                                   hs_spheres=False)

# # SHrink the hotspot maps to the binding site
# lig_paths = [list(Path(p).parent.glob('*_ligands.sdf'))[0] for ens in ensemble_dict.keys() for p in ensemble_dict[ens]['input_protein_paths']]
lig_paths = list(Path("/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/revisions/brd7_brd9/scored_chembl_compounds").glob("*/*_aligned.sdf"))

# Calculate the ensemble maps
for ens in ensemble_dict.keys():
   # shrunk_hotspots, shrunk_hot_paths = shrink_bs_maps(hotspot_paths=ensemble_dict[ens]['fullsize_hotspots'],
   #                                                    ligand_paths=lig_paths)
#    shrunk_hotspots = [HotspotReader(str(Path(f, 'binding_site_maps', 'out'))).read() for f in Path(ensemble_dict[ens]['ensemble_df_path']).parent.glob('hotspot_results/*/fullsize_hotspots_3000')]
   shrunk_hotspots, shrunk_hot_paths = shrink_bs_maps(hotspot_paths= ensemble_dict[ens]['fullsize_hotspots'],
                                                      ligand_paths=lig_paths,
                                                      padding=5)

   ensemble_settings = EnsembleResult.Settings()
   ensemble_settings.combine_mode = 'median'
   ensemble_settings.apolar_frequency_threshold = None
   ensemble_settings.polar_frequency_threshold = 20.0

   ensemble = EnsembleResult(hs_results_list=shrunk_hotspots,
                             ensemble_id=ensemble_dict[ens]['target_name'],
                             settings=ensemble_settings)

   ensemble.make_ensemble_maps(save_grid_ensembles=False)
   ensemble_hs_result = ensemble.ensemble_hotspot_result
   ens_path = Path(Path(ensemble_dict[ens]['ensemble_df_path']).parent, 'ensemble_maps')

   with HotspotWriter(str(ens_path), grid_extension=".ccp4", zip_results=True) as w:
       w.write(ensemble_hs_result)
   ensemble_dict[ens]['ensemble_maps'] = Path(ens_path, 'out.zip')

# Calculate the selectivity maps
selectivity_settings = SelectivityResult.Settings()
selectivity_settings.cluster_distance_cutoff = 3.0
selectivity_settings.minimal_cluster_score = 10.0

with HotspotReader(ensemble_dict['on_target']['ensemble_maps']) as onhr:
    tar_res = onhr.read()

with HotspotReader(ensemble_dict['off_target']['ensemble_maps']) as offhr:
    off_res = offhr.read()

on_over_off = SelectivityResult(target_result=tar_res,
                                       other_result=off_res
                                )
on_over_off.make_selectivity_maps()

bromo_dir = ens_dir.parent

on_target = ensemble_dict['on_target']['target_name']
off_target = ensemble_dict['off_target']['target_name']
with HotspotWriter(str(Path(bromo_dir, f"selectivity_{on_target}_over_{off_target}")),grid_extension=".ccp4", zip_results=False) as w:
    w.write(on_over_off.selectivity_result)

out_file = Path(ens_dir.parent, 'ensemble_data.txt')

# with open(str(out_file), 'w') as out:
#     json.dump(ensemble_dict, out, indent=4)