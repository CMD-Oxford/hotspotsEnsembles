import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import random
from ccdc.io import MoleculeReader, MoleculeWriter
from ccdc.protein import Protein
from hotspots.hs_io import HotspotWriter, HotspotReader
from hotspots.hs_ensembles import EnsembleResult, SelectivityResult
from utils import get_subset, process_KLIFS_pdbs, find_bs_residues, shrink_bs_maps, KLIFS_molecular_weight
from run_hotspots_job import run_parallel_hotspot_jobs


def run_example(ensemble_dict, rseed=123):

    for case_study in ensemble_dict:

        on_target_dir = ensemble_dict[case_study]['on_target_dir']
        off_target_dir = ensemble_dict[case_study]['off_target_dir']
        on_target = on_target_dir.name
        off_target = off_target_dir.name

        # For the reference ensemble (p38a), get a 'super ligand' from all the bound ligands and use that to define the binding site
        reference_struct = Protein.from_file(str(ensemble_dict[case_study]['ref_structure_path']))
        ligands = [MoleculeReader(str(m))[0] for m in list(on_target_dir.glob('*/ligand.mol2'))]
        ligands.extend([MoleculeReader(str(m))[0] for m in list(on_target_dir.glob('*/ligand.mol2'))])

        bs_residues = find_bs_residues(prot=reference_struct, list_of_ligands=ligands)
        reference_binding_site = Protein.BindingSiteFromListOfResidues(protein=reference_struct, list_of_residues=bs_residues)
        chain_superposition = Protein.ChainSuperposition()
        sets = chain_superposition.Settings()
        sets.superposition_atoms = "CALPHA"

        path_dict = {}

        for ens_dir in [on_target_dir, off_target_dir]:
            path_dict[ens_dir.name] = {}
            in_prots = list(Path(ens_dir).glob('*/protein.mol2'))
            for p in in_prots:
                c_prot = Protein.from_file(str(p))
                (rmsd, transformation) = chain_superposition.superpose(reference_struct.chains[0], c_prot.chains[0], reference_binding_site)
                report_string = f"Structure {p.parent.name}, RMSD: {rmsd} \n"

                with MoleculeWriter(str(p).replace('.mol2', '_aligned.mol2')) as wr:
                    wr.write(c_prot)

                lig = MoleculeReader(str(p).replace('protein.mol2', 'ligand.mol2'))[0]
                lig.transform(transformation)

                with MoleculeWriter(str(p).replace('protein.mol2', 'ligand_aligned.mol2')) as wl:
                    wl.write(lig)

                complex = Protein.from_file(str(p).replace('protein.mol2', 'complex.mol2'))
                complex.transform(transformation)
                with MoleculeWriter(str(p).replace('protein.mol2', 'complex_aligned.mol2')) as wc:
                    wc.write(complex)

                new_rmsd, new_trans =  chain_superposition.superpose(reference_struct.chains[0], c_prot.chains[0], reference_binding_site)
                report_string += f"Aligned structure {p.parent.name}, RMSD: {new_rmsd}"
                repor = Path(str(p).replace('protein.mol2', 'alignement_log.txt'))
                repor.write_text(report_string)

            KLIFS_molecular_weight(Path(ens_dir, 'overview.csv')) # get the mw of ligands in the ensemble dataframe
            ensemble_df = Path(ens_dir, 'overview.csv')
            subset_df = get_subset(ensemble_df, data_source='KLIFS', unique_ligands=True)
            subset_path = Path(ensemble_df.parent, 'subset_overview.csv')
            subset_df.to_csv(subset_path)

            ens_paths = process_KLIFS_pdbs(subset_path)
            #hot_paths = list(Path(ens_dir, 'hotspot_results').glob('*/fullsize_hotspots_3000/out.zip'))
            random.seed(rseed)
            hot_paths = run_parallel_hotspot_jobs(str_paths=[str(p) for p in ens_paths],
                                                  hs_charged=False,
                                                  hs_rotations=3000,
                                                  hs_protonated=False,
                                                  hs_spheres=False)

            path_dict[ens_dir.name]['protein_paths'] = ens_paths
            path_dict[ens_dir.name]['fullsize_hotspots'] = hot_paths

        # Calculate the ensemble maps
        for ens_dir in [on_target_dir, off_target_dir]:
            hot_paths = list(Path(ens_dir, 'hotspot_results').glob('*/fullsize_hotspots_3000/out.zip'))
            lig_paths = [str(Path(hp.parents[1], 'ligand.mol2')) for hp in hot_paths]
            print(lig_paths)
            shrunk_hotspots, shrunk_hot_paths = shrink_bs_maps(hotspot_paths=path_dict[ens_dir.name]['fullsize_hotspots'],
                                                              ligand_paths=lig_paths)

            ensemble_settings = EnsembleResult.Settings()
            ensemble_settings.combine_mode = 'median'
            ensemble_settings.apolar_frequency_threshold = None
            ensemble_settings.polar_frequency_threshold = 20.0
            ensemble = EnsembleResult(hs_results_list=shrunk_hotspots,
                                     ensemble_id=ens_dir.name,
                                     settings=ensemble_settings)
            ensemble.make_ensemble_maps(save_grid_ensembles=False)
            ensemble_hs_result = ensemble.ensemble_hotspot_result
            ens_path = Path(ens_dir, 'ensemble_maps')

            with HotspotWriter(str(ens_path), grid_extension=".ccp4", zip_results=True) as w:
               w.write(ensemble_hs_result)
            path_dict[ens_dir.name]['ensemble_maps'] = Path(ens_path, 'out.zip')

        # Calculate the selectivity maps
        selectivity_settings = SelectivityResult.Settings()
        selectivity_settings.apolar_percentile_threshold = 90.0
        selectivity_settings.cluster_distance_cutoff = 3.0
        selectivity_settings.minimal_cluster_score = 10.0

        with HotspotReader(path_dict[on_target]['ensemble_maps']) as tar_reader:
            tar_res = tar_reader.read()

        with HotspotReader(path_dict[off_target]['ensemble_maps']) as off_reader:
            off_res = off_reader.read()

        on_over_off = SelectivityResult(target_result= tar_res,
                                        other_result=off_res
                                        )
        on_over_off.settings = selectivity_settings
        on_over_off.make_selectivity_maps()

        case_study_path = on_target_dir.parent
        
        on_target = on_target_dir.name.split('_')[1]
        off_target = off_target_dir.name.split('_')[1]

        with HotspotWriter(str(Path(case_study_path, f"selectivity_{on_target}_over_{off_target}")),grid_extension=".ccp4", zip_results=False) as w:
            w.write(on_over_off.selectivity_result)

if __name__ == '__main__':
    # e_dict = {'case_study_1': {'on_target_dir': Path('../case_studies/KLIFS_kinases_2021/KLIFS_p38a_fragments_09_02_2021'),
    #                                   'off_target_dir': Path('../case_studies/KLIFS_kinases_2021/KLIFS_ERK2_fragments_09_02_2021'),
    #                                   'ref_structure_path': Path(
    #                                       '../case_studies/KLIFS_kinases_2021/KLIFS_p38a_fragments_09_02_2021/1w7h_chainA/protein.mol2')
    #                                   },
    #                  'case_study_2': {'on_target_dir': Path('../case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020'),
    #                                   'off_target_dir': Path('../case_studies/KLIFS_kinases_2021/KLIFS_PIM1_01_10_2020'),
    #                                   'ref_structure_path': Path(
    #                                       '../case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020/1pjk_altA_chainA/protein.mol2')
    #                                   }
    #                  }

    e_dict = {'case_study_2':
                  {'on_target_dir': Path('../case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020'),
                   'off_target_dir': Path('../case_studies/KLIFS_kinases_2021/KLIFS_PIM1_01_10_2020'),
                   'ref_structure_path': Path('../case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020/1pjk_altA_chainA/protein.mol2')
                   }
              }

    run_example(e_dict, rseed=123)



