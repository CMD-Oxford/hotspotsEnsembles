from pathlib import Path
import pandas as pd
import numpy as np
import utils
from ccdc.protein import Protein
from ccdc.io import MoleculeReader
from hotspots.hs_utilities import Helper


info_dict = {'BRD1': Path('../case_studies/bromodomains_protoss_csd2021/BRD1/subset_ensemble_data.csv'),
             'BRPF1': Path('../case_studies/bromodomains_protoss_csd2021/BRPF1/subset_ensemble_data.csv'),
             'ERK2': Path(
                 '../case_studies/KLIFS_kinases_2021/KLIFS_ERK2_fragments_09_02_2021/subset_overview.csv'),
             'p38a': Path(
                 '../case_studies/KLIFS_kinases_2021/KLIFS_p38a_fragments_09_02_2021/subset_overview.csv'),
             'CK2a': Path('../case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020/subset_overview.csv'),
             'PIM1': Path('../case_studies/KLIFS_kinases_2021/KLIFS_PIM1_01_10_2020/subset_overview.csv')
             }
bs_dict = {}

for ens, ens_df_path in info_dict.items():
    print(ens)
    ens_df = pd.read_csv(ens_df_path)
    if (ens == 'BRD1') or (ens == 'BRPF1'):
        ens_calpha_coords = []
        for idx, row in ens_df.iterrows():
            bs_residues = row['Binding site']
            bs_residues = bs_residues.replace("[", '').replace("]", '').split(',')
            bs_residues = [x.replace("'", '').strip() for x in bs_residues]
            pname = Path(row['Filename']).stem
            prot_path = Path(ens_df_path.parent, 'hotspot_results', pname, 'fullsize_hotspots_3000', 'hotspots_input_protein.pdb')
            prot = Protein.from_file(str(prot_path))
            c_alpha_coords = []

            for bs in bs_residues:
                res = [r for r in prot.residues if r.identifier==bs][0]
                c_alpha_coords.append(np.array(res.c_alpha.coordinates))
            ens_calpha_coords.append(np.array(c_alpha_coords))

    else:
        ref_prot = None
        ref_bs_residues = []
        ligand_paths = list(Path(ens_df_path.parent, 'hotspot_results').glob('*/ligand.mol2'))
        ligands = [MoleculeReader(str(l))[0] for l in ligand_paths]

        ens_calpha_coords = []
        for idx, row in ens_df.iterrows():
            if row['alt'] == ' ':
                pname = "{}_chain{}".format(row['pdb'].strip(), row['chain'].strip())
            else:
                pname = "{}_alt{}_chain{}".format(row['pdb'].strip(), row['alt'].strip(), row['chain'].strip())

            prot = Protein.from_file(str(Path(ens_df_path.parent, 'hotspot_results', pname, 'fullsize_hotspots_3000', 'hotspots_input_protein.pdb')))
            print(pname)

            if ref_prot is None:
                ref_prot = prot
                ref_bs_residues = utils.find_bs_residues(ref_prot, ligands)
                ens_calpha_coords.append(np.array([np.array(r.c_alpha.coordinates) for r in ref_bs_residues]))
                continue

            nearest_res = []
            for res in ref_bs_residues:
                tar_resis = [pr for pr in prot.residues if ((res.three_letter_code == pr.three_letter_code) and (pr.c_alpha.coordinates is not None))]
                print([(tr.identifier, tr.c_alpha.coordinates) for tr in tar_resis])
                dists =np.array([Helper.get_distance(res.c_alpha.coordinates, tr.c_alpha.coordinates) for tr in tar_resis])
                nearest_idx = np.argmin(dists)
                print(res.identifier, tar_resis[nearest_idx].identifier, dists[nearest_idx])
                if res.three_letter_code == tar_resis[nearest_idx].three_letter_code:
                    nearest_res.append(tar_resis[nearest_idx])
                else:
                    print('PROBLEM MATCHING RESIDUES')
                    break
            ens_calpha_coords.append(np.array([np.array(nr.c_alpha.coordinates) for nr in nearest_res]))

    # Get the mean coordinates of the binding site residues
    coords_ar = np.array(ens_calpha_coords)
    mean_ar = np.mean(coords_ar, axis=0)
    rmsds = [np.sqrt(np.sum((mean_ar - ens_calpha_coords[i]) ** 2) / mean_ar.shape[0]) for i in
             range(len(ens_calpha_coords))]
    bs_dict[ens] = np.array(rmsds)

df =  pd.DataFrame(columns=['Ensemble', 'RMSD', 'StDev', 'num_structures'])
ensembles = list(bs_dict.keys())
mean_rmsd = [round(np.mean(bs_dict[d]),2) for d in ensembles]
stdev_rmsd = [round(np.std(bs_dict[d]),2) for d in ensembles]
num_struct = [len(bs_dict[d]) for d in ensembles]

df['Ensemble'] = ensembles
df['RMSD'] = mean_rmsd
df['StDev'] = stdev_rmsd
df['num_structures'] = num_struct

df.to_csv('ensemble_flexibility_data.csv')