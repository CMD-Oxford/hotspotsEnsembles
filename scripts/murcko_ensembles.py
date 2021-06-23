from pathlib import Path
import numpy as np
import pandas as pd
from hdbscan import HDBSCAN
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity
from rdkit.Chem.Scaffolds import MurckoScaffold
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.hs_ensembles import EnsembleResult
from parameter_run_analysis import visualise_ensmap_cluster_contribs, reduce_clusters
from calculate_ensemble_maps import ensemble_cluster_summary

def get_similar_molecules(mols, cutoff=0.9):
    """
    given a list of molecules as rdkit_mols or smiles, groups anything with > cutoff Tanimoto similarity and returns dict of indices
    :param mols: list of smiles strings
    :return:
    """
    if all(type(m) is str for m in mols):
        print("Received molecules as strings, assuming smiles")
        calc_mols = [Chem.MolFromSmiles(s) for s in mols]

    elif all(type(m) is Chem.rdchem.Mol for m in mols):
        print("Received molecules as RDKit molecules, continuing")
        calc_mols = mols

    else:
        print(
            "Unrecognised molecule type. For converting ccdc molecules, use rd_mol = rdkit.Chem.MolFromMol2Block(ccdc_mol.to_string())")
        return

    # Caclulate the similarity matrix:
    fps = [AllChem.GetMorganFingerprintAsBitVect(m, 3, 1024) for m in calc_mols]

    distances = []
    nfps = len(fps)
    for i in range(nfps):
        sims = BulkTanimotoSimilarity(fps[i], fps)
        distances.append(sims)

    dists = np.array(distances)
    # The matrix is symmetric, and 1 along the diagonal
    dists = np.triu(dists, 1)

    print(f"Length of mols list: {len(mols)}, number of fingerprints: {nfps} ")
    passed_indxs = []
    similarity_dict = {}
    for mol_ind in range(nfps): #
        if mol_ind not in passed_indxs:
            dist_slice = dists[mol_ind]
            similar_compound_idxs = [x for x in np.where(dist_slice>=cutoff)[0]]
            similar_compound_idxs.append(mol_ind)
            passed_indxs.extend(similar_compound_idxs)
            similarity_dict[mol_ind] = similar_compound_idxs # This is all returning the indexes, not actual mols

    return similarity_dict

def plot_murcko_ligands(clust_dict, murcko_image_save_dir, frag_sub_df):
    for clust_ind, comp_idxs in clust_dict.items():
        # murcko_scaff = Chem.MolFromSmiles(murcko_smiles[clust_ind])
        murcko_scaff = murcko_smiles_generic[clust_ind]
        AllChem.Compute2DCoords(murcko_scaff)
        Draw.MolToFile(murcko_scaff, str(Path(murcko_image_save_dir, f"cluster_{clust_ind}_murcko.png")))
        cluster_ligands = [Chem.MolFromSmiles(smiles[a]) for a in comp_idxs]
        cluster_ligands_scaffs = [murcko_smiles_generic[a] for a in comp_idxs]
        master_ligs = []
        master_labels = []
        for it, ci in enumerate(comp_idxs):
            master_ligs.extend([cluster_ligands[it], cluster_ligands_scaffs[it]])
            master_labels.extend([frag_sub_df['PDB code'].values[ci]]*2)

        for cl in master_ligs:
            AllChem.Compute2DCoords(cl)
        print(master_labels)
        # labels = [frag_sub_df['PDB code'].values[b] for b in comp_idxs]
        img = Chem.Draw.MolsToGridImage(master_ligs, molsPerRow=2, subImgSize=(200, 200),
                                        legends=master_labels)
        img.save(str(Path(murcko_image_save_dir, f"cluster_{clust_ind}_ligands.png")))

def hotspot_list_from_subset(subset_df, ens_df_path):
    hotspot_list = []

    for idx, row in subset_df.iterrows():
        # Find the protein file:
        if (ens_name == 'BRD1') or (ens_name == 'BRPF1'):  # bromodomains have slightly different overviews
            pname = Path(row['Filename']).stem
            print(pname)
        else:  # assume KLIFS summary

            if row['alt'] == ' ':
                pname = "{}_chain{}".format(row['pdb'].strip(), row['chain'].strip())
            else:
                pname = "{}_alt{}_chain{}".format(row['pdb'].strip(), row['alt'].strip(), row['chain'].strip())

        hs_path = Path(ens_df_path.parent, 'hotspot_results', pname, 'fullsize_hotspots_3000', 'binding_site_maps',
                       'out.zip')

        if hs_path.exists():
            with HotspotReader(str(hs_path)) as hsr:
                res = hsr.read()
                res.protein.identifier = pname
                hotspot_list.append(res)
            # print('ok')
        else:
            print(f"No hotspot maps found for {pname}. Re-run the case study script")
    return hotspot_list

def make_subset_ensemble(clust_subdf, ens_df_path, clust_save_dir):
    cluster_hotspots = hotspot_list_from_subset(clust_subdf, ens_df_path)
    ensemble_settings = EnsembleResult.Settings()
    ensemble_settings.combine_mode = 'median'
    ensemble_settings.apolar_frequency_threshold = None
    ensemble_settings.polar_frequency_threshold = 20.0
    ensemble = EnsembleResult(hs_results_list=cluster_hotspots,
                              ensemble_id=f"{ens_name}_murcko",
                              settings=ensemble_settings)
    ensemble.make_ensemble_maps(save_grid_ensembles=True)
    ensemble_hs_result = ensemble.ensemble_hotspot_result
    ens_path = Path(clust_save_dir, f'ensemble_maps')
    if not ens_path.exists():
        ens_path.mkdir()
    try:
        ens_df = ensemble_cluster_summary(ensemble, save_dir=ens_path)
        ens_df.to_csv(str(Path(ens_path, "ensemble_cluster_summary.csv")))

    except ValueError as ve:
        print(ve)
        return

    with HotspotWriter(str(ens_path), grid_extension=".ccp4", zip_results=True) as w:
        w.write(ensemble_hs_result)
    return Path(ens_path, "ensemble_cluster_summary.csv")

def make_ens_map_cluster_clusters(df_paths, ens_name, ens_dir):
    df_list = []

    for dcd in df_paths:
        df = pd.read_csv(Path(dcd, 'ensemble_cluster_summary.csv'), index_col=0)

        df['probe_type'] = [i.split('_')[0] for i in df['clust_id'].values]
        df['on_target'] = ens_name

        df_list.append(df)
    master_df = pd.concat(df_list)

    targets = list(set(master_df['on_target']))
    probes = list(set(master_df['probe_type']))
    table_path = Path(ens_dir, 'parameter_tables')
    if not table_path.exists():
        table_path.mkdir()
    save_paths = []
    for t in targets:
        for p in probes:
            print(t, p)
            sub_df = master_df[master_df['on_target'] == t]
            sub_df = sub_df[sub_df['probe_type'] == p]
            center_list = [x.split('[')[1].split(']')[0].split(',') for x in sub_df['centroid'].values]
            center_list = np.array([list(map(float, x)) for x in center_list])
            if len(center_list) > 0:
                clusterer = HDBSCAN(min_cluster_size=2)
                clusterer.fit(center_list)
                labels = reduce_clusters(center_list=center_list, labels=clusterer.labels_, threshold=0.5)
                sub_df['cluster_cluster'] = labels
                sub_df = sub_df.sort_values(by='cluster_cluster')
                s_path = Path(table_path, f"{t}_{p}.csv")
                sub_df.to_csv(s_path)
                save_paths.append(s_path)
    return save_paths

if __name__ == "__main__":
    info_dict = {'BRD1': Path('../case_studies/bromodomains_protoss_csd2021/BRD1/subset_ensemble_data.csv'),
                 'BRPF1': Path('../case_studies/bromodomains_protoss_csd2021/BRPF1/subset_ensemble_data.csv'),
                 'ERK2': Path(
                     '../case_studies/KLIFS_kinases_2021/KLIFS_ERK2_fragments_09_02_2021/subset_overview.csv'),
                 'p38a': Path(
                     '../case_studies/KLIFS_kinases_2021/KLIFS_p38a_fragments_09_02_2021/subset_overview.csv'),
                 'CK2a': Path('../case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020/subset_overview.csv'),
                 'PIM1': Path('../case_studies/KLIFS_kinases_2021/KLIFS_PIM1_01_10_2020/subset_overview.csv')
                 }


    fragment_summary_path = Path('../case_studies/ensemble_summary_fragments_new.csv')
    frag_df = pd.read_csv(fragment_summary_path, index_col=0)

    for ens_name in info_dict.keys():
        print(ens_name)
        ens_path = info_dict[ens_name].parent
        murcko_save_dir = Path(ens_path, 'Murcko_ensembles')
        if not murcko_save_dir.exists():
            murcko_save_dir.mkdir()
        murcko_image_save_dir = Path(murcko_save_dir, 'cluster_images')
        if not murcko_image_save_dir.exists():
            murcko_image_save_dir.mkdir()
        # Now get the ensemble fragment smiles
        frag_sub_df = frag_df[frag_df['Target'] == ens_name]
        smiles = frag_sub_df['Ligand SMILES'].values
        murcko_smiles = [MurckoScaffold.MurckoScaffoldSmiles(sm) for sm in smiles]
        murcko_smiles_generic = [MurckoScaffold.MakeScaffoldGeneric(Chem.MolFromSmiles(s)) for s in murcko_smiles]
        clust_dict = get_similar_molecules(murcko_smiles_generic, cutoff=1.0)
        plot_murcko_ligands(clust_dict, murcko_image_save_dir, frag_sub_df)
        frag_sub_df['Murcko_smiles'] = [Chem.MolToSmiles(mm) for mm in murcko_smiles_generic]
        murcko_reverse_clusters = {val:key for key, x in clust_dict.items() for val in x}
        frag_sub_df['Murcko_cluster'] = [murcko_reverse_clusters[i] for i in range(len(murcko_smiles_generic))]

        ens_df = pd.read_csv(info_dict[ens_name], index_col=0)
        murcko_clusters = []
        for i, row in ens_df.iterrows():
            if (ens_name == 'BRD1') or (ens_name == 'BRPF1'):
                pdb_id = row['ID'].split('-')[0]
                chain = row['ID'].split('-')[1]

            else:
                pdb_id = row['pdb'].upper()
                chain = row['chain']
            tar_row = frag_sub_df[(frag_sub_df['PDB code'] == pdb_id) & (frag_sub_df['Chain'] == chain)]
            if len(tar_row) > 1:
                print(tar_row)
                #tar_row = tar_row.iloc[0]
            murcko_clusters.append(tar_row['Murcko_cluster'].values[0])

        ens_df['Murcko_cluster']= murcko_clusters
        frag_sub_df.to_csv(Path(murcko_save_dir, f"{ens_name}_fragment_subset_murckos.csv"))
        ens_df.to_csv(Path(murcko_save_dir, f"ensemble_subset_hotspot_input.csv"))

        clust_ensemble_summaries = []
        for c in clust_dict.keys(): # Make ensembles based on the Murcko clusters
            if len(ens_df[ens_df['Murcko_cluster'] == c])>= 4:
                clust_subdf = ens_df[ens_df['Murcko_cluster'] == c]
                clust_save_dir = Path(murcko_save_dir, f"murcko_cluster_{c}_ensemble")
                if not clust_save_dir.exists():
                    clust_save_dir.mkdir()
                ens_summary_path = make_subset_ensemble(clust_subdf, info_dict[ens_name], clust_save_dir)
                if ens_summary_path.exists():
                    clust_ensemble_summaries.append(ens_summary_path)
        visualise_ensmap_cluster_contribs(clust_ensemble_summaries, on_target=ens_name)

        #Now do the random Murckos
        # np.random.seed = 3
        # ens_df_clust = {c: list(np.where(ens_df['Murcko_cluster'] == c)[0]) for c in clust_dict.keys()}
        # print(ens_df_clust)
        # rand_clust_ensemble_summaries = []
        # for samp in range(5):
        #     sampled_idxs = []
        #     for k, v in ens_df_clust.items():
        #         sampled_idxs.append(v[np.random.randint(0, len(v))])
        #     print(sampled_idxs)
        #     rand_clust_df = ens_df.iloc[sampled_idxs, :]
        #     rand_clust_save_dir = Path(murcko_save_dir, f"random_murcko_cluster_{samp}")
        #     if not rand_clust_save_dir.exists():
        #         rand_clust_save_dir.mkdir()
        #     rand_clust_df.to_csv(Path(rand_clust_save_dir, f"murcko_random_cluster_{samp}_ensemble"))
        #     rand_clust_summary_path = make_subset_ensemble(rand_clust_df, info_dict[ens_name], rand_clust_save_dir)
        #     rand_clust_ensemble_summaries.append(rand_clust_summary_path)
        # visualise_ensmap_cluster_contribs(rand_clust_ensemble_summaries, on_target=ens_name)