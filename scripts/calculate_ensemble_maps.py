from pathlib import Path
import pandas as pd
from utils import get_subset, get_clusters_centre_mass
from hotspots.hs_io import HotspotReader, HotspotWriter
from hotspots.hs_ensembles import EnsembleResult
from hotspots.grid_extension import _GridEnsemble
import numpy as np
import matplotlib.pyplot as plt

def get_grid_dic_histograms(grid_dic, out_dir, prot_name, suffix=None):
    """
    Same as histograms function in main hotspots code, but for any grid_dic
    :param grid_dic: Python dictionary {[probe]:ccdc.Utilities.Grid}
    :param out_dir: str(path to out_dir)
    :param prot_name: str
    :param suffix: str
    :return:
    """
    data = {}
    for g in grid_dic.keys():
        grd = grid_dic[g]
        nx, ny, nz = grd.nsteps
        data[g] = np.array([grd.value(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz) if
                            round(grd.value(i, j, k)) != 0])
        # data[g] = grd.to_vector()

    plt.figure(1)
    for n, key in enumerate(data.keys()):

        plt.subplot(3, 1, (n + 1))
        # hotspot_result._histogram_info(data, key, n)
        colour_dict = {"acceptor": "r", "donor": "b", "apolar": "y"}
        hist, bin_edges = np.histogram(data[key], bins=range(0, 40), normed=True)
        plt.bar(bin_edges[:-1], hist, width=1, color=colour_dict[key])
        plt.xlim(min(bin_edges), max(bin_edges))
        plt.ylim(0, 0.35)
        plt.yticks([])
        if n == 0:
            plt.title("Fragment hotspot Maps")
        if n < 2:
            plt.xticks([])
        if n == 2:
            plt.xlabel("Fragment hotspot score")
        if n == 1:
            plt.ylabel("Frequency")
        plt.annotate(f"Points in map: {len(data[key])}", (0.70, 0.8), xycoords='axes fraction')
    if suffix != None:
        plt.savefig(Path(out_dir, (prot_name + suffix)))
    else:
        plt.savefig(Path(out_dir, prot_name))
    plt.close()

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


    for ens_name, ens_df_path in info_dict.items():

        subset_df = pd.read_csv(ens_df_path)

        hotspot_list = []

        for idx, row in subset_df.iterrows():
            # Find the protein file:
            if (ens_name == 'BRD1') or (ens_name == 'BRPF1'): # bromodomains have slightly different overviews
                pname = Path(row['Filename']).stem
                print(pname)
            else: # assume KLIFS summary

                if row['alt'] == ' ':
                    pname = "{}_chain{}".format(row['pdb'].strip(), row['chain'].strip())
                else:
                    pname = "{}_alt{}_chain{}".format(row['pdb'].strip(), row['alt'].strip(), row['chain'].strip())

            hs_path = Path(ens_df_path.parent, 'hotspot_results', pname, 'fullsize_hotspots_3000', 'binding_site_maps', 'out.zip')

            if hs_path.exists():
                with HotspotReader(str(hs_path)) as hsr:
                    res = hsr.read()
                    res.protein.identifier = pname
                    hotspot_list.append(res)
                #print('ok')
            else:
                print(f"No hotspot maps found for {pname}. Re-run the case study script")

        # Now make the ensembles

        thresholds = [None, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        # thresholds = [None, 0]

        image_dir = Path(ens_df_path.parent, 'ensemble_map_histograms')
        if not image_dir.exists():
            image_dir.mkdir()
        # break

        for t in thresholds:
            ensemble_settings = EnsembleResult.Settings()
            ensemble_settings.combine_mode = 'median'
            if t is not None:
                ensemble_settings.apolar_frequency_threshold = float(t)
                ensemble_settings.polar_frequency_threshold = float(t)
            else:
                print('t is None')
                ensemble_settings.apolar_frequency_threshold = None
                ensemble_settings.polar_frequency_threshold = None
            ensemble = EnsembleResult(hs_results_list=hotspot_list,
                                      ensemble_id=f"{ens_name}_{t}",
                                      settings=ensemble_settings)
            ensemble.make_ensemble_maps(save_grid_ensembles=True)
            ensemble_hs_result = ensemble.ensemble_hotspot_result
            if t is not None:
                ens_path = Path(ens_df_path.parent, f'ensemble_maps_threshold_{t}')
            else:
                ens_path = Path(ens_df_path.parent, f'ensemble_maps_threshold_{ensemble_settings.combine_mode}_zeros')

            try:
                ens_df = ensemble_cluster_summary(ensemble, save_dir=ens_path)
                ens_df.to_csv(str(Path(ens_path, "ensemble_cluster_summary.csv")))

            except ValueError as ve:
                print(ve)

            with HotspotWriter(str(ens_path), grid_extension=".ccp4", zip_results=True) as w:
                w.write(ensemble_hs_result)

