from pathlib import Path
import numpy as np
import pandas as pd
from utils import get_clusters_centre_mass
from hotspots.hs_io import  HotspotReader, HotspotWriter
from hotspots.hs_ensembles import SelectivityResult
from hotspots.grid_extension import Grid, _GridEnsemble

def as_grid(origin_coords, far_corner_coords, array, spacing=0.5):
    """
    Given an array, outputs a grid with the dimensions of the GridEnsemble
    :param array: 3D numpy array, usually containing processed ensemble data
    :return: a :class: 'ccdc.utilities.Grid' instance
    """
    # Initialise the Grid
    grid = Grid(origin=origin_coords,
                far_corner=far_corner_coords,
                spacing=spacing,
                default=0.0,
                _grid=None)

    # Get the nonzero indices and values of the array
    nonz = array.nonzero()
    values = array[nonz]
    # Get indices per value
    as_triads = zip(*nonz)

    # Fill in the grid
    for (i, j, k), v in zip(as_triads, values):
        grid._grid.set_value(int(i), int(j), int(k), v)

    return grid

def ensemble_cluster_summary(hs_res, save_dir, polar_min_clust=3, apolar_min_clust=20):
    clust_ids = []
    clust_size = []
    total_points_map = []
    median_clust_value = []
    clust_centroid = []


    for probe in hs_res.super_grids.keys():
        # Get the contributing grids for all the clusters
        try:
            # Convert into numpy array:
            probe_array = hs_res.super_grids[probe].get_array()

            # Find the hotspot features using HDBSCAN
            if probe == 'apolar':
                ensemble_probe_clusters = _GridEnsemble.HDBSCAN_cluster(probe_array,
                                                                        min_cluster_size=apolar_min_clust,
                                                                        allow_single_cluster=True)
            else:
                ensemble_probe_clusters = _GridEnsemble.HDBSCAN_cluster(probe_array,
                                                                        min_cluster_size=polar_min_clust,
                                                                        allow_single_cluster=True)


            probe_clust_grid = as_grid(hs_res.super_grids[probe].bounding_box[0],
                                       hs_res.super_grids[probe].bounding_box[1],
                                       ensemble_probe_clusters,
                                       hs_res.super_grids[probe].spacing)
            probe_clust_grid.write(str(Path(save_dir, f"{probe}_clusters_ranges.ccp4")))
            coords = get_clusters_centre_mass(ensemble_probe_clusters, probe_array)
            clust_medians = []

            for c in coords.keys():
                clust_ids.append(f"{probe}_{c}")
                clust_size.append(len(ensemble_probe_clusters[ensemble_probe_clusters == c]))
                total_points_map.append(probe_clust_grid.count_grid())
                clust_centroid.append(coords[c])
                clust_medians.append(np.median(probe_array[ensemble_probe_clusters==c]))
        except ValueError as ve:
            print(ve)


    probe_df = pd.DataFrame()
    probe_df['clust_id'] = clust_ids
    probe_df['clust_size'] = clust_size
    probe_df['clust_median'] = clust_medians
    probe_df['num_points_map'] = total_points_map
    probe_df['centroid'] = clust_centroid

    return probe_df

if __name__ == "__main__":

    info_dict = {
    'BRD1': Path(f'../ensembles/BRD1/ensemble_maps/out.zip'),
    'BRPF1': Path(f'../ensembles/BRPF1/ensemble_maps/out.zip'),
    'BRD7': Path(f'../ensembles/BRD7/ensemble_maps/out.zip'),
    'BRD2': Path(f'../ensembles/BRD2/ensemble_maps/out.zip'),
    'BRD4' : Path('../ensembles/BRD4/ensemble_maps/out.zip'),
    'BRD9' : Path('../ensembles/BRD9/ensemble_maps/out.zip')
    }
    
    on_target = 'BRD1'
    for off_target in info_dict.keys():
        if on_target != off_target:
            on_target_path = info_dict[on_target]
            off_target_path = info_dict[off_target]
            with  HotspotReader(str(on_target_path)) as otr:
                on_target_ens = otr.read()

            with HotspotReader(str(off_target_path)) as oftr:
                off_target_ens = oftr.read()

            for c in [10.0]:
                for d in [1.5]:

                    selectivity_settings = SelectivityResult.Settings()
                    selectivity_settings.cluster_distance_cutoff = d
                    selectivity_settings.minimal_cluster_score = c

                    on_over_off = SelectivityResult(target_result=on_target_ens,
                                                    other_result=off_target_ens,
                                                    settings=selectivity_settings
                                                    )
                    on_over_off.make_selectivity_maps()

                    bromo_dir = Path('../results/select_maps_all')
                    if not bromo_dir.exists(): bromo_dir.mkdir()

                    bc_path = Path(bromo_dir, f"selectivity_{on_target}_over_{off_target}_score_{c}_distance_{d}")

                    with HotspotWriter(str(bc_path),grid_extension=".ccp4", zip_results=False) as w:
                        w.write(on_over_off.selectivity_result)

#                         bc_df = ensemble_cluster_summary(on_over_off.selectivity_result, save_dir=bc_path)
#                         bc_df.to_csv(Path(bc_path, "ensemble_cluster_summary.csv"))