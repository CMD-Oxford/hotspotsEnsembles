from hotspots.hs_utilities import Helper
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hdbscan import HDBSCAN
import collections


Coordinates = collections.namedtuple('Coordinates', ['x', 'y', 'z'])

def reduce_clusters(center_list, labels, threshold=0.5):
    """
    Given a list of centroids and labels, put any clusters below threshold distance apart in the same cluster.
    :param center_list:
    :param labels:
    :return:
    """
    cluster_dict = {l:c for l, c in zip(labels, center_list)}
    distances = []
    for c in center_list:
        to_others=[]
        coords1 = Coordinates(c[0], c[1], c[2])
        for d in center_list:
            coords2 = Coordinates(d[0], d[1], d[2])
            to_others.append(Helper.get_distance(coords1, coords2))
        distances.append(to_others)
    c_arr = np.array(distances)
    new_arr = (c_arr > 0.0) & (c_arr <= threshold)
    new_arr = np.triu(new_arr, k=0)
    maps1 = [labels[i] for i in new_arr.nonzero()[0]]
    maps2 = [labels[x] for x in new_arr.nonzero()[1]]
    mapping = {m1: m2 for m1, m2 in zip(maps1, maps2)}
    new_labels = labels.copy()
    passed_clusts = []

    for m in mapping.keys():
        for ite, n in enumerate(new_labels):
            if (m not in passed_clusts) and (m == n):
                new_labels[ite] = mapping[m]
        passed_clusts.append(m)

    return new_labels


def make_selmap_clusters(df_paths, ens_dir):
    df_list = []

    for dcd in df_paths:
        df = pd.read_csv(Path(dcd, 'ensemble_cluster_summary.csv'), index_col=0)
        vals = dcd.name.split('_')
        df['probe_type'] = [i.split('_')[0] for i in df['clust_id'].values]
        df['on_target'] = vals[1]
        df['off_target'] = vals[3]
        df[vals[4]] = vals[5]
        df[vals[6]] = vals[7]

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

def visualise_selmap_clusters(master_df_paths, distances, scores):
    colour_dict = {'acceptor': 'r',
                   'donor' : 'b',
                   'apolar': 'y'}
    for mp in master_df_paths:
        mdf = pd.read_csv(mp, index_col=0)
        image_save_dir = Path(mp.parent, 'cluster_images')
        if not image_save_dir.exists():
            image_save_dir.mkdir()
        clusts = set(mdf['cluster_cluster'])
        on_target = list(set(mdf['on_target']))[0]
        off_target = list(set(mdf['off_target']))[0]
        probe = list(set(mdf['probe_type']))[0]
        for c in clusts:
            clust_df = mdf[mdf['cluster_cluster'] == c]
            scores_dict = {score: i for i, score in enumerate(scores)}
            dists_dict = {dist: ind for ind, dist in enumerate(distances)}
            center_list = [x.split('[')[1].split(']')[0].split(',') for x in clust_df['centroid'].values]
            center_list = np.array([list(map(float, x)) for x in center_list])
            center_list = np.mean(center_list, axis=0)
            print(center_list)
            for i, r in clust_df.iterrows():
                plt.scatter(int(dists_dict[r['distance']]) + 1,
                            int(scores_dict[r['score']]) + 1,
                            marker="o",
                            s=r['clust_size'] / r['num_points_map'] * 500,
                            alpha=0.6,
                            c=colour_dict[r['probe_type']])

            ax = plt.gca()
            xticks = distances
            yticks = scores
            ax.set_xticks([0, 1, 2, 3, 4, 5, 6])
            ax.set_yticks([0, 1, 2, 3])
            ax.set_xticklabels([0.0] + xticks)
            ax.set_yticklabels([0.0] + yticks)
            ax.set_xlabel('Distance to closest off-target cluster')
            ax.set_ylabel('Cluster median cutoff (Hotspot score)')
            title_str = f"{on_target} over {off_target} {probe} cluster {c} \n centroid: x={round(center_list[0],2)}, y={round(center_list[1],2)}, z={round(center_list[2],2)}"
            ax.set_title(title_str)
            plt.savefig(Path(image_save_dir, f"{on_target}_{off_target}_{probe}_cluster_{c}.png" ))
            plt.close()

def visualise_ensmap_cluster_contribs(master_df_paths,on_target):
    # colour_dict = {'acceptor': 'Reds',
    #                'donor' : 'Blues',
    #                'apolar': 'YlOrBr'}
    colour_dict = {'acceptor': 'r',
                   'donor' : 'b',
                   'apolar': 'y'}

    for mp in master_df_paths:
        mdf = pd.read_csv(mp, index_col=0)
        image_save_dir = Path(mp.parents[1], 'contrib_images')
        if not image_save_dir.exists():
            image_save_dir.mkdir()
        contribs = [x.split('[')[1].split(']')[0].split(',') for x in mdf['contributing_structures'].values]
        contribs = list(set([a.split('_')[0].replace(("'"), "").strip() for x in contribs for a in x]))
        contribs_dict = {a:i for i,a in enumerate(contribs)}

        plt.figure(figsize=(12, 8))
        for it, row in mdf.iterrows():

            probe = row['clust_id'].split('_')[0]
            contributing_structures = row['contributing_structures'].split('[')[1].split(']')[0].split(',')
            contributing_structures = [y.split('_')[0].replace(("'"), "").strip() for y in contributing_structures]
            row_contrib_idxs = [contribs_dict[name] for name in contributing_structures]
            print(row)
            print(contributing_structures)
            contributions = [int(al) for al in row['structure_contributions'].split('[')[1].split(']')[0].split(',')]

            plt.scatter(row_contrib_idxs,
                        [it] * len(row_contrib_idxs),
                        marker="o",
                        s=[c/max(contributions)*500 for c in contributions],
                        c=colour_dict[probe],
                        alpha=0.6)
                        # c=[row['clust_median']] * len(row_contrib_idxs),
                        # cmap=colour_dict[probe])


        ax = plt.gca()
        xticks = contribs
        yticks = mdf['clust_id']
        ax.set_xticks(range(len(xticks)))
        ax.set_yticks(range(len(yticks)))
        ax.set_xticklabels(xticks, rotation=45)
        ax.set_yticklabels(yticks)
        ax.set_xlabel('PDB code')
        ax.set_ylabel('Cluster name')
        title_str = mp.parent.name
        ax.set_title(title_str)
        print(image_save_dir)
        plt.savefig(Path(image_save_dir, f"{mp.parent.name}_contributions.png" ))
        plt.close()

def ti_si_guz(summary_paths):
    guz_dict = {}
    for p in summary_paths:
        if p.parent.name.split('_')[-1]=='zeros':
            continue
        threshold = int(p.parent.name.split('_')[-1])
        print(threshold)
        pdf = pd.read_csv(p)
        probes = [x.split("_")[0] for x in pdf['clust_id']]
        tar_idxs = []
        for i, p in enumerate(probes):
            if p == 'donor' or p == 'acceptor':
                tar_idxs.append(i)
        polar_df = pdf.iloc[tar_idxs, :]
        guz_dict[threshold] = polar_df
    return guz_dict

if __name__ == "__main__":
    np.random.seed=3
    # ens_dir = Path('../case_studies/select_maps_param_run/')
    # dist_clust_dfs = list(ens_dir.glob('selectivity_*'))
    # table_paths = make_selmap_clusters(dist_clust_dfs)
    # print(table_paths)
    # scores_list = [5.0, 10.0, 15.0]
    # dists_list = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    # visualise_selmap_clusters(table_paths, distances=dists_list, scores=scores_list)
    info_dict = {
    'BRD1': Path(f'/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/case_studies/bromodomains_protoss_csd2021/BRD1'),
    # 'BRPF1': Path(f'/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/case_studies/bromodomains_protoss_csd2021/BRPF1'),
    # 'CK2alpha': Path(f'/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/case_studies/KLIFS_kinases_2021/KLIFS_CK2alpha_01_10_2020'),
    # 'PIM1': Path(f'/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/case_studies/KLIFS_kinases_2021/KLIFS_PIM1_01_10_2020'),
    # 'p38a' : Path('/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/case_studies/KLIFS_kinases_2021/KLIFS_p38a_fragments_09_02_2021'),
    # 'ERK2' : Path('/home/jin76872/Desktop/Mih/Data/JCIM_selectivity_paper/case_studies/KLIFS_kinases_2021/KLIFS_ERK2_fragments_09_02_2021')
    }
    for ens_name, ensemble_dir in info_dict.items():
        summary_paths = list(ensemble_dir.glob("ensemble_maps_threshold_*/ensemble_cluster_summary.csv"))
        # visualise_ensmap_cluster_contribs(summary_paths, ens_name)
        gdict = ti_si_guz(summary_paths)