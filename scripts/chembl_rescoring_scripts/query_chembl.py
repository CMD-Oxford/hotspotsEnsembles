from pathlib import Path
import requests
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET 
import json

# Get all Chembl activities per target
chembl_ids = {'BRD9': 'CHEMBL3108640',
         'BRD1': 'CHEMBL2176774',
         'BRPF1': 'CHEMBL3132741',
         'BRD7': 'CHEMBL3085622',
         'BRD2': 'CHEMBL1293289',
         'BRD4': 'CHEMBL1163125'}

def query_chembl(entity_name):
    chembl_query = "https://www.ebi.ac.uk/chembl/api/data/activity/search?q="
    #chembl_query = "https://www.ebi.ac.uk/chembl/api/data/molecule/"
    res = requests.get(chembl_query + entity_name)
    tree = ET.fromstring(res.text)
    meta = list(tree)[1]
    activities = list(tree)[0]
    activities_list = []
    for a in activities:
        activities_list.append({a[i].tag:a[i].text for i in range(len(a))})
    
    while type(meta[1].text) is str:
        new_page = requests.get("https://www.ebi.ac.uk" + meta[1].text)
        new_tree = ET.fromstring(new_page.text)
        new_activities = list(new_tree)[0]
        for na in new_activities:
            activities_list.append({na[i].tag:na[i].text for i in range(len(na))})
        meta = list(new_tree)[1]
    query_df = pd.DataFrame.from_records(activities_list, columns=list(activities_list[0].keys()).sort())
    return query_df

def remove_second_bromodomain_assays(sub_df):
    to_drop = []
    for ai, row in sub_df.iterrows():
        if 'bromodomain 2' in row['assay_description']:
            to_drop.append(ai)
        elif 'bromodomain2' in row['assay_description']:
            to_drop.append(ai)
        elif 'bromodomain 1/2' in row['assay_description']:
            to_drop.append(ai)
        elif 'BD2' in row['assay_description']:
            to_drop.append(ai)
        elif 'BRPF1A' in row['assay_description']:
            to_drop.append(ai)
        elif 'BD12' in row['assay_description']:
            to_drop.append(ai)
    sub_df = sub_df.drop(index=to_drop)
    return sub_df

def curate_chembl_data(chembl_df):
    curated_df = chembl_df[master_df['standard_relation'] == '=']
    curated_df = curated_df[curated_df['standard_flag'] == True]
    curated_df = curated_df[pd.isna(curated_df['assay_variant_mutation'])]
    curated_df = curated_df[pd.isna(curated_df['data_validity_comment'])]
    curated_df = curated_df[curated_df['assay_type'] == 'B']
    curated_df = curated_df[curated_df['src_id'] == 1]
    curated_df = curated_df[curated_df['standard_units'] == 'nM']
    curated_df = curated_df[pd.isna(curated_df['potential_duplicate'])]
    curated_df = curated_df[(curated_df['standard_type'] == 'IC50') | (curated_df['standard_type'] == 'Kd')]
    curated_df = remove_second_bromodomain_assays(curated_df)
    return curated_df
    
def get_chembl_quality_score(assay_id):
    chembl_query = "https://www.ebi.ac.uk/chembl/api/data/assay/" + assay_id
    res = requests.get(chembl_query)
    tree = ET.fromstring(res.text)
    tags = {a.tag: a.text for a in list(tree)}
    return float(tags['confidence_score'])

def get_target_structures(chembl_target_id):
    chembl_query = "https://www.ebi.ac.uk/chembl/api/data/target/" + chembl_target_id
    res = requests.get(chembl_query)
    tree = ET.fromstring(res.text)
    PDB_codes = {}
    t_comps = [x for x in list(tree) if x.tag == 'target_components']
    for t_comp in t_comps:
        for tc in t_comp:
            for t in tc:
                if t.tag == 'target_component_xrefs':
                    pdbs = [list(m)[0].text for m in list(t) if list(m)[2].text=='PDBe']
                    PDB_codes[[x for x in tc if x.tag == 'component_id'][0].text] = pdbs
    return PDB_codes

def get_compound_structures(chembl_compound_id):
    query_dict = {
      "query": {
        "type": "terminal",
        "service": "chemical",
        "parameters": {
          "value": "",
          "type": "descriptor",
          "descriptor_type": "InChI",
          "match_type": "graph-strict"
        }
      },
      "return_type": "entry"
    }

    chembl_query = "https://www.ebi.ac.uk/chembl/api/data/molecule/" + chembl_compound_id
    pdb_ids = []
    res = requests.get(chembl_query)
    tree = ET.fromstring(res.text)
    for e in list(tree):
        if e.tag == 'molecule_structures':
            for ms in list(e):
                if ms.tag == 'standard_inchi':
                    query_dict["query"]["parameters"]["value"] = ms.text
                    json_query = json.dumps(query_dict)
                    pdb_query = " https://search.rcsb.org/rcsbsearch/v1/query?json=" + json_query
                    pdb_res = requests.get(pdb_query)
                    print('query_result', pdb_res.text)
                    if len(pdb_res.text) > 0:
                        resp = json.loads(pdb_res.text)
                        pdb_ids.extend([resp["result_set"][i]['identifier'] for i in range(len(resp['result_set'])) if resp["result_set"][i]['score'] == 1])
    return pdb_ids

def structures_targets_compounds_overlap(with_structs, all_pdbs_direct):
    """
    param: with_structs: dict of compounds with PDB codes of all cocrystallised targets
    param: with_structs: dict of targets with PDB codes of all structures
    """
    bromodomain_structures = {key: [] for key in chembl_ids.keys()}
    with_bromo_structs = []

    for compound, pdb_codes in with_structs.items():
        for target in chembl_ids.keys():
            tar_codes = []
            for pdb in pdb_codes:
                if pdb in all_pdbs_direct[target]:
                    tar_codes.append(pdb)

            if len(tar_codes) == 0:
                bromodomain_structures[target].append(np.nan)
            elif len(tar_codes) == 1:
                bromodomain_structures[target].append(tar_codes[0])
            else:
                print(tar_codes)
                bromodomain_structures[target].append(tar_codes)
        with_bromo_structs.append(compound)
        
    struct_df  = pd.DataFrame()
    struct_df['compound_chembl_id'] = with_bromo_structs

    for target in bromodomain_structures.keys():
        struct_df[target] = bromodomain_structures[target]
    bromo_structs = struct_df.dropna(subset=list(bromodomain_structures.keys()), how='all')
    return bromo_structs

def get_activities(bromo_structs, curated_df, mode='median'):
    all_compounds_activity_dict = {}
    for cid in bromo_structs['compound_chembl_id'].values:
        activities_dict = {key: {} for key in chembl_ids.keys()}
        comp_df = curated_df[curated_df['molecule_chembl_id'] == cid]
        for target in activities_dict.keys():
            tar_comp_df = comp_df[comp_df['target_chembl_id'] == chembl_ids[target]]
            for binding_type in ['IC50', 'Kd']:
                activities_dict[target][binding_type] = []
                btar_comp_df = tar_comp_df[tar_comp_df['standard_type'] == binding_type]

                if len(btar_comp_df) == 0:
                    activities_dict[target][binding_type].append(np.nan)
                elif len(btar_comp_df) == 1:
                    activities_dict[target][binding_type].append(btar_comp_df['standard_value'].values[0])
                else:
                    if mode == 'all':
                        activities_dict[target][binding_type].append(btar_comp_df['standard_value'].values)
                    elif mode == 'median':
                        activities_dict[target][binding_type].append(np.median(btar_comp_df['standard_value'].values))
                    elif mode == 'min':
                        activities_dict[target][binding_type].append(np.min(btar_comp_df['standard_value'].values))
        all_compounds_activity_dict[cid] = activities_dict
    return all_compounds_activity_dict

def assemble_final_dataset(bromo_structs, activities_dict, targets):
    final_df = pd.DataFrame()
    comp_ids = []
    for x in bromo_structs['compound_chembl_id'].values:
        comp_ids.extend([x,x])
        
    final_df['compound_chembl_id'] = comp_ids
    for target in targets:
        target_structs = []
        target_activities = []
        target_activity_types = []
        for i, c in enumerate(bromo_structs['compound_chembl_id'].values):
            for stype, act_vals in activities_dict[c][target].items():
                target_structs.append(bromo_structs[target].values[i])
                if len(act_vals) > 1:
                    target_activities.append(act_vals)
                else:
                    target_activities.append(act_vals[0])
                target_activity_types.append(stype)
                
        final_df[f"{target}_struct"] = target_structs
        final_df[f"{target}_activity"] = target_activities
        final_df[f"{target}_standard_type"] = target_activity_types
     
    

    ite_list = []
    
    for ite, row in final_df.iterrows():
        activities = [row[x] for x in final_df.columns if 'activity' in x]
        activities = np.array([x[0] if type(x)==np.ndarray else x for x in activities])
        print(activities)

        activities = activities[np.invert(np.isnan(activities))]
        print(len(activities))
        if len(activities) > 1:
            ite_list.append(ite)
        
    two_actives_df = final_df.iloc[ite_list]
    return two_actives_df

def selratio_dict_from_final_dataset(compound_df, targets, save_dir):
    activity_types = ['IC50', 'Kd']
    on_count = {}
    for compound_chembl in set(compound_df['compound_chembl_id'].values):
        compound_save_dir = Path(save_dir, compound_chembl)   
        comp_df = compound_df[compound_df['compound_chembl_id'] == compound_chembl]

        on_count[compound_chembl] = {}
        for assay in activity_types:
            on_count[compound_chembl][assay] = {}
            sel_arr = np.zeros((len(targets), len(targets)))
            for t in range(len(targets)):
                for m in range(len(targets)):
                    if m != t:
                        sub_df = comp_df[(comp_df[f'{targets[t]}_standard_type'] == assay) & (comp_df[f'{targets[m]}_standard_type'] == assay)]
                        print(sub_df)
                        on_value = sub_df[f'{targets[t]}_activity']
                        off_value = sub_df[f'{targets[m]}_activity']

                        on_over_off = [float(y)/float(v) for v in on_value.values for y in off_value.values]
                        on_count[compound_chembl][assay][f"{targets[t]}_{targets[m]}"] = on_over_off
    jcount = json.dumps(on_count, indent=4)
    Path(save_dir, "activities_min_val.json").write_text(jcount)
    return on_count 

def selratio_df_from_dict(on_count):
    compound_sel_list = []
    activities_types_list = []
    selratio_list = []
    on_target_list = []
    off_target_list = []

    for compound in on_count.keys():
        for ac_type in on_count[compound].keys():
            for relation in on_count[compound][ac_type].keys():
                #if type(on_count[compound][ac_type][relation][0]) is list:
                if type(on_count[compound][ac_type][relation]) is list:
                    s_c_act = np.mean([np.mean(a) for a in on_count[compound][ac_type][relation]])
                else:
                    s_c_act = np.mean(on_count[compound][ac_type][relation])
                compound_sel_list.append(compound)
                activities_types_list.append(ac_type)
                selratio_list.append(s_c_act)
                on_target_list.append(relation.split('_')[0])
                off_target_list.append(relation.split('_')[1])

    selratio_df = pd.DataFrame()
    selratio_df['compound_chembl'] = compound_sel_list
    selratio_df['on_target'] = on_target_list
    selratio_df['off_target'] = off_target_list
    selratio_df['activity_type'] = activities_types_list
    selratio_df['activity_ratio'] = selratio_list
    return selratio_df

if __name__ == "__main__":
    targets = [ 'BRD1', 'BRPF1', 'BRD2', 'BRD4', 'BRD7', 'BRD9']
    save_dir = Path('../results/chembl_compounds')
    if not save_dir.exists():
        save_dir.mkdir()

#     for chembl_id in chembl_ids.keys():
#         qdf = query_chembl(chembl_ids[chembl_id])
#         qdf = qdf[qdf['target_chembl_id'] == chembl_ids[chembl_id]]
#         qdf.to_csv(str(Path(save_dir, f"{chembl_id}_all_chembl.csv")))
        
#     all_activities_paths = list(save_dir.glob("*_all_chembl.csv"))
#     master_df = pd.concat([pd.read_csv(str(aap), index_col=0) for aap in all_activities_paths])
    
#     # Curate the data available so far
#     cura_df = curate_chembl_data(master_df)
    
#     # Get the chembl quality score and retain only 9s
#     chembl_assay_ids = cura_df['assay_chembl_id'].values
#     chembl_confidence = [get_chembl_quality_score(ai) for ai in chembl_assay_ids]
#     cura_df['confidence_score'] = chembl_confidence
#     cura_df = cura_df[cura_df['confidence_score'] == 9]
#     cura_df.to_csv(str(Path(save_dir, "chembl_all_targets_curated.csv")))
    
#     #cura_df = pd.read_csv(str(Path(save_dir, "chembl_all_targets_curated.csv", index_col=0)))
    
#     compound_ids = set(cura_df['molecule_chembl_id'])
    
#     # Now get the PDB structures for all the targets available
#     all_pdbs = {}
#     for prot_name, chem_id in chembl_ids.items():
#         all_pdbs[prot_name] = get_target_structures(chem_id)
#     target_all_pdbs = {key: val[list(val.keys())[0]] for key, val in all_pdbs.items()}
#     with open(str(Path(save_dir, "targets_all_pdbs.json")), "w") as ofile:
#         ofile.write(json.dumps(target_all_pdbs, indent=4))
  
#     target_all_pdbs = json.loads(Path(save_dir, "targets_all_pdbs.json").read_text())
    
#     #Get the PDB structures for each compound
#     compound_pdbs = {cid:get_compound_structures(cid) for cid in list(compound_ids)}
#     with open(str(Path(save_dir, "compound_any_pdbs.json")), "w") as outfile:
#         outfile.write(json.dumps(compound_pdbs, indent=4))
#     compound_pdbs = json.loads(Path(save_dir, "compound_any_pdbs.json").read_text())

#     compounds_with_structs = {cid: val for cid, val in compound_pdbs.items() if len(val) > 0} # list of compounds with structures
    
#     # Now figure out the overlap between target structures and compound structures
#     compound_structures_df = structures_targets_compounds_overlap(compounds_with_structs, target_all_pdbs)
#     compound_structures_df.to_csv(str(Path(save_dir, "chembl_with_pdb_structure.csv")))
    
#     compound_activities = get_activities(compound_structures_df, cura_df, mode='min')
#     final_dataset = assemble_final_dataset(compound_structures_df, compound_activities, targets)
#     final_dataset.to_csv(str(Path(save_dir, "chembl_final_two_actives.csv")))

    final_dataset = pd.read_csv(str(Path(save_dir, "chembl_final_two_actives.csv")), index_col=0)
    selratio_dict = selratio_dict_from_final_dataset(compound_df=final_dataset, targets=targets, save_dir=save_dir)
    selratio_df = selratio_df_from_dict(on_count=selratio_dict)
    selratio_df = selratio_df.dropna(axis=0, how='any', subset=['activity_ratio'])
    selratio_df.to_csv(Path(save_dir, 'selectivities_chembl_compounds.csv'))