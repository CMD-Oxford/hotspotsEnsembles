from pathlib import Path
import pandas as pd
import numpy as np
import requests
import json
from ccdc.protein import Protein
from ccdc.molecule import Molecule
from ccdc.io import MoleculeReader, MoleculeWriter
from hotspots.hs_io import HotspotReader
import xml.etree.ElementTree as ET 
import matplotlib.pyplot as plt

def get_pdb_protein_structure(pdb_id):
    res = requests.get(f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb")
    return res.text

def get_compound_pdb_id(chembl_compound_id):
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
      "return_type": "mol_definition"
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

def get_pdb_compound_structure(ccdc_prot, pdb_prot_id, compound_pdb_id):
    lident = [l.identifier for  l in ccdc_prot.ligands if compound_pdb_id in l.identifier]
    lig_dict = {}
    for lid in lident:
        chain = lid.split(':')[0]
        authseq = lid.split(':')[1][3:]
        res = requests.get(f"https://models.rcsb.org/v1/{pdb_prot_id}/ligand?auth_asym_id={chain}&auth_seq_id={authseq}&encoding=sdf")
        lig_dict[f"{pdb_prot_id}_{compound_pdb_id}_{chain}"] = res.text.split('>')[0]
    return lig_dict

def align_to_ref_prot(ref_pdb_path, compound_df, save_dir, chembl_target_ids):
    ref_structure = Protein.from_file(str(ref_pdb_path))
    sup = Protein.ChainSuperposition()
    sets = sup.Settings()
    sets.superposition_atoms = "CALPHA"
    tar_binding_site = Protein.BindingSiteFromMolecule(ref_structure,ref_structure.ligands[0], 6.0)
    
    structure_columns = [x for x in compound_df.columns if 'struct' in x]
    
    for compound_chembl in compound_df['compound_chembl_id'].values:
        compound_save_dir = Path(save_dir, compound_chembl)
    
        if not compound_save_dir.exists():
            compound_save_dir.mkdir()

        target_rows = compound_df[compound_df['compound_chembl_id'] == compound_chembl]
        pdb_structs = {col:target_rows[col].values for col in structure_columns}

        for target, vals in pdb_structs.items():
            if not pd.isna(vals[0]):
                strs = [x.split("'") for x in vals][0]
                if len(strs) > 1:
                    strs = [strs[z] for z in range(1, len(strs), 2)]

                for pdb_code in strs:
                    pdb_save_path = Path(compound_save_dir, f"{pdb_code}_{target.split('_struct')[0]}.pdb")

                    if not pdb_save_path.exists():
                        pdb_text = get_pdb_protein_structure(pdb_code)
                        pdb_save_path.write_text(pdb_text)
                    cur_prot = Protein.from_file(str(pdb_save_path))

                    comp_name = get_compound_pdb_id(compound_chembl)[0]
                    sdf_comps = get_pdb_compound_structure(cur_prot, pdb_code, comp_name)

                    for sdf_name, sdf_text in sdf_comps.items():
                        ccdc_prot = cur_prot.copy()
                        comp_save_path = Path(compound_save_dir, f"{target}_{sdf_name}.sdf")
                        if not comp_save_path.exists():
                            comp_save_path.write_text(sdf_text)

                        chain_name = sdf_name.split('_')[-1].split('.sdf')[0]
                        tar_chain = [ch for ch in ccdc_prot.chains if ch.identifier==chain_name][0]

                        rmsd, transformation = sup.superpose(ref_structure.chains[0], tar_chain, tar_binding_site)
                        aligned_pdb_path = Path(str(pdb_save_path).replace(".pdb", f"_{chain_name}_aligned.pdb"))

                        with MoleculeWriter(str(aligned_pdb_path)) as protwr:
                            protwr.write(ccdc_prot)

                        try:
                            lig_mols = list(MoleculeReader(str(comp_save_path)))
                            for lig_mol in lig_mols:
                                lig_mol.transform(transformation)
                        except:
                            print(sdf_name)
                            print([l.identifier.split(':')[1][:3] for l in ccdc_prot.ligands])
                            lig_mols = [l for l in ccdc_prot.ligands if (l.identifier.split(':')[1][:3] == sdf_name.split("_")[-2] and l.identifier.split(':')[0] == chain_name)]
                            print(lig_mols)
    #                         for lig_mol in lig_mols:
    #                             lig_mol.transform(transformation)

                        new_comp_save_path = str(comp_save_path).replace(".sdf", "_aligned.sdf")

                        with MoleculeWriter(new_comp_save_path) as molwr:
                            for lig_mol in lig_mols:
                                molwr.write(lig_mol)

# if __name __ == "__main__":
    
chembl_targets = {'BRD9': 'CHEMBL3108640',
         'BRD1': 'CHEMBL2176774',
         'BRPF1': 'CHEMBL3132741',
         'BRD7': 'CHEMBL3085622',
         'BRD2': 'CHEMBL1293289',
         'BRD4': 'CHEMBL1163125'}

compounds = pd.read_csv("../results/chembl_compounds/chembl_final_two_actives.csv", index_col=0)

save_path = Path("../results/aligned_chembl_compounds")
if not save_path.exists(): save_path.mkdir()

ref_struct_path = Path("../../../case_studies/bromodomains_protoss_csd2021/BRD1/pdb_files/5AMF_7.pdb")
    
align_to_ref_prot(ref_pdb_path=ref_struct_path,
                  compound_df=compounds,
                  save_dir=save_path,
                  chembl_target_ids=chembl_targets)
