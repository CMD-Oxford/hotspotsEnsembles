from pathlib import Path
import pandas as pd
from ccdc.molecule import Molecule
from ccdc.protein import Protein
from ccdc.io import MoleculeWriter
from hotspots.wrapper_siena import Ensemble, Search
from rdkit.Chem.Scaffolds import MurckoScaffold

from utils import find_bs_ligands, get_pdb_resolution, get_subset

def query_Siena(pdb_code, ligand_name, out_dir=""):
    """
    Use Pete's script to query the SIENA restful API; change the search parameters as needed.
    :param pdb_code:
    :param ligand_name:
    :param chain:
    :return:
    """
    ssettings = Search.Settings()
    ssettings.data["siena"]["pdbCode"] = pdb_code
    ssettings.data["siena"]["ligand"] = ligand_name
    ssettings.data["siena"]["minimalSiteIdentity"] = '1.0'
    ssettings.data["siena"]["resolution"] = '2.5'
    #ssettings.data["siena"]["unique_ligands"] = 'true'
    ssettings.data["siena"]["holo_only"] = ''
    ssettings.data["siena"]["complete_residues_only"] = 'true'
    print(ssettings.data["siena"])
    searcher = Search(settings=ssettings)
    results_url = searcher._post()
    ensemble = Ensemble(searcher._get(results_url))

    ensemble.save(out_dir=out_dir)

class SienaQuery():
    def __init__(self, target_name, pdb_code, ligand_name, siena_dir=''):
        self.ensemble_name = target_name
        self.query_pdb_code = pdb_code
        self.query_ligand_name = ligand_name
        self.siena_dir = siena_dir
        self.extracted_ligand_dir = Path(self.siena_dir, 'extracted_ligands')
        if not self.extracted_ligand_dir.exists(): self.extracted_ligand_dir.mkdir()

        print(self.query_pdb_code)
        self.siena_resultstat = Path(self.siena_dir, "result_table", "resultStatistic.csv")

        if not self.siena_resultstat.exists():
            query_Siena(self.query_pdb_code, self.query_ligand_name, self.siena_dir)

        self.aa_dic = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE", "K" : "LYS",
                       "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
                       "W": "TRP", "Y": "TYR"}

    def read_alignment(self):
        """
        Reformats the SIENA output file in a way that can be easily read by pandas. Outputs the alignment as a dataframe.
        :return: pandas.dataframe object
        """
        # Want to avoid modifying the original output file, so we create a copy
        new_path = Path(self.siena_dir, "pd_readable_alignment.csv")
        self.alignment_file = new_path

        with open(Path(self.siena_dir, "alignment", "alignment.txt"), "r") as input:
            with open(self.alignment_file, "w") as output:
                l = input.readlines()
                # The first line says "Alignment:" - we want to get rid of it
                l_out = l[1:]
                [output.write(l_o.replace(" ", "")) for l_o in l_out]

        ali_df = pd.read_csv(self.alignment_file,  sep="|")
        # removes the "Nan" column created at the end of the dataframe by the final separator.
        ali_df = ali_df.iloc[:, :-1]

        return ali_df

    def _get_aa_from_alignment(self, s):
        """
        Each cell of the alignment_df is formatted "<1 letter amino acid code>-<residue number>-<chain>"
        This function returns the amino acid, converted to a 3-letter code (for compatibility with CSD API).
        :param str s: contents of cell in the alignment dataframe. If not correct type of cell, returns nothing.
        :return: str or None
        """
        # Checks that the string is in the format described above.
        if len(s.split("-")) > 2:
            # get the 1 letter aa code from the string
            aa_letter = s.split("-")[0].strip()
            # Convert to 3-letter aa code
            aa_full = self.aa_dic[aa_letter]
            return aa_full
        else:
            return

    def get_chains(self, row):
        """
        The SIENA output also has a "chains" section, but this way we get the chains involved in the binding site only
        :param row: a row in the alignment dataframe
        :return: list
        """
        # row.values[0] is the PDB code
        vals = row[1:]
        # get a list of chains for each amino acid in the binding site
        allchains = [s.split("-")[2].strip() for s in vals]
        chains = list(set(allchains))
        return chains

    def get_binding_site_residues(self, row):
        """
        Converts the binding site information in the SIENA output to a format that can be passed to the CSD API
        :param row: a row in the alignment dataframe
        :return:
        """
        # row.values[0] is the PDB code
        vals = row[1:]
        # get the amino acid position in the structure
        pos = [s.split("-")[1].strip() for s in vals]
        # get the amino acids and convert to 3-letter code
        aas = [self._get_aa_from_alignment(s) for s in vals]
        # get the chain for each amino acid
        chains = [s.split("-")[2].strip() for s in vals]
        # format
        bs_res = ["{0}:{1}{2}".format(chains[i], aas[i], pos[i]) for i in range(len(vals))]

        return bs_res

    def get_binding_site_ligands(self, row, fname):
        protein = Protein.from_file(str(fname))
        print('loading protein')
        bs_residues = self.get_binding_site_residues(row)
        bs_ligands = find_bs_ligands(protein, bs_residues, threshold=4.0)
        lig_path = Path(self.extracted_ligand_dir, str(fname.stem) +'.sdf')

        if len(bs_ligands)>0:
            if len(bs_ligands)>1:
                new_mol = Molecule()
                for b in bs_ligands:
                    new_mol.add_molecule(b)
                bs_ligands=[new_mol]

            with MoleculeWriter(lig_path) as wr:
                for b in bs_ligands:
                    wr.write(b)

            return bs_ligands[0], lig_path

        else:
            return None

    def make_row_dic(self, idx, row):
        """
        Compiles each row of the dataframe that gets passed to the EnsembleResult.
        :param int idx: row index
        :param pandas.Series row:
        :return: dictionary
        """
        try:
            chains = self.get_chains(row)
            fname = Path(self.siena_dir, "pdb_files",  "{}_{}.pdb".format(row[0], idx+1))
            #print(self.get_binding_site_ligands(row, fname))
            row_dic = {"Filename": fname,
                       "ID": "{}-{}".format(row[0], "".join(chains)),
                       "PDB ID": row[0],
                       "Chains": chains,
                       "Binding site": str(self.get_binding_site_residues(row))}
            return row_dic
        except AttributeError:
            return


    def get_ensemble_protein_df(self):
        """
        Creates the dataframe that holds the ensemble information and will be propagated to the EnsembleResult.
        :return:
        """
        ali_df = self.read_alignment()
        #ali_df = ali_df.iloc[:, :-1]

        # Create a dataframe with the column names needed by EnsembleResult
        new_rows = ["ID", "Filename", "PDB ID", "Chains", "Binding site"]
        a = pd.DataFrame(columns=new_rows, dtype="object")
        ali_df_vals = ali_df.values
        # Fill dataframe
        for idx, r in enumerate(ali_df_vals):
            print(idx, r)
            a.loc[len(a)] = self.make_row_dic(idx, r)

        # All structures should have the same ensemble ID (e.g. target name)
        a["Ensemble ID"] = self.ensemble_name
        lig_tups = [self.get_binding_site_ligands(row, Path(self.siena_dir, "pdb_files",  "{}_{}.pdb".format(row[0], idx+1))) for idx, row in enumerate(ali_df_vals)]
        a["Binding site ligands"] = [a[1] if a is not None else None for a in lig_tups]
        a["Ligand MW"]  = [a[0].molecular_weight if a is not None else None for a in lig_tups]
        a["Ligand SMILES"] = [a[0].smiles if a is not None else None for a in lig_tups]
        try:
            a["Ligand Murcko"] =[MurckoScaffold.MurckoScaffoldSmilesFromSmiles(l[0].smiles) if l is not None  else None for l in lig_tups]
        except Exception as e:
            a["Ligand Murcko"] = str(e)
        a['Resolution'] = [get_pdb_resolution(pdb_file) for pdb_file in a['Filename'].values]

        # Get the information on the reference structure to pass on to EnsembleResult
        a.to_csv(str(Path(self.siena_dir, "ensemble_data.csv")))

        return a

