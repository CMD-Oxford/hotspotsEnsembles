# Generic imports
from pathlib import Path
from shutil import make_archive, rmtree, unpack_archive
import argparse
import luigi
import tempfile

#CCDC imports
from ccdc.io import MoleculeWriter
from ccdc.protein import Protein

# Hotspots imports
from hotspots.atomic_hotspot_calculation import _AtomicHotspotResult
from hotspots.calculation import Runner
from hotspots.grid_extension import Grid
from hotspots import hs_io

tempfile.tempdir =  "/home/jin76872/Desktop/Mih/Data/tmp_superstar_ghecom"

class RunHotspots:

    def __init__(self, pdb_file, number_rotations,chains=None, protonated=False, charged=False, spheres=False):
        self.in_file = Path(pdb_file)
        self.number_rotations = number_rotations
        self.protonated=protonated
        self.charged = charged
        self.chains = chains
        self.spheres=spheres
        if self.spheres:
            self.out_dir = Path((self.in_file).parent, "fullsize_hotspots_{}_spheres".format(self.number_rotations))
        else:
            self.out_dir = Path((self.in_file).parent, "fullsize_hotspots_{}".format(self.number_rotations))

        if not self.out_dir.exists(): self.out_dir.mkdir()

    def prepare_protein(self):

        prot = Protein.from_file(str(self.in_file))
        prot.identifier = self.in_file.stem
        if self.chains:
            for ch in prot.chains:
                if ch.identifier not in self.chains:
                    prot.remove_chain(ch.identifier)
        if not self.protonated:
            prot.add_hydrogens()
        #prot.remove_all_waters()
        prot.detect_ligand_bonds()

        for ligand in prot.ligands:
            prot.remove_ligand(ligand.identifier)
        prot.remove_all_metals()
        with MoleculeWriter(str(Path(self.out_dir, "hotspots_input_protein.pdb"))) as wr:
            wr.write(prot)

        return prot

    def _save_superstar_grids(self, hs_runner):
        """
        Saves and Zips the SuperStar grids from the hotspots.calculation.Runner.
        Need to also save the buriedness grid.
        :param hs_runner: 
        :return: 
        """
        ss_dir = Path(self.out_dir.parent, "superstar_grids")
        if not ss_dir.exists(): ss_dir.mkdir()

        # Need to save the buriedness grid as well as the superstar grids
        hs_runner.buriedness.write(str(Path(ss_dir, "buriedness.ccp4").resolve()))
        for s in hs_runner.superstar_grids:
            s.grid.write(str(Path(ss_dir, "superstar_{}.ccp4".format(s.identifier)).resolve()))
        make_archive(ss_dir, 'zip', ss_dir)
        rmtree(ss_dir)

    def create_atomic_hotspots(self, superstar_grids_dir):
        """
        
        :param superstar_grids_dir: path to where the superstar grids are stored
        :return: 
        """
        atomic_hotspots = []

        # Read in the SuperStar and Buriedness info
        probes = ['donor', 'acceptor', 'apolar', 'positive', 'negative']
        b_grid = Grid.from_file(str(Path(superstar_grids_dir, 'buriedness.ccp4').resolve()))

        for p in probes:
            g_path = Path(superstar_grids_dir, f'superstar_{p}.ccp4')
            if g_path.exists():
                print(f" found grid for probe of type {p}")
                p_grid = Grid.from_file(str(g_path.resolve()))
                ahs = _AtomicHotspotResult(identifier=p, grid=p_grid, buriedness=b_grid)
                atomic_hotspots.append(ahs)
            else:
                continue

        return atomic_hotspots

    def run_hotspot_calculation(self,  method="ghecom"):
        """
        Runs the hotspots calculation on the specified PDB structure
        :return: 
        """
        h = Runner()
        settings = h.Settings()
        settings.nrotations = self.number_rotations
        settings.apolar_translation_threshold = 15
        settings.polar_translation_threshold = 15
        settings.sphere_maps = self.spheres


        # Check if SuperStar and Ghecom have already been run.
        super_archive_path = Path(self.out_dir.parent, "superstar_grids.zip")

        if super_archive_path.exists():
            super_tmp_path = Path(self.out_dir.parent, super_archive_path.stem)

            if not super_tmp_path.exists(): super_tmp_path.mkdir()
            unpack_archive(super_archive_path, super_tmp_path, 'zip')
            b_grid = Grid.from_file(str(Path(super_tmp_path, 'buriedness.ccp4').resolve()))

            result = h.from_superstar(protein=self.prepare_protein(),
                                            superstar_grids=self.create_atomic_hotspots(super_tmp_path),
                                            buriedness=b_grid,
                                            charged_probes=self.charged,
                                            settings=settings,
                                            clear_tmp=True)
            rmtree(super_tmp_path)

        else:

            result = h.from_protein(protein=self.prepare_protein(),
                                    charged_probes=self.charged,
                                    probe_size=7,
                                    buriedness_method=method,
                                    cavities=None,
                                    nprocesses=1,
                                    settings=settings)

            # Save and zip the SuperStar Grids:
            self._save_superstar_grids(h)

        # Save and zip the Results
        with hs_io.HotspotWriter(str(self.out_dir.resolve()), grid_extension=".ccp4", zip_results=True) as writer:
            writer.write(result)

        print(f"out_file: {str(Path(self.out_dir, 'out.zip').resolve())}")

        return Path(self.out_dir, 'out.zip')

class lRunner(luigi.Task):

    in_file = luigi.parameter.Parameter()
    nrot = luigi.parameter.IntParameter()
    chains = luigi.parameter.Parameter()
    charged = luigi.parameter.BoolParameter()
    protonated = luigi.parameter.BoolParameter()
    spheres = luigi.parameter.BoolParameter()

    def run(self):
        RunHotspots(self.in_file, number_rotations=self.nrot, chains=self.chains, charged=self.charged, protonated=self.protonated, spheres=self.spheres).run_hotspot_calculation()

    def output(self):
        parent_path = Path(self.in_file).parent

        if self.spheres == True:
            tar = Path(parent_path, "fullsize_hotspots_{}_spheres".format(self.nrot), "out.zip")
        else:
            tar = Path(parent_path, "fullsize_hotspots_{}".format(self.nrot), "out.zip")

        return luigi.LocalTarget(str(tar.resolve()))


class ParalleliselRunner(luigi.WrapperTask):

    in_list = luigi.parameter.ListParameter()
    nrot = luigi.parameter.IntParameter()
    chains = luigi.parameter.Parameter()
    charged = luigi.parameter.BoolParameter()
    protonated = luigi.parameter.BoolParameter()
    spheres = luigi.parameter.BoolParameter()

    def requires(self):
        for k in self.in_list:
            yield lRunner(k, self.nrot, self.chains,  self.charged, self.protonated, self.spheres)

def run_parallel_hotspot_jobs(str_paths, hs_rotations=3000, hs_chains= None, hs_charged=False, hs_protonated=False, hs_spheres=False):
    luigi.build(
        [ParalleliselRunner(str_paths, hs_rotations, hs_chains, hs_charged, hs_protonated, hs_spheres)],
        local_scheduler=True,
        workers=25)

    if hs_spheres == True:
        return [Path(Path(str_p).parent, "fullsize_hotspots_{}_spheres".format(hs_rotations), "out.zip") for str_p in str_paths]
    else:
        return [Path(Path(str_p).parent,  "fullsize_hotspots_{}".format(hs_rotations), "out.zip") for str_p in str_paths]



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_path", dest="pdb_file", help="Path (str) to the input pdb file")

    parser.add_argument("--number_rotations", dest="number_rotations", type=int,
                        help="Number (int) of probe rotations used by Hotspots algorithm")

    parser.add_argument("-p", "--protonated", dest="protonated", default=False, type=bool,
                        help="Whether the input protein has been protonated. If False (default), the CCDC API protonation is used")

    parser.add_argument("-s", "--sphere_maps", dest="spheres", default=False, type=bool,
                        help="Whether to use the sphere_maps option in the hotspots alogorithm. Defaults to False")

    parser.add_argument("-c", "--charged", dest="charged", default=False, type=bool,
                        help="Whether to use charged probes (defaults to False")


    args = parser.parse_args()
    print("Arguments", args)

    rh = RunHotspots(pdb_file=args.pdb_file,
                     number_rotations=args.number_rotations,
                     protonated=args.protonated,
                     charged=args.charged,
                     spheres=args.spheres)
    rh.run_hotspot_calculation()

if __name__ == "__main__":
    main()

