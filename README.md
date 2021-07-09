# hotspotsEnsembles

Scripts used to generate the data in the Case Studies Section of the publication Fragment Hotspot Mapping to Identify Selectivity-Determining Regions between Related Proteins.
The starting case study files and resulting hotspot maps are supplied as a zip file in the Supplementary Information of the publication.

The code in this repository uses the CSD Python API and the Hotspots API.

CSD Python API:
https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/

Hotspots API, including installaton instructions:
https://github.com/prcurran/hotspots

In addition, the scripts/run_hotspots_job.py depends on the Luigi module to manage the hotspot runs
https://luigi.readthedocs.io/en/stable/

To visualise the results, only Pymol is needed. 

The paper_figures folder contains Pymol sessions with loaded views of the paper figures.

Each indiviudal hotspot map folder also has a 'pymol_file.py', which can be used to view those particular maps (housed in 'out/hotspot' or 'out.zip'.
To load, open Pymol, then make sure the normalize_ccp4_maps option is set to 'off'/0 (set normalize_ccp4_maps, 0). Then, load the file in Pymol to view the maps interactively.

Map thresholds of 14 and above are generally used for the ensemble and individual hotspot maps, unless specified otherwise in the publication.
Map threhsholds of 10 and above are generally used for the selectivity maps, unless specified otherwise in the publication.
