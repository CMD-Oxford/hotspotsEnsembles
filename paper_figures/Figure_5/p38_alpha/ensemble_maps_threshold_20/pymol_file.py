
from os.path import join
import tempfile
import tkinter as tk
import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

wd = None
finish_launching()

dirpath = tempfile.mkdtemp()
zip_dir = "out.zip"
wd = os.getcwd()
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

os.chdir(dirpath)
cmd.set("normalize_ccp4_maps", "0")
cmd.load("hotspot/donor.ccp4", "donor_hotspot")
cmd.isosurface(name="surface_donor_hotspot", map="donor_hotspot", level="5")

cmd.color("blue", "surface_donor_hotspot")
cmd.set("transparency", 0.2, "surface_donor_hotspot")
cmd.load("hotspot/acceptor.ccp4", "acceptor_hotspot")
cmd.isosurface(name="surface_acceptor_hotspot", map="acceptor_hotspot", level="5")

cmd.color("red", "surface_acceptor_hotspot")
cmd.set("transparency", 0.2, "surface_acceptor_hotspot")
cmd.load("hotspot/apolar.ccp4", "apolar_hotspot")
cmd.isosurface(name="surface_apolar_hotspot", map="apolar_hotspot", level="5")

cmd.color("yellow", "surface_apolar_hotspot")
cmd.set("transparency", 0.2, "surface_apolar_hotspot")
cmd.group("hotspot", members="donor_hotspot")
cmd.group("hotspot", members="surface_donor_hotspot")
cmd.group("hotspot", members="acceptor_hotspot")
cmd.group("hotspot", members="surface_acceptor_hotspot")
cmd.group("hotspot", members="apolar_hotspot")
cmd.group("hotspot", members="surface_apolar_hotspot")
cmd.group("hotspot", members="buriedness_hotspot")
cmd.group("hotspot", members="surface_buriedness_hotspot")
cmd.pseudoatom(object="PS_donor_hotspot_0", pos=(1.0, 18.0, 47.0), color=(1, 1, 1), label=18.0)

cmd.pseudoatom(object="PS_donor_hotspot_1", pos=(1.5, 22.0, 42.0), color=(1, 1, 1), label=17.7)

cmd.pseudoatom(object="PS_donor_hotspot_2", pos=(4.0, 18.5, 39.0), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_donor_hotspot_3", pos=(5.5, 22.5, 36.0), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_donor_hotspot_4", pos=(2.0, 13.5, 34.5), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_donor_hotspot_5", pos=(6.5, 19.0, 40.0), color=(1, 1, 1), label=14.7)

cmd.pseudoatom(object="PS_donor_hotspot_6", pos=(-1.5, 15.5, 40.0), color=(1, 1, 1), label=12.6)

cmd.pseudoatom(object="PS_donor_hotspot_7", pos=(-1.0, 13.5, 40.0), color=(1, 1, 1), label=12.3)

cmd.pseudoatom(object="PS_donor_hotspot_8", pos=(-1.0, 14.0, 44.5), color=(1, 1, 1), label=12.2)

cmd.pseudoatom(object="PS_donor_hotspot_9", pos=(0.5, 14.5, 37.5), color=(1, 1, 1), label=11.8)

cmd.pseudoatom(object="PS_donor_hotspot_10", pos=(6.0, 11.0, 37.0), color=(1, 1, 1), label=10.4)

cmd.pseudoatom(object="PS_donor_hotspot_11", pos=(-2.5, 27.5, 47.5), color=(1, 1, 1), label=10.4)

cmd.pseudoatom(object="PS_donor_hotspot_12", pos=(-3.0, 13.5, 45.0), color=(1, 1, 1), label=10.0)

cmd.group("label_donor_hotspot", members="PS_donor_hotspot_0")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_0")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_1")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_1")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_2")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_2")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_3")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_3")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_4")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_4")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_5")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_5")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_6")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_6")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_7")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_7")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_8")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_8")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_9")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_9")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_10")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_10")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_11")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_11")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_12")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_12")
cmd.pseudoatom(object="PS_acceptor_hotspot_0", pos=(1.5, 20.5, 44.5), color=(1, 1, 1), label=18.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_1", pos=(6.0, 22.5, 37.5), color=(1, 1, 1), label=14.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_2", pos=(0.0, 13.5, 44.5), color=(1, 1, 1), label=13.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_3", pos=(2.5, 15.0, 34.5), color=(1, 1, 1), label=13.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_4", pos=(9.0, 19.0, 34.0), color=(1, 1, 1), label=12.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_5", pos=(7.0, 17.5, 38.0), color=(1, 1, 1), label=12.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_6", pos=(5.5, 26.0, 39.0), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_7", pos=(5.5, 19.0, 34.0), color=(1, 1, 1), label=11.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_8", pos=(3.5, 22.5, 35.5), color=(1, 1, 1), label=6.5)

cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_0")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_0")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_1")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_1")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_2")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_2")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_3")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_3")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_4")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_4")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_5")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_5")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_6")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_6")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_7")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_7")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_8")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_8")
cmd.pseudoatom(object="PS_apolar_hotspot_0", pos=(4.0, 21.5, 40.0), color=(1, 1, 1), label=18.1)

cmd.pseudoatom(object="PS_apolar_hotspot_1", pos=(5.5, 23.5, 39.5), color=(1, 1, 1), label=17.5)

cmd.pseudoatom(object="PS_apolar_hotspot_2", pos=(6.0, 23.5, 41.0), color=(1, 1, 1), label=17.5)

cmd.pseudoatom(object="PS_apolar_hotspot_3", pos=(3.0, 19.0, 43.5), color=(1, 1, 1), label=17.1)

cmd.pseudoatom(object="PS_apolar_hotspot_4", pos=(2.0, 19.0, 41.0), color=(1, 1, 1), label=16.8)

cmd.pseudoatom(object="PS_apolar_hotspot_5", pos=(1.0, 15.5, 45.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_6", pos=(-0.5, 18.5, 41.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_apolar_hotspot_7", pos=(-3.0, 24.0, 43.0), color=(1, 1, 1), label=14.4)

cmd.pseudoatom(object="PS_apolar_hotspot_8", pos=(9.0, 13.0, 34.0), color=(1, 1, 1), label=14.2)

cmd.pseudoatom(object="PS_apolar_hotspot_9", pos=(-2.0, 15.0, 43.0), color=(1, 1, 1), label=13.3)

cmd.pseudoatom(object="PS_apolar_hotspot_10", pos=(6.0, 17.0, 34.0), color=(1, 1, 1), label=12.8)

cmd.pseudoatom(object="PS_apolar_hotspot_11", pos=(4.0, 13.5, 34.0), color=(1, 1, 1), label=12.4)

cmd.pseudoatom(object="PS_apolar_hotspot_12", pos=(5.0, 15.5, 37.0), color=(1, 1, 1), label=11.4)

cmd.pseudoatom(object="PS_apolar_hotspot_13", pos=(3.5, 12.0, 37.5), color=(1, 1, 1), label=10.6)

cmd.pseudoatom(object="PS_apolar_hotspot_14", pos=(0.5, 12.0, 39.5), color=(1, 1, 1), label=8.4)

cmd.pseudoatom(object="PS_apolar_hotspot_15", pos=(-3.0, 23.5, 49.0), color=(1, 1, 1), label=7.2)

cmd.pseudoatom(object="PS_apolar_hotspot_16", pos=(-3.0, 25.5, 49.0), color=(1, 1, 1), label=7.2)

cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_0")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_0")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_1")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_1")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_2")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_2")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_3")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_3")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_4")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_4")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_5")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_5")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_6")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_6")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_7")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_7")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_8")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_8")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_9")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_9")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_10")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_10")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_11")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_11")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_12")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_12")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_13")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_13")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_14")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_14")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_15")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_15")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_16")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_16")
cmd.group("labels_hotspot", members="label_donor_hotspot")
cmd.group("labels_hotspot", members="label_acceptor_hotspot")
cmd.group("labels_hotspot", members="label_apolar_hotspot")
cmd.load("hotspot/protein.pdb", "protein_hotspot")
cmd.show("cartoon", "protein_hotspot")
cmd.hide("line", "protein_hotspot")
cmd.show("sticks", "organic")


class IsoLevel(tk.Variable):
    def __init__(self, master, name, level):
        tk.Variable.__init__(self, master, value=level)
        self.name = name
        self.trace('w', self.callback)

    def callback(self, *args):
        cmd.isolevel(self.name, self.get())

    def increment(self, event=None, delta=0.1):
        self.set(round(float(self.get()) + delta, 2))

    def decrement(self, event=None):
        self.increment(None, -0.1)


surface_list = {'hotspot': {'fhm': ['surface_donor_hotspot', 'surface_acceptor_hotspot', 'surface_apolar_hotspot']}}
surface_max_list = {'hotspot': {'fhm': 18.3}}

top = tk.Toplevel(plugins.get_tk_root())

master = tk.Frame(top, padx=10, pady=10)
master.pack(fill="both", expand=1)

for child in list(master.children.values()):
    child.destroy()


row_counter = 0
for identifier, component_dic in surface_list.items():
    # add calculation identifier
    tk.Label(master, text=identifier).grid(row=row_counter, column=0, sticky="w")
    row_counter += 1
    
    for component_id, surfaces in component_dic.items():
        # add collection label, e.g. superstar or hotspot etc.
        tk.Label(master, text=component_id).grid(row=row_counter, column=1, sticky='w')
        row_counter += 1
        
        for i, surface in enumerate(surfaces):
            # add grid type label
            probe = surface.split("_")[-2]
            tk.Label(master, text=probe).grid(row=row_counter, column=2, sticky="w")
            
            # slider code 
            v = IsoLevel(master, surface, 5)
            e = tk.Scale(master, orient=tk.HORIZONTAL, from_=0, to=surface_max_list[identifier][component_id],
                         resolution=0.1, showvalue=0, variable=v)
            e.grid(row=row_counter, column=3, sticky="ew")

            e = tk.Entry(master, textvariable=v, width=4)
            e.grid(row=row_counter, column=4, sticky="e")
            master.columnconfigure(3, weight=1)
            row_counter += 1



cmd.bg_color("white")
if wd:
    os.chdir(wd)