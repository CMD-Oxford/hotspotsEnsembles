
from os.path import join
import tempfile
import tkinter as tk
import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()

dirpath = tempfile.mkdtemp()
zip_dir = "out.zip"
wd = os.getcwd()
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

os.chdir(dirpath)
cmd.load("hotspot/acceptor.ccp4", "acceptor_hotspot")
cmd.isosurface(name="surface_acceptor_hotspot", map="acceptor_hotspot", level="5")

cmd.color("red", "surface_acceptor_hotspot")
cmd.set("transparency", 0.2, "surface_acceptor_hotspot")
cmd.load("hotspot/donor.ccp4", "donor_hotspot")
cmd.isosurface(name="surface_donor_hotspot", map="donor_hotspot", level="5")

cmd.color("blue", "surface_donor_hotspot")
cmd.set("transparency", 0.2, "surface_donor_hotspot")
cmd.load("hotspot/apolar.ccp4", "apolar_hotspot")
cmd.isosurface(name="surface_apolar_hotspot", map="apolar_hotspot", level="5")

cmd.color("yellow", "surface_apolar_hotspot")
cmd.set("transparency", 0.2, "surface_apolar_hotspot")
cmd.group("hotspot", members="acceptor_hotspot")
cmd.group("hotspot", members="surface_acceptor_hotspot")
cmd.group("hotspot", members="donor_hotspot")
cmd.group("hotspot", members="surface_donor_hotspot")
cmd.group("hotspot", members="apolar_hotspot")
cmd.group("hotspot", members="surface_apolar_hotspot")
cmd.load("hotspot/buriedness.ccp4", "buriedness_hotspot")
cmd.isosurface(name="surface_buriedness_hotspot", map="buriedness_hotspot", level="3")

cmd.color("gray", "surface_buriedness_hotspot")
cmd.set("transparency", 0.2, "surface_buriedness_hotspot")
cmd.group("hotspot", members="buriedness_hotspot")
cmd.group("hotspot", members="surface_buriedness_hotspot")
cmd.pseudoatom(object="PS_acceptor_hotspot_0", pos=(-9.5, 11.0, 1.5), color=(1, 1, 1), label=22.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_1", pos=(-9.5, 9.0, 4.0), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_2", pos=(-4.5, 8.0, 3.0), color=(1, 1, 1), label=18.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_3", pos=(8.5, -12.5, 11.5), color=(1, 1, 1), label=18.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_4", pos=(-7.0, 7.5, 3.0), color=(1, 1, 1), label=17.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_5", pos=(-6.5, 8.0, 6.5), color=(1, 1, 1), label=13.2)

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
cmd.pseudoatom(object="PS_donor_hotspot_0", pos=(-10.0, 12.5, 0.0), color=(1, 1, 1), label=22.0)

cmd.pseudoatom(object="PS_donor_hotspot_1", pos=(-5.5, 7.0, 4.0), color=(1, 1, 1), label=18.9)

cmd.pseudoatom(object="PS_donor_hotspot_2", pos=(-6.0, 10.5, 0.5), color=(1, 1, 1), label=18.4)

cmd.group("label_donor_hotspot", members="PS_donor_hotspot_0")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_0")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_1")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_1")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_2")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_2")
cmd.pseudoatom(object="PS_apolar_hotspot_0", pos=(-9.5, 11.5, 2.5), color=(1, 1, 1), label=19.6)

cmd.pseudoatom(object="PS_apolar_hotspot_1", pos=(-8.25, 10.0, 2.25), color=(1, 1, 1), label=16.7)

cmd.pseudoatom(object="PS_apolar_hotspot_2", pos=(-7.5, 10.25, 0.25), color=(1, 1, 1), label=16.6)

cmd.pseudoatom(object="PS_apolar_hotspot_3", pos=(-6.75, 8.5, 1.5), color=(1, 1, 1), label=17.4)

cmd.pseudoatom(object="PS_apolar_hotspot_4", pos=(9.5, -14.0, 11.5), color=(1, 1, 1), label=17.2)

cmd.pseudoatom(object="PS_apolar_hotspot_5", pos=(11.0, -14.5, 9.5), color=(1, 1, 1), label=16.5)

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
cmd.group("labels_hotspot", members="label_acceptor_hotspot")
cmd.group("labels_hotspot", members="label_donor_hotspot")
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


surface_list = {'hotspot': {'fhm': ['surface_acceptor_hotspot', 'surface_donor_hotspot', 'surface_apolar_hotspot'], 'buriedness': ['surface_buriedness_hotspot']}}
surface_max_list = {'hotspot': {'fhm': 22.2, 'buriedness': 8}}

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