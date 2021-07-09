
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
cmd.pseudoatom(object="PS_donor_hotspot_0", pos=(-6.5, 12.5, 41.5), color=(1, 1, 1), label=24.2)

cmd.pseudoatom(object="PS_donor_hotspot_1", pos=(-1.0, 15.0, 45.0), color=(1, 1, 1), label=24.1)

cmd.pseudoatom(object="PS_donor_hotspot_2", pos=(2.5, 14.5, 39.5), color=(1, 1, 1), label=23.8)

cmd.pseudoatom(object="PS_donor_hotspot_3", pos=(0.5, 16.5, 44.0), color=(1, 1, 1), label=22.5)

cmd.pseudoatom(object="PS_donor_hotspot_4", pos=(0.0, 22.0, 41.5), color=(1, 1, 1), label=22.1)

cmd.pseudoatom(object="PS_donor_hotspot_5", pos=(-3.0, 11.5, 44.0), color=(1, 1, 1), label=21.6)

cmd.pseudoatom(object="PS_donor_hotspot_6", pos=(-0.5, 13.5, 37.5), color=(1, 1, 1), label=19.7)

cmd.pseudoatom(object="PS_donor_hotspot_7", pos=(-4.5, 8.5, 34.0), color=(1, 1, 1), label=19.5)

cmd.pseudoatom(object="PS_donor_hotspot_8", pos=(1.5, 12.5, 41.5), color=(1, 1, 1), label=19.4)

cmd.pseudoatom(object="PS_donor_hotspot_9", pos=(-8.5, 14.0, 38.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_donor_hotspot_10", pos=(-2.0, 15.0, 40.0), color=(1, 1, 1), label=18.9)

cmd.pseudoatom(object="PS_donor_hotspot_11", pos=(0.0, 12.0, 46.0), color=(1, 1, 1), label=18.8)

cmd.pseudoatom(object="PS_donor_hotspot_12", pos=(-0.5, 17.5, 48.5), color=(1, 1, 1), label=18.3)

cmd.pseudoatom(object="PS_donor_hotspot_13", pos=(9.0, 19.5, 30.5), color=(1, 1, 1), label=17.8)

cmd.pseudoatom(object="PS_donor_hotspot_14", pos=(1.0, 12.0, 44.5), color=(1, 1, 1), label=17.7)

cmd.pseudoatom(object="PS_donor_hotspot_15", pos=(-8.5, 6.5, 34.0), color=(1, 1, 1), label=16.4)

cmd.pseudoatom(object="PS_donor_hotspot_16", pos=(5.5, 15.0, 36.5), color=(1, 1, 1), label=16.3)

cmd.pseudoatom(object="PS_donor_hotspot_17", pos=(9.0, 18.0, 33.0), color=(1, 1, 1), label=16.2)

cmd.pseudoatom(object="PS_donor_hotspot_18", pos=(-13.5, 6.5, 36.5), color=(1, 1, 1), label=16.2)

cmd.pseudoatom(object="PS_donor_hotspot_19", pos=(7.5, 14.5, 34.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_donor_hotspot_20", pos=(2.0, 10.5, 30.0), color=(1, 1, 1), label=15.9)

cmd.pseudoatom(object="PS_donor_hotspot_21", pos=(5.0, 13.5, 32.0), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_donor_hotspot_22", pos=(-8.0, 8.5, 37.5), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_donor_hotspot_23", pos=(0.5, 10.0, 39.0), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_donor_hotspot_24", pos=(1.0, 10.5, 37.0), color=(1, 1, 1), label=15.0)

cmd.pseudoatom(object="PS_donor_hotspot_25", pos=(-11.5, 8.5, 34.5), color=(1, 1, 1), label=14.0)

cmd.pseudoatom(object="PS_donor_hotspot_26", pos=(-2.5, 9.5, 41.5), color=(1, 1, 1), label=13.7)

cmd.pseudoatom(object="PS_donor_hotspot_27", pos=(-6.5, 7.0, 37.5), color=(1, 1, 1), label=9.4)

cmd.pseudoatom(object="PS_donor_hotspot_28", pos=(-6.5, 23.5, 47.5), color=(1, 1, 1), label=7.9)

cmd.pseudoatom(object="PS_donor_hotspot_29", pos=(-6.5, 22.5, 49.0), color=(1, 1, 1), label=5.4)

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
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_13")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_13")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_14")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_14")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_15")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_15")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_16")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_16")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_17")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_17")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_18")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_18")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_19")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_19")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_20")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_20")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_21")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_21")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_22")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_22")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_23")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_23")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_24")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_24")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_25")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_25")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_26")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_26")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_27")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_27")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_28")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_28")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_29")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_29")
cmd.pseudoatom(object="PS_acceptor_hotspot_0", pos=(-8.0, 14.0, 41.5), color=(1, 1, 1), label=23.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_1", pos=(-2.0, 20.0, 44.0), color=(1, 1, 1), label=21.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_2", pos=(-10.0, 15.0, 35.5), color=(1, 1, 1), label=20.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_3", pos=(-11.0, 12.5, 40.0), color=(1, 1, 1), label=20.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_4", pos=(-16.0, 15.5, 38.0), color=(1, 1, 1), label=20.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_5", pos=(-19.5, 25.5, 43.0), color=(1, 1, 1), label=20.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_6", pos=(-1.0, 14.5, 46.0), color=(1, 1, 1), label=20.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_7", pos=(3.5, 19.0, 37.5), color=(1, 1, 1), label=20.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_8", pos=(-12.0, 12.0, 37.5), color=(1, 1, 1), label=20.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_9", pos=(2.5, 14.5, 39.0), color=(1, 1, 1), label=19.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_10", pos=(0.0, 19.0, 43.5), color=(1, 1, 1), label=19.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_11", pos=(1.5, 16.5, 41.0), color=(1, 1, 1), label=19.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_12", pos=(0.0, 16.0, 44.5), color=(1, 1, 1), label=19.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_13", pos=(-15.5, 15.0, 35.5), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_14", pos=(-2.0, 13.5, 43.0), color=(1, 1, 1), label=18.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_15", pos=(-6.5, 7.5, 34.0), color=(1, 1, 1), label=17.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_16", pos=(2.0, 11.0, 30.0), color=(1, 1, 1), label=17.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_17", pos=(-2.0, 16.0, 47.5), color=(1, 1, 1), label=17.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_18", pos=(-17.5, 7.5, 36.5), color=(1, 1, 1), label=17.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_19", pos=(-3.0, 17.5, 45.5), color=(1, 1, 1), label=17.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_20", pos=(-5.5, 13.5, 42.0), color=(1, 1, 1), label=16.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_21", pos=(8.0, 14.5, 31.5), color=(1, 1, 1), label=16.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_22", pos=(5.5, 12.0, 30.5), color=(1, 1, 1), label=16.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_23", pos=(-2.0, 14.0, 48.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_24", pos=(-0.5, 18.0, 48.0), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_25", pos=(9.0, 14.0, 30.0), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_26", pos=(-11.5, 11.5, 34.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_27", pos=(-5.5, 12.0, 39.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_28", pos=(-9.0, 7.5, 34.0), color=(1, 1, 1), label=15.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_29", pos=(10.5, 19.5, 31.5), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_30", pos=(2.5, 8.0, 30.0), color=(1, 1, 1), label=15.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_31", pos=(7.5, 14.5, 34.0), color=(1, 1, 1), label=15.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_32", pos=(-5.0, 8.5, 33.0), color=(1, 1, 1), label=15.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_33", pos=(-1.0, 13.0, 40.0), color=(1, 1, 1), label=15.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_34", pos=(7.0, 11.0, 32.0), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_35", pos=(-13.5, 12.0, 40.0), color=(1, 1, 1), label=15.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_36", pos=(1.0, 11.5, 41.5), color=(1, 1, 1), label=15.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_37", pos=(4.0, 19.5, 40.0), color=(1, 1, 1), label=14.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_38", pos=(-12.0, 8.5, 34.5), color=(1, 1, 1), label=14.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_39", pos=(0.5, 13.0, 43.5), color=(1, 1, 1), label=14.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_40", pos=(-14.5, 6.5, 36.5), color=(1, 1, 1), label=14.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_41", pos=(1.5, 11.0, 35.5), color=(1, 1, 1), label=14.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_42", pos=(11.0, 18.0, 32.5), color=(1, 1, 1), label=13.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_43", pos=(-3.5, 7.5, 31.5), color=(1, 1, 1), label=13.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_44", pos=(-4.0, 10.0, 39.0), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_45", pos=(-12.0, 8.5, 36.5), color=(1, 1, 1), label=13.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_46", pos=(-5.0, 10.5, 41.5), color=(1, 1, 1), label=13.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_47", pos=(-11.5, 7.0, 38.5), color=(1, 1, 1), label=13.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_48", pos=(6.5, 13.5, 36.0), color=(1, 1, 1), label=13.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_49", pos=(7.0, 11.5, 35.5), color=(1, 1, 1), label=12.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_50", pos=(-4.5, 6.5, 36.5), color=(1, 1, 1), label=12.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_51", pos=(-4.0, 9.5, 34.5), color=(1, 1, 1), label=12.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_52", pos=(1.5, 14.0, 33.0), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_53", pos=(-6.0, 10.5, 44.0), color=(1, 1, 1), label=11.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_54", pos=(-6.0, 23.5, 45.0), color=(1, 1, 1), label=11.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_55", pos=(4.0, 7.5, 31.5), color=(1, 1, 1), label=10.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_56", pos=(7.0, 13.5, 39.5), color=(1, 1, 1), label=9.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_57", pos=(5.5, 7.0, 30.0), color=(1, 1, 1), label=9.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_58", pos=(4.5, 16.0, 31.0), color=(1, 1, 1), label=8.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_59", pos=(-6.5, 7.0, 36.5), color=(1, 1, 1), label=8.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_60", pos=(3.5, 9.5, 40.0), color=(1, 1, 1), label=8.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_61", pos=(8.0, 18.5, 34.0), color=(1, 1, 1), label=6.6)

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
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_9")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_9")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_10")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_10")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_11")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_11")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_12")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_12")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_13")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_13")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_14")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_14")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_15")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_15")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_16")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_16")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_17")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_17")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_18")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_18")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_19")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_19")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_20")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_20")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_21")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_21")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_22")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_22")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_23")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_23")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_24")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_24")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_25")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_25")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_26")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_26")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_27")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_27")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_28")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_28")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_29")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_29")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_30")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_30")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_31")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_31")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_32")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_32")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_33")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_33")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_34")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_34")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_35")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_35")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_36")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_36")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_37")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_37")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_38")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_38")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_39")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_39")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_40")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_40")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_41")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_41")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_42")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_42")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_43")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_43")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_44")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_44")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_45")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_45")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_46")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_46")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_47")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_47")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_48")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_48")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_49")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_49")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_50")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_50")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_51")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_51")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_52")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_52")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_53")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_53")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_54")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_54")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_55")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_55")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_56")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_56")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_57")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_57")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_58")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_58")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_59")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_59")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_60")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_60")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_61")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_61")
cmd.pseudoatom(object="PS_apolar_hotspot_0", pos=(-12.75, 14.5, 36.5), color=(1, 1, 1), label=18.2)

cmd.pseudoatom(object="PS_apolar_hotspot_1", pos=(-11.0, 13.5, 36.0), color=(1, 1, 1), label=19.4)

cmd.pseudoatom(object="PS_apolar_hotspot_2", pos=(1.5, 18.0, 41.0), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_3", pos=(1.0, 20.5, 41.0), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_apolar_hotspot_4", pos=(-11.0, 10.0, 38.0), color=(1, 1, 1), label=18.3)

cmd.pseudoatom(object="PS_apolar_hotspot_5", pos=(-19.5, 26.0, 43.0), color=(1, 1, 1), label=17.0)

cmd.pseudoatom(object="PS_apolar_hotspot_6", pos=(-8.0, 14.0, 42.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_apolar_hotspot_7", pos=(9.0, 16.5, 30.0), color=(1, 1, 1), label=14.8)

cmd.pseudoatom(object="PS_apolar_hotspot_8", pos=(7.5, 12.5, 30.0), color=(1, 1, 1), label=14.7)

cmd.pseudoatom(object="PS_apolar_hotspot_9", pos=(-6.0, 9.0, 39.0), color=(1, 1, 1), label=14.2)

cmd.pseudoatom(object="PS_apolar_hotspot_10", pos=(5.0, 11.5, 31.5), color=(1, 1, 1), label=13.9)

cmd.pseudoatom(object="PS_apolar_hotspot_11", pos=(-8.5, 7.0, 33.5), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_apolar_hotspot_12", pos=(-5.5, 9.0, 43.0), color=(1, 1, 1), label=13.0)

cmd.pseudoatom(object="PS_apolar_hotspot_13", pos=(-6.0, 7.5, 33.5), color=(1, 1, 1), label=12.6)

cmd.pseudoatom(object="PS_apolar_hotspot_14", pos=(-14.5, 6.5, 37.5), color=(1, 1, 1), label=12.4)

cmd.pseudoatom(object="PS_apolar_hotspot_15", pos=(3.5, 7.5, 30.0), color=(1, 1, 1), label=11.9)

cmd.pseudoatom(object="PS_apolar_hotspot_16", pos=(7.5, 8.0, 31.0), color=(1, 1, 1), label=11.4)

cmd.pseudoatom(object="PS_apolar_hotspot_17", pos=(9.0, 20.5, 30.0), color=(1, 1, 1), label=11.0)

cmd.pseudoatom(object="PS_apolar_hotspot_18", pos=(-6.5, 25.5, 44.0), color=(1, 1, 1), label=10.9)

cmd.pseudoatom(object="PS_apolar_hotspot_19", pos=(0.0, 8.0, 41.5), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_apolar_hotspot_20", pos=(-7.5, 24.0, 45.0), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_apolar_hotspot_21", pos=(-16.0, 26.0, 30.0), color=(1, 1, 1), label=8.6)

cmd.pseudoatom(object="PS_apolar_hotspot_22", pos=(9.5, 17.5, 35.5), color=(1, 1, 1), label=6.0)

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
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_17")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_17")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_18")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_18")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_19")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_19")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_20")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_20")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_21")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_21")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_22")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_22")
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
surface_max_list = {'hotspot': {'fhm': 24.2}}

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