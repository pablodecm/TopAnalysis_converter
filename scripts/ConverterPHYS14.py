#!/usr/bin/env python 

import ROOT
from ROOT import TChain
from ROOT import ConverterPHYS14
from os.path import split

file_dir = "/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/"

stop_files = [ "Tree_T2tt_425LSP325.root",
              "Tree_T2tt_500LSP325.root",
              "Tree_T2tt_650LSP325.root",
              "Tree_T2tt_850LSP100.root"]  

ttbar_files = ["Tree_TTJets_MadSpin_{}.root".format(n) for n in range(8)]

rootfiles = stop_files + ttbar_files

o_dir = "/gpfs/csic_projects/cms/pablodcm/data/mut/PHYS14_dilepton/"              
outfiles = []

for f in rootfiles:
    converter = ConverterPHYS14()
    tchain = TChain("demo/Tree")
    outfiles.append(o_dir+split(f)[1]) 
    print "Add: " + f
    tchain.Add(file_dir + f)
    tchain.Process(converter,"ofile="+outfiles[-1])


    
    
