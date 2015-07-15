#!/usr/bin/env python 

import ROOT
from ROOT import TChain, ConverterPHYS14
from os.path import split

file_dir = "/gpfs/csic_projects/tier3data/TreesPHYS14/PU20bx25/"

rootfiles = [ "Tree_T2tt_425LSP325.root" ]

o_dir = "./"              
outfiles = []


for f in rootfiles:
    print "Converting file:" + f
    converter = ConverterPHYS14()
    outfiles.append(o_dir+split(f)[1]) 
    tchain = TChain("demo/Tree")
    tchain.Add(file_dir + f)
    tchain.Process(converter,"ofile="+outfiles[-1])


    
    
