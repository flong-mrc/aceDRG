# Python script
#
#
## The date of last modification: 06/11/2017
#

from __future__ import print_function
import os,os.path,sys
import glob,shutil
import re,string

class PeriodicTab(dict):

    def __init__(self):

        self["H"]            = {}
        self["H"]["row"]     = 1
        self["H"]["group"]   = 1
        self["H"]["matType"] = 1
        self["H"]["atomNum"] = 1
        self["H"]["val"]     = 1
        
        self["D"]            = {}
        self["D"]["row"]     = 1
        self["D"]["group"]   = 1
        self["D"]["matType"] = 1
        self["D"]["atomNum"] = 1
        self["D"]["val"]     = 1
        
        # Organic set or Non-metal
  
        self["C"]            = {}
        self["C"]["row"]     = 2
        self["C"]["group"]   = 14
        self["C"]["matType"] = 2
        self["C"]["atomNum"] = 6 
        self["C"]["val"]     = 4 
        
        self["N"]              = {}
        self["N"]["row"]       = 2 
        self["N"]["group"]     = 15 
        self["N"]["matType"]   = 2 
        self["N"]["atomNum"]   = 7 
        self["N"]["val"]       = 3 
        self["N"]["extraVal"]  = []
        self["N"]["extraVal"].append(5)
      
         
        self["O"]              = {}
        self["O"]["row"]       = 2 
        self["O"]["group"]     = 16 
        self["O"]["matType"]   = 2 
        self["O"]["atomNum"]   = 8 
        self["O"]["val"]       = 2 
        
        self["P"]              = {}
        self["P"]["row"]       = 3 
        self["P"]["group"]     = 15 
        self["P"]["matType"]   = 2 
        self["P"]["atomNum"]   = 15 
        self["P"]["val"]       = 5 
        self["P"]["extraVal"]  = []
        self["P"]["extraVal"].append(3)
        
        self["S"]              = {}
        self["S"]["row"]       = 3 
        self["S"]["group"]     = 16 
        self["S"]["matType"]   = 2 
        self["S"]["atomNum"]   = 16 
        self["S"]["val"]       = 6 
        self["S"]["extraVal"]  = []
        self["S"]["extraVal"].append(2)
        self["S"]["extraVal"].append(4)
        
        self["Se"]             = {}
        self["Se"]["row"]      = 4 
        self["Se"]["group"]    = 16 
        self["Se"]["matType"]  = 2 
        self["Se"]["atomNum"]  = 34 
        self["Se"]["val"]      = 6 
        self["Se"]["extraVal"] = []
        self["Se"]["extraVal"].append(2)
        self["Se"]["extraVal"].append(4)
        
        self["B"]               = {}
        self["B"]["row"]        = 2 
        self["B"]["group"]      = 13 
        self["B"]["matType"]    = 7 
        self["B"]["atomNum"]    = 5 
        self["B"]["val"]        = 3 
        
        # Halogens
        self["F"]            = {}
        self["F"]["row"]     = 2 
        self["F"]["group"]   = 17 
        self["F"]["matType"] = 8 
        self["F"]["atomNum"] = 9 
        self["F"]["val"]     = 1 
        
        self["Cl"]              = {}
        self["Cl"]["row"]     = 3 
        self["Cl"]["group"]   = 17 
        self["Cl"]["matType"] = 8 
        self["Cl"]["atomNum"] = 17 
        self["Cl"]["val"]     = 1 
        
        self["Br"]            = {}
        self["Br"]["row"]     = 4 
        self["Br"]["group"]   = 17 
        self["Br"]["matType"] = 8 
        self["Br"]["atomNum"] = 35 
        self["Br"]["val"]     = 1 
        
        self["I"]             = {}
        self["I"]["row"]      = 5 
        self["I"]["group"]    = 17 
        self["I"]["matType"]  = 8 
        self["I"]["atomNum"]  = 53 
        self["I"]["val"]      = 1 
        
        self["At"]             = {}
        self["At"]["row"]      = 6 
        self["At"]["group"]    = 17 
        self["At"]["matType"]  = 8 
        self["At"]["atomNum"]  = 85 
        self["At"]["val"]      = 1  
        
        # Alkali-Metals
        self["Li"]             = {}
        self["Li"]["row"]      = 2 
        self["Li"]["group"]    = 1 
        self["Li"]["matType"]  = 3 
        self["Li"]["atomNum"]  = 3 
        self["Li"]["val"]      = 1 
        
        self["Na"]             = {}
        self["Na"]["row"]      = 3 
        self["Na"]["group"]    = 1 
        self["Na"]["matType"]  = 3 
        self["Na"]["atomNum"]  = 11 
        self["Na"]["val"]      = 1 
       
        self["K"]              = {} 
        self["K"]["row"]       = 4 
        self["K"]["group"]     = 1 
        self["K"]["matType"]   = 3 
        self["K"]["atomNum"]   = 19 
        self["K"]["val"]       = 1 
        
        self["Rb"]             = {}
        self["Rb"]["row"]      = 5 
        self["Rb"]["group"]    = 1 
        self["Rb"]["matType"]  = 3 
        self["Rb"]["atomNum"]  = 37 
        self["Rb"]["val"]      = 1 
        
        self["Cs"]             = {}
        self["Cs"]["row"]      = 6 
        self["Cs"]["group"]    = 1 
        self["Cs"]["matType"]  = 3 
        self["Cs"]["atomNum"]  = 55 
        self["Cs"]["val"]      = 1 
        self["Cs"]["vdw"]      = 2.00 
                
        self["Fr"]             = {}
        self["Fr"]["row"]      = 7 
        self["Fr"]["group"]    = 1 
        self["Fr"]["matType"]  = 3 
        self["Fr"]["atomNum"]  = 87 
        self["Fr"]["val"]      = 1 
        
        # Alkaline-earth-Metals
        self["Be"]             = {}
        self["Be"]["row"]      = 2 
        self["Be"]["group"]    = 2 
        self["Be"]["matType"]  = 4 
        self["Be"]["atomNum"]  = 4 
        self["Be"]["val"]      = 2 
        
        self["Mg"]             = {}
        self["Mg"]["row"]      = 3 
        self["Mg"]["group"]    = 2 
        self["Mg"]["matType"]  = 4 
        self["Mg"]["atomNum"]  = 12 
        self["Mg"]["val"]      = 2 
        
        self["Ca"]             = {} 
        self["Ca"]["row"]      = 4 
        self["Ca"]["group"]    = 2 
        self["Ca"]["matType"]  = 4 
        self["Ca"]["atomNum"]  = 20 
        self["Ca"]["val"]      = 2 
        
        self["Sr"]             = {}
        self["Sr"]["row"]      = 5 
        self["Sr"]["group"]    = 2 
        self["Sr"]["matType"]  = 4 
        self["Sr"]["atomNum"]  = 38 
        self["Sr"]["val"]      = 2 
        
        self["Ba"]             = {}
        self["Ba"]["row"]      = 6 
        self["Ba"]["group"]    = 2 
        self["Ba"]["matType"]  = 4 
        self["Ba"]["atomNum"]  = 56 
        self["Ba"]["val"]      = 2 
        
        self["Ra"]             = {}
        self["Ra"]["row"]      = 7 
        self["Ra"]["group"]    = 2 
        self["Ra"]["matType"]  = 4 
        self["Ra"]["atomNum"]  = 88 
        self["Ra"]["val"]      = 2 
        
        
        # Transition-Metals

        self["Sc"]             = {}
        self["Sc"]["row"]      = 4 
        self["Sc"]["group"]    = 3 
        self["Sc"]["matType"]  = 5 
        self["Sc"]["atomNum"]  = 21 
        self["Sc"]["val"]      = 3 
        
        self["Y"]              = {}
        self["Y"]["row"]       = 5 
        self["Y"]["group"]     = 3 
        self["Y"]["matType"]   = 5 
        self["Y"]["atomNum"]   = 39  
        self["Y"]["val"]       = 3 
        
        self["Ti"]             = {}
        self["Ti"]["row"]      = 4 
        self["Ti"]["group"]    = 4 
        self["Ti"]["matType"]  = 5 
        self["Ti"]["atomNum"]  = 22 
        self["Ti"]["val"]      = 4 
        
        self["Zr"]             = {}
        self["Zr"]["row"]      = 5 
        self["Zr"]["group"]    = 4 
        self["Zr"]["matType"]  = 5 
        self["Zr"]["atomNum"]  = 40 
        self["Zr"]["val"]      = 4 
        
        self["Hf"]             = {}
        self["Hf"]["row"]      = 6 
        self["Hf"]["group"]    = 4 
        self["Hf"]["matType"]  = 5 
        self["Hf"]["atomNum"]  = 72 
        self["Hf"]["val"]      = 4 
        
        self["Rf"]             = {}
        self["Rf"]["row"]      = 7 
        self["Rf"]["group"]    = 4 
        self["Rf"]["matType"]  = 5 
        self["Rf"]["atomNum"]  = 104 
        self["Rf"]["val"]      = 4 
        
        self["V"]              = {}
        self["V"]["row"]       = 4 
        self["V"]["group"]     = 5 
        self["V"]["matType"]   = 5 
        self["V"]["atomNum"]   = 23 
        self["V"]["val"]       = 5 
        
        self["Nb"]             = {}
        self["Nb"]["row"]      = 5 
        self["Nb"]["group"]    = 5 
        self["Nb"]["matType"]  = 5 
        self["Nb"]["atomNum"]  = 41 
        self["Nb"]["val"]      = 5 
        
        self["Ta"]              = {}
        self["Ta"]["row"]       = 6 
        self["Ta"]["group"]     = 5 
        self["Ta"]["matType"]   = 5 
        self["Ta"]["atomNum"]   = 73 
        self["Ta"]["val"]       = 5 
        
        self["Db"]              = {}
        self["Db"]["row"]       = 7 
        self["Db"]["group"]     = 5 
        self["Db"]["matType"]   = 5 
        self["Db"]["atomNum"]   = 105 
        self["Db"]["val"]       = 5 
        
        self["Cr"]              = {}
        self["Cr"]["row"]       = 4 
        self["Cr"]["group"]     = 6 
        self["Cr"]["matType"]   = 5 
        self["Cr"]["atomNum"]   = 24 
        self["Cr"]["val"]       = 6 
        
        self["Mo"]              = {}
        self["Mo"]["row"]       = 5 
        self["Mo"]["group"]     = 6 
        self["Mo"]["matType"]   = 5 
        self["Mo"]["atomNum"]   = 42 
        self["Mo"]["val"]       = 6 
        
        self["W"]               = {}
        self["W"]["row"]        = 6 
        self["W"]["group"]      = 6 
        self["W"]["matType"]    = 5 
        self["W"]["atomNum"]    = 74 
        self["W"]["val"]        = 6 
        
        self["Sg"]              = {}
        self["Sg"]["row"]       = 7 
        self["Sg"]["group"]     = 6 
        self["Sg"]["matType"]   = 5 
        self["Sg"]["atomNum"]   = 106 
        self["Sg"]["val"]       = 6 

        self["Mn"]              = {}
        self["Mn"]["row"]       = 4 
        self["Mn"]["group"]     = 7 
        self["Mn"]["matType"]   = 5 
        self["Mn"]["atomNum"]   = 25 
        self["Mn"]["val"]       = 7 
        
        self["Tc"]              = {}
        self["Tc"]["row"]       = 5 
        self["Tc"]["group"]     = 7 
        self["Tc"]["matType"]   = 5 
        self["Tc"]["atomNum"]   = 43 
        self["Tc"]["val"]       = 7 
        
        self["Re"]              = {}
        self["Re"]["row"]       = 6 
        self["Re"]["group"]     = 7 
        self["Re"]["matType"]   = 5 
        self["Re"]["atomNum"]   = 75 
        self["Re"]["val"]       = 7 
        
        self["Bh"]              = {}
        self["Bh"]["row"]       = 7 
        self["Bh"]["group"]     = 7 
        self["Bh"]["matType"]   = 5 
        self["Bh"]["atomNum"]   = 107 
        self["Bh"]["val"]       = 7 
        
        self["Fe"]              = {}
        self["Fe"]["row"]       = 4 
        self["Fe"]["group"]     = 8 
        self["Fe"]["matType"]   = 5 
        self["Fe"]["atomNum"]   = 26 
        self["Fe"]["val"]       = 6 
      
        self["Ru"]              = {}
        self["Ru"]["row"]       = 5 
        self["Ru"]["group"]     = 8 
        self["Ru"]["matType"]   = 5 
        self["Ru"]["atomNum"]   = 44 
        self["Ru"]["val"]       = 8 
        
        self["Os"]              = {}
        self["Os"]["row"]       = 6 
        self["Os"]["group"]     = 8 
        self["Os"]["matType"]   = 5 
        self["Os"]["atomNum"]   = 76 
        self["Os"]["val"]       = 8 
        
        self["Hs"]              = {}
        self["Hs"]["row"]       = 7 
        self["Hs"]["group"]     = 8 
        self["Hs"]["matType"]   = 5 
        self["Hs"]["atomNum"]   = 108 
        self["Hs"]["val"]       = 8 
        
        self["Co"]              = {}
        self["Co"]["row"]       = 4 
        self["Co"]["group"]     = 9 
        self["Co"]["matType"]   = 5 
        self["Co"]["atomNum"]   = 27 
        self["Co"]["val"]       = 5 
       
        self["Rh"]              = {} 
        self["Rh"]["row"]       = 5 
        self["Rh"]["group"]     = 9 
        self["Rh"]["matType"]   = 5 
        self["Rh"]["atomNum"]   = 45 
        
        self["Ir"]              = {}
        self["Ir"]["row"]       = 6 
        self["Ir"]["group"]     = 9 
        self["Ir"]["matType"]   = 5 
        self["Ir"]["atomNum"]   = 77 
        self["Ir"]["val"]       = 8 
        
        self["Mt"]              = {}
        self["Mt"]["row"]       = 7 
        self["Mt"]["group"]     = 9 
        self["Mt"]["matType"]   = 5 
        self["Mt"]["atomNum"]   = 109 
        self["Mt"]["val"]       = -1 
        
        self["Ni"]              = {}
        self["Ni"]["row"]       = 4 
        self["Ni"]["group"]     = 10 
        self["Ni"]["matType"]   = 5 
        self["Ni"]["atomNum"]   = 28 
        self["Ni"]["val"]       = 4 
        
        self["Pd"]              = {}
        self["Pd"]["row"]       = 5 
        self["Pd"]["group"]     = 10 
        self["Pd"]["matType"]   = 5 
        self["Pd"]["atomNum"]   = 46 
        self["Pd"]["val"]       = 4 
        
        self["Pt"]              = {}
        self["Pt"]["row"]       = 6 
        self["Pt"]["group"]     = 10 
        self["Pt"]["matType"]   = 5 
        self["Pt"]["atomNum"]   = 78 
        self["Pt"]["val"]       = 6 
        
        self["Ds"]              = {}
        self["Ds"]["row"]       = 7 
        self["Ds"]["group"]     = 10 
        self["Ds"]["matType"]   = 5 
        self["Ds"]["atomNum"]   = 110 
        self["Ds"]["val"]       = -1 
        
        self["Cu"]              = {}
        self["Cu"]["row"]       = 4 
        self["Cu"]["group"]     = 11 
        self["Cu"]["matType"]   = 5 
        self["Cu"]["atomNum"]   = 29 
        self["Cu"]["val"]       = 2 
        
        self["Ag"]              = {}
        self["Ag"]["row"]       = 5 
        self["Ag"]["group"]     = 11 
        self["Ag"]["matType"]   = 5 
        self["Ag"]["atomNum"]   = 47 
        self["Ag"]["val"]       = 4 
        
        self["Au"]              = {}
        self["Au"]["row"]       = 6 
        self["Au"]["group"]     = 11 
        self["Au"]["matType"]   = 5 
        self["Au"]["atomNum"]   = 79 
        self["Au"]["val"]       = 5 
        
        self["Rg"]              = {}
        self["Rg"]["row"]       = 7 
        self["Rg"]["group"]     = 11 
        self["Rg"]["matType"]   = 5 
        self["Rg"]["atomNum"]   = 111 
        self["Rg"]["val"]       = -1 
        
        self["Zn"]              = {}
        self["Zn"]["row"]       = 4 
        self["Zn"]["group"]     = 12 
        self["Zn"]["matType"]   = 5 
        self["Zn"]["atomNum"]   = 30 
        self["Zn"]["val"]       = 2 
        
        self["Cd"]              = {}
        self["Cd"]["row"]       = 5 
        self["Cd"]["group"]     = 12 
        self["Cd"]["matType"]   = 5 
        self["Cd"]["atomNum"]   = 48 
        self["Cd"]["val"]       = 2 
        
        self["Hg"]              = {}
        self["Hg"]["row"]       = 6 
        self["Hg"]["group"]     = 12 
        self["Hg"]["matType"]   = 5 
        self["Hg"]["atomNum"]   = 80 
        self["Hg"]["val"]       = 2 
        
        self["Cn"]              = {}
        self["Cn"]["row"]       = 7 
        self["Cn"]["group"]     = 12 
        self["Cn"]["matType"]   = 5 
        self["Cn"]["atomNum"]   = 112 
        self["Cn"]["val"]       = -1 
        
        # Other-Metal
        self["Al"]              = {}
        self["Al"]["row"]       = 3 
        self["Al"]["group"]     = 13 
        self["Al"]["matType"]   = 6 
        self["Al"]["atomNum"]   = 13 
        self["Al"]["val"]       = 3 
        
        self["Ga"]              = {}
        self["Ga"]["row"]       = 4 
        self["Ga"]["group"]     = 13 
        self["Ga"]["matType"]   = 6 
        self["Ga"]["atomNum"]   = 31 
        self["Ga"]["val"]       = 3 
        
        self["In"]              = {}
        self["In"]["row"]       = 5 
        self["In"]["group"]     = 13 
        self["In"]["matType"]   = 6 
        self["In"]["atomNum"]   = 49 
        self["In"]["val"]       = 3 
        
        self["Ti"]              = {}
        self["Ti"]["row"]       = 6 
        self["Ti"]["group"]     = 13 
        self["Ti"]["matType"]   = 6 
        self["Ti"]["atomNum"]   = 81 
        self["Ti"]["val"]       = 3 
        
        self["Sn"]              = {}
        self["Sn"]["row"]       = 5 
        self["Sn"]["group"]     = 14 
        self["Sn"]["matType"]   = 6 
        self["Sn"]["atomNum"]   = 50 
        self["Sn"]["val"]       = 4 
        
        self["Pb"]              = {}
        self["Pb"]["row"]       = 6 
        self["Pb"]["group"]     = 14 
        self["Pb"]["matType"]   = 6 
        self["Pb"]["atomNum"]   = 82 
        self["Pb"]["val"]       = 4 
        
        self["Bi"]              = {}
        self["Bi"]["row"]       = 6 
        self["Bi"]["group"]     = 15 
        self["Bi"]["matType"]   = 6 
        self["Bi"]["atomNum"]   = 83 
        self["Bi"]["val"]       = 5 
        

        # Semimetallics
        
        self["Si"]              = {}
        self["Si"]["row"]       = 3 
        self["Si"]["group"]     = 14 
        self["Si"]["matType"]   = 7 
        self["Si"]["atomNum"]   = 14 
        self["Si"]["val"]       = 4 
        
        self["Ge"]              = {}
        self["Ge"]["row"]       = 4 
        self["Ge"]["group"]     = 14 
        self["Ge"]["matType"]   = 7 
        self["Ge"]["atomNum"]   = 32 
        self["Ge"]["val"]       = 4 
        
        self["As"]              = {}
        self["As"]["row"]       = 4 
        self["As"]["group"]     = 15 
        self["As"]["matType"]   = 7 
        self["As"]["atomNum"]   = 33 
        self["As"]["val"]       = 5 
        
        self["Sb"]              = {}
        self["Sb"]["row"]       = 5 
        self["Sb"]["group"]     = 15 
        self["Sb"]["matType"]   = 7 
        self["Sb"]["atomNum"]   = 51 
        self["Sb"]["val"]       = 5 
        
        self["Te"]              = {}
        self["Te"]["row"]       = 5 
        self["Te"]["group"]     = 16 
        self["Te"]["matType"]   = 7 
        self["Te"]["atomNum"]   = 52 
        self["Te"]["val"]       = 6 
        
        self["Po"]              = {}
        self["Po"]["row"]       = 6 
        self["Po"]["group"]     = 16 
        self["Po"]["matType"]   = 7 
        self["Po"]["atomNum"]   = 84 
        self["Po"]["val"]       = 6 
        
        # Rare earth metals
        
        self["La"]              = {}
        self["La"]["row"]       = 6 
        self["La"]["group"]     = 3 
        self["La"]["matType"]   = 9 
        self["La"]["atomNum"]   = 57 
        self["La"]["val"]       = 3 
        
        self["Ce"]              = {}
        self["Ce"]["row"]       = 6 
        self["Ce"]["group"]     = 3 
        self["Ce"]["matType"]   = 9 
        self["Ce"]["atomNum"]   = 58 
        self["Ce"]["val"]       = 4 
        
        self["Pr"]              = {}
        self["Pr"]["row"]       = 6 
        self["Pr"]["group"]     = 3 
        self["Pr"]["matType"]   = 9 
        self["Pr"]["atomNum"]   = 59 
        self["Pr"]["val"]       = 4 
        
        self["Nd"]              = {}
        self["Nd"]["row"]       = 6 
        self["Nd"]["group"]     = 3 
        self["Nd"]["matType"]   = 9 
        self["Nd"]["atomNum"]   = 60 
        self["Nd"]["val"]       = 3 
        
        self["Pm"]              = {}
        self["Pm"]["row"]       = 6 
        self["Pm"]["group"]     = 3 
        self["Pm"]["matType"]   = 9 
        self["Pm"]["atomNum"]   = 61 
        
        self["Sm"]              = {}
        self["Sm"]["row"]       = 6 
        self["Sm"]["group"]     = 3 
        self["Sm"]["matType"]   = 9 
        self["Sm"]["atomNum"]   = 62 
        self["Sm"]["val"]       = 3 
        
        self["Eu"]              = {}
        self["Eu"]["row"]       = 6 
        self["Eu"]["group"]     = 3 
        self["Eu"]["matType"]   = 9 
        self["Eu"]["atomNum"]   = 63 
        self["Eu"]["val"]       = 6 
        
        self["Gd"]              = {}
        self["Gd"]["row"]       = 6 
        self["Gd"]["group"]     = 3 
        self["Gd"]["matType"]   = 9 
        self["Gd"]["atomNum"]   = 64 
        self["Gd"]["val"]       = 3 
        
        self["Tb"]              = {}
        self["Tb"]["row"]       = 6 
        self["Tb"]["group"]     = 3 
        self["Tb"]["matType"]   = 9 
        self["Tb"]["atomNum"]   = 65 
        self["Tb"]["val"]       = 4 
        
        self["Dy"]              = {} 
        self["Dy"]["row"]       = 6 
        self["Dy"]["group"]     = 3 
        self["Dy"]["matType"]   = 9 
        self["Dy"]["atomNum"]   = 66 
        self["Dy"]["val"]       = 3 
        
        self["Ho"]              = {}
        self["Ho"]["row"]       = 6 
        self["Ho"]["group"]     = 3 
        self["Ho"]["matType"]   = 9 
        self["Ho"]["atomNum"]   = 67 
        self["Ho"]["val"]       = 3 
       
        self["Er"]              = {}
        self["Er"]["row"]       = 6 
        self["Er"]["group"]     = 3 
        self["Er"]["matType"]   = 9 
        self["Er"]["atomNum"]   = 68 
        self["Er"]["val"]       = 3 
        
        self["Tm"]              = {}
        self["Tm"]["row"]       = 6 
        self["Tm"]["group"]     = 3 
        self["Tm"]["matType"]   = 9 
        self["Tm"]["atomNum"]   = 69 
        self["Tm"]["val"]       = 3 
        
        self["Yb"]              = {}
        self["Yb"]["row"]       = 6 
        self["Yb"]["group"]     = 3 
        self["Yb"]["matType"]   = 9 
        self["Yb"]["atomNum"]   = 70 
        self["Yb"]["val"]       = 1 
        
        self["Lu"]              = {}
        self["Lu"]["row"]       = 6 
        self["Lu"]["group"]     = 3 
        self["Lu"]["matType"]   = 9 
        self["Lu"]["atomNum"]   = 71 
        self["Lu"]["val"]       = 1 
        
        self["Ac"]              = {}
        self["Ac"]["row"]       = 7 
        self["Ac"]["group"]     = 3 
        self["Ac"]["matType"]   = 9 
        self["Ac"]["atomNum"]   = 89 
        self["Ac"]["val"]       = 3 
        
        self["Th"]              = {}
        self["Th"]["row"]       = 7 
        self["Th"]["group"]     = 3 
        self["Th"]["matType"]   = 9 
        self["Th"]["atomNum"]   = 90 
        self["Th"]["val"]       = 4 
        
        self["Pa"]              = {}
        self["Pa"]["row"]       = 7 
        self["Pa"]["group"]     = 3 
        self["Pa"]["matType"]   = 9 
        self["Pa"]["atomNum"]   = 91 
        self["Pa"]["val"]       = 4 
        
        self["U"]              = {}
        self["U"]["row"]       = 7 
        self["U"]["group"]     = 3 
        self["U"]["matType"]   = 9 
        self["U"]["atomNum"]   = 92 
        self["U"]["val"]       = 6 
        
        self["Np"]              = {}
        self["Np"]["row"]       = 7 
        self["Np"]["group"]     = 3 
        self["Np"]["matType"]   = 9 
        self["Np"]["atomNum"]   = 93 
        self["Np"]["val"]       = 4 
        
        self["Pu"]              = {}
        self["Pu"]["row"]       = 7 
        self["Pu"]["group"]     = 3 
        self["Pu"]["matType"]   = 9 
        self["Pu"]["atomNum"]   = 94 
        self["Pu"]["val"]       = 6 
        
        self["Am"]              = {}
        self["Am"]["row"]       = 7 
        self["Am"]["group"]     = 3 
        self["Am"]["matType"]   = 9 
        self["Am"]["atomNum"]   = 95 
        self["Am"]["val"]       = 6 
        
        self["Cm"]              = {}
        self["Cm"]["row"]       = 7 
        self["Cm"]["group"]     = 3 
        self["Cm"]["matType"]   = 9 
        self["Cm"]["atomNum"]   = 96 
        self["Cm"]["val"]       = 4 
        
        self["Bk"]              = {}
        self["Bk"]["row"]       = 7 
        self["Bk"]["group"]     = 3 
        self["Bk"]["matType"]   = 9 
        self["Bk"]["atomNum"]   = 97 
        self["Bk"]["val"]       = 4 
        
        self["Cf"]              = {}
        self["Cf"]["row"]       = 7 
        self["Cf"]["group"]     = 3 
        self["Cf"]["matType"]   = 9 
        self["Cf"]["atomNum"]   = 98 
        self["Cf"]["val"]       = 4 
        
        self["Es"]              = {}
        self["Es"]["row"]       = 7 
        self["Es"]["group"]     = 3 
        self["Es"]["matType"]   = 9 
        self["Es"]["atomNum"]   = 99 
        self["Es"]["val"]       = 3 
        
        self["Fm"]              = {}
        self["Fm"]["row"]       = 7 
        self["Fm"]["group"]     = 3 
        self["Fm"]["matType"]   = 9 
        self["Fm"]["atomNum"]   = 100 
        self["Fm"]["val"]       = 3 
        
        self["Md"]              = {}
        self["Md"]["row"]       = 7 
        self["Md"]["group"]     = 3 
        self["Md"]["matType"]   = 9 
        self["Md"]["atomNum"]   = 101 
        self["Md"]["val"]       = 3 
        
        self["No"]              = {}
        self["No"]["row"]       = 7 
        self["No"]["group"]     = 3 
        self["No"]["matType"]   = 9 
        self["No"]["atomNum"]   = 102 
        self["No"]["val"]       = 3 
        
        self["Lr"]              = {}
        self["Lr"]["row"]       = 7 
        self["Lr"]["group"]     = 3 
        self["Lr"]["matType"]   = 9 
        self["Lr"]["atomNum"]   = 103 
        self["Lr"]["val"]       = 3 


    def standlizeSymbol(self, tStr):
   
        aRet = ""
        tStr = tStr.strip()
        if len(tStr)==1:
            aRet = tStr.upper()
        elif len(tStr)==2:
            aRet = tStr[0].upper() + tStr[1].lower()
        # for a string of length other than 1 or 2, return a empty standlized string
        return aRet

    def getElementSymb(self, tStr):

        aEl = tStr.strip()
        if len(tStr) == 1:
            aEl = tStr.upper()
        elif len(tStr) == 2:
            aEl = tStr[0] + tStr[1].lower()
        else :
             print("Wrong element symbol ", tStr)
        # for a string of length other than 1 or 2, return the original string 
        return aEl

                
                 
                
                
                
 
    

