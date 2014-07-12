#!/usr/bin/env python
# Python script
#
import os, sys
import glob, re, shutil
import string, time
import math


def getMod(aV):
    
    t_sqrt = 0.0
    for a_i in aV:
        t_sqrt +=(a_i*a_i)
    t_sqrt = math.sqrt(t_sqrt)
    
    return t_sqrt
    
def getDotP(aV1, aV2):

    t_dot_p = 0.0

    for i in range(len(aV1)):
        t_dot_p +=(aV1[i]*aV2[i])

    return t_dot_p
    
def CountMolsCifs(in_file_name, out_file_name):
    """counting molecules and cif files in a bonds.txt, angles.txt etc."""
         
    try :
        in_file = open(in_file_name, "r")
    except IOError :
        print "%s can not be opened for reading "%in_file_name
        sys.exit(1)
    else :
        try :
            out_file = open(out_file_name, "a")
        except IOError :
            print "%s can not be opened for writing "%out_file_name
            sys.exit(1)
        else:
            allCifMols = {}
            a_line = in_file.readline()
            i_line = 0
            while len(a_line):
                i_line += 1
                str_grp = a_line.strip().split()
                if len(str_grp) == 7:
                    t_cif = str_grp[5].strip()
                    t_mol = str_grp[6].strip()
                    if not allCifMols.has_key(t_cif):
                        allCifMols[t_cif] = []
                    if not t_mol in allCifMols[t_cif]:
                        allCifMols[t_cif].append(t_mol)
                a_line = in_file.readline()

            in_file.close()

            i_cif = 0
            i_mol = 0

            for a_key in allCifMols.keys():
                i_cif += 1
                i_mol += len(allCifMols[a_key])

            out_file.write("Number of cif file in %s is : %d \n"%(in_file_name, i_cif))
            out_file.write("Number of molecules in %s is : %d \n"%(in_file_name, i_mol))
            out_file.write("Number of bonds in %s is : %d \n"%(in_file_name, i_line))

            out_file.close()


def delMolsInCif(in_file_name, out_kept_name, out_deleted_name, crit_dict, num_cols):
                             
    # Now scan the input file see if items in crit_list exist
    try:
        in_file = open(in_file_name, "r")
    except IOError :
        print "%s can not be opened for reading "%in_file_name
        sys.exit(1)
    else :
        try :
            out_file  = open(out_kept_name, "w")
            out_file2 = open(out_deleted_name, "w")
        except IOError :
            print "%s or %s can not be opened for writing "%(out_file_name, out_deleted_name)
            sys.exit(1)
        else:
            a_line = in_file.readline()
            while len(a_line):
                str_grp = a_line.strip().split()
                if len(str_grp) == num_cols:
                    t_cif = str_grp[7].strip()
                    t_mol = str_grp[8].strip()
                    if crit_dict.has_key(t_cif):
                        if not t_mol in crit_dict[t_cif]:
                            out_file.write(a_line)
                        else:
                            out_file2.write(a_line)
                    else:
                        out_file.write(a_line)
                else:
                    out_file.write(a_line)
                a_line = in_file.readline()

            in_file.close()
            out_file.close()
            out_file2.close()

def allUpper(in_file_name, out_file_name):

    try :
        in_file = open(in_file_name, "r")
    except IOError :
        print "%s can not be opened for reading "%in_file_name
        sys.exit(1)
    else :
        try :
            out_file = open(out_file_name, "w")
        except IOError :
            print "%s can not be opened for writing "%out_file_name
            sys.exit(1)
        else:
            a_line = in_file.readline()
            while len(a_line):
                if len(a_line) !=0 :
                    a_line = a_line.upper()
                    out_file.write(a_line)
                a_line = in_file.readline()

            in_file.close()
            out_file.close()

# comparison functions for sorting
def listComp(a_list, b_list):

    if a_list[-1] > b_list[-1]:
        return 1
    elif a_list[-1] == b_list[-1]:
        return 0
    elif a_list[-1] < b_list[-1]:
        return -1

def listComp2(a_list, b_list):

    if a_list[1] > b_list[1]:
        return 1
    elif a_list[1] < b_list[1]:
        return -1

    if len(a_list[0]) < len(b_list[0]):
        return 1
    elif len(a_list[0]) > len(b_list[0]):
        return -1

    if a_list[0] > b_list[0]:
        return 1
    elif a_list[0] == b_list[0]:
        return 0
    elif a_list[0] < b_list[0]:
        return -1


def dictObjCmp(dict_a, dict_b):
    
    aKey = [len(dict_a["cod_type"]), dict_a["cod_type"]]
    bKey = [len(dict_b["cod_type"]), dict_b["cod_type"]]

    return cmp(aKey, bKey)

def getMax(t_list):
    
    tMax = -1000000000000.0
    if len(t_list):
        tMax =t_list[0]
        for a_v in t_list:
            if a_v > tMax:
                tMax = a_v

    if tMax <= -1000000000000.0:
        print "Max value = ", tMax
        print "Check your data "
    return tMax

def getMin(t_list):
    
    tMin = 1000000000000.0
    if len(t_list):
        tMin =t_list[0]
        for a_v in t_list:
            if a_v < tMin:
                tMin = a_v

    if tMin >= 1000000000000.0:
        print "Min value = ", tMin
        print "Check your data "
    return tMin

def getMean(t_list):
    t_mean = 0.0
    if len(t_list) > 0:
        t_mean = math.fsum(t_list)/len(t_list)
    return t_mean

def getMedian(t_list):

    t_inc = False

    # get median,
    e_len = len(t_list)
    if math.fmod(e_len, 2) != 0.0:
        t_pos =int((e_len-1)/2)
        t_median  = t_list[t_pos]
    else:
        t_pos =int((e_len)/2)
        sum2 = math.fsum([t_list[t_pos-1],
                            t_list[t_pos]])
        t_median = sum2/2.0
        t_inc = True

    return t_median, t_pos, t_inc


def getPopulationVar(t_list, t_mean):

    t_sum     = 0.0
    t_var     = 0.0
    if len(t_list) > 0:
        for a_etr in t_list:
            t_sum +=(a_etr-t_mean)*(a_etr-t_mean)
            t_var     = math.sqrt(t_sum/len(t_list)) 
    
    return t_var

def getSampleVar(t_list, t_mean):

    t_sum     = 0.0
    t_var     = 0.0
    if len(t_list) > 1:
        for a_etr in t_list:
            t_sum +=(a_etr-t_mean)*(a_etr-t_mean)
            t_var     = math.sqrt(t_sum/(len(t_list)-1)) 
    
    return t_var
    

# The following several functions are used 
# for extract information of atoms, neighbour atoms

def getStdChemID(t_id):

    std_id = t_id
    if std_id.find("[") != -1:
        std_id = std_id.strip().split("[")[0].strip()
    if std_id[0].upper() !=std_id[0]:
        std_id = std_id[0].upper() + std_id[1:]

    return std_id.strip()

def getSmallFamily(a_str):

    p_str   = ""
    ch_list = []
    if a_str.find("]") != -1:
        t_strs = a_str.split("]")
        p_str = t_strs[0] + "]"
        ch_list.append(p_str)
        c_rep  = ""
        i = 0
        for a_char in t_strs[1]:
            if a_char != a_char.upper():
                c_rep += a_char
            elif a_char.isdigit():
                for i in range(int(a_char)):
                    ch_list.append(c_rep)
                c_rep  = ""
            else:
                if len(c_rep) == 0:
                    c_rep = a_char
                else:
                    ch_list.append(c_rep)
                    c_rep = a_char
        if len(c_rep):
            ch_list.append(c_rep)
        return p_str, ch_list
    else:
        c_rep  = ""
        i = 0
        for a_char in a_str:
            if a_char != a_char.upper():
                c_rep += a_char
            elif a_char.isdigit():
                for i in range(int(a_char)):
                    ch_list.append(c_rep)
                c_rep  = ""
            else:
                if len(c_rep) == 0:
                    c_rep = a_char
                else:
                    ch_list.append(c_rep)
                    c_rep = a_char
        if len(c_rep):
            ch_list.append(c_rep)
        return  ch_list[0],  ch_list

def getRootAndNB(a_str):

    ats = a_str.strip().split("(") 
    a_root = ats[0].strip()
    if a_root.find("[") != -1:
        a_root = a_root.strip().split("[")[0].strip()
    i_nb = 0
    if len(ats) > 1:
        for a_grp in ats[1:]:
            if a_grp.find(")") != -1:
               # repeat units of neighbors
               s_na = a_grp.strip().split(")")[-1]
               if len(s_na):
                   i_na = int(s_na.strip())
               i_nb = i_nb + i_na
          
    return a_root, i_nb  

def getAtomSet(a_str):
    # get atoms in a set of neighbour and next neighbour atoms

    p_str   = ""
    ch_list = []
    if a_str.find("]") != -1:
        t_strs = a_str.split("]")
        p_str = t_strs[0].strip().split("[")[0].strip() 
        ch_list.append(p_str)
        c_rep  = ""
        i = 0
        for a_char in t_strs[1]:
            if a_char != a_char.upper():
                c_rep += a_char
            elif a_char.isdigit():
                for i in range(int(a_char)):
                    ch_list.append(c_rep)
                c_rep  = ""
            else:
                if len(c_rep) == 0:
                    c_rep = a_char
                else:
                    ch_list.append(c_rep)
                    c_rep = a_char
        if len(c_rep):
            ch_list.append(c_rep)
        return ch_list
    else:
        c_rep  = ""
        i = 0
        for a_char in a_str:
            if a_char != a_char.upper():
                c_rep += a_char
            elif a_char.isdigit():
                for i in range(int(a_char)):
                    ch_list.append(c_rep)
                c_rep  = ""
            else:
                if len(c_rep) == 0:
                    c_rep = a_char
                else:
                    ch_list.append(c_rep)
                    c_rep = a_char
        if len(c_rep):
            ch_list.append(c_rep)
        return ch_list


# Chemistry related:
def organicOnly(t_alist):
 
    organic_sec = ["AT", "At", "at", "B", "b", "BR", "Br", "br", "C", "c", "CL", "Cl", "cl", 
                   "F", "f", "H", "h", "I", "i", "N","n",  "O", "o", "P", "p", "S", "s", "Se", "se"]

    for a_at in t_alist:
        if not a_at in organic_sec:
            return False
    return True

def isOrganic(t_elem):

    organic_sec = ["AT", "At", "at", "B", "b", "BR", "Br", "br", "C", "c", "CL", "Cl", "cl", 
                   "F", "f", "H", "h", "I", "i", "N","n",  "O", "o", "P", "p", "S", "s", "Se", "se"]
    if not t_elem in organic_sec:
        return False
    else:
        return True

def containMetal(t_alist):

    # metal_sec including rare-earth elememts
    metal_sec     = ['Li', 'li', 'Na', 'na', 'K',  'k',  'Rb', 'rb', 'Cs', 'cs', 'Fr', 'fr',
                     'Be', 'be', 'Mg', 'mg', 'Ca', 'ca', 'Sr', 'sr', 'Ba', 'ba', 'Ra', 'ra',
                     'Sc', 'sc', 'Y',  'y',
                     'B',  'b', 'Si',  'si', 'Ge', 'ge', 'As', 'as', 'Sb', 'sb', 'Te', 'te', 'Po', 'po', 
                     'Ti', 'ti', 'Zr', 'zr', 'Hf', 'hf', 'Rf', 'rf',
                     'V',  'v'   'Nb', 'nb', 'Ta', 'ta', 'Db', 'db', 
                     'Cr', 'cr', 'Mo', 'mo', 'W',  'w',  'Sg', 'sg', 
                     'Mn', 'mn', 'Tc', 'tc', 'Re', 're', 'Bh', 'bh',  
                     'Fe', 'fe', 'Ru', 'ru', 'Os', 'os', 'Hs', 'hs',   
                     'Co', 'co', 'Rh', 'rh', 'Ir', 'ir', 'Mt', 'mt',  
                     'Ni', 'ni', 'Pd', 'pd', 'Pt', 'pt', 'Ds', 'ds',  
                     'Cu', 'cu', 'Ag', 'ag', 'Au', 'au', 'Rg', 'rg',   
                     'Zn', 'zn', 'Cd', 'cd', 'Hg', 'hg',   
                     'Al', 'al', 'Ga', 'ga', 'In', 'in', 'Tl', 'tl', 
                     'Sn', 'sn', 'Pb', 'pb', 'Bi', 'bi', 
                     'La', 'la', 'Ce', 'ce', 'Pr', 'pr', 'Nd', 'nd', 
                     'Pm', 'pm', 'Sm', 'sm', 'Eu', 'eu', 'Gd', 'gd',
                     'Tb', 'tb', 'Dy', 'dy', 'Ho', 'ho', 'Er', 'er',
                     'Tm', 'tm', 'Yb', 'yb', 'Lu', 'lu', 'Ac', 'ac', 
                     'Th', 'th', 'Pa', 'pa', 'U',  'u',  'Np', 'np',
                     'Pu', 'pu', 'Am', 'am', 'Cm', 'cm', 'Bk', 'bk',
                     'Cf', 'cf', 'Es', 'es', 'Fm', 'fm', 'Md', 'md',
                     'No', 'no', 'Lr', 'lr']

    for a_at in t_alist:
        if  a_at in metal_sec:
            return True
    return False
 
def isMetal(t_elem):

    metal_sec     = ['Li', 'li', 'Na', 'na', 'K',  'k',  'Rb', 'rb', 'Cs', 'cs', 'Fr', 'fr',
                     'Be', 'be', 'Mg', 'mg', 'Ca', 'ca', 'Sr', 'sr', 'Ba', 'ba', 'Ra', 'ra',
                     'Sc', 'sc', 'Y',  'y',
                     'B',  'b',  'Si', 'si', 'Ge', 'ge', 'As', 'as', 'Sb', 'sb', 'Te', 'te', 'Po', 'po', 
                     'Ti', 'ti', 'Zr', 'zr', 'Hf', 'hf', 'Rf', 'rf',
                     'V',  'v'   'Nb', 'nb', 'Ta', 'ta', 'Db', 'db', 
                     'Cr', 'cr', 'Mo', 'mo', 'W',  'w',  'Sg', 'sg', 
                     'Mn', 'mn', 'Tc', 'tc', 'Re', 're', 'Bh', 'bh',  
                     'Fe', 'fe', 'Ru', 'ru', 'Os', 'os', 'Hs', 'hs',   
                     'Co', 'co', 'Rh', 'rh', 'Ir', 'ir', 'Mt', 'mt',  
                     'Ni', 'ni', 'Pd', 'pd', 'Pt', 'pt', 'Ds', 'ds',  
                     'Cu', 'cu', 'Ag', 'ag', 'Au', 'au', 'Rg', 'rg',   
                     'Zn', 'zn', 'Cd', 'cd', 'Hg', 'hg',   
                     'Al', 'al', 'Ga', 'ga', 'In', 'in', 'Tl', 'tl', 
                     'Sn', 'sn', 'Pb', 'pb', 'Bi', 'bi', 
                     'La', 'la', 'Ce', 'ce', 'Pr', 'pr', 'Nd', 'nd', 
                     'Pm', 'pm', 'Sm', 'sm', 'Eu', 'eu', 'Gd', 'gd',
                     'Tb', 'tb', 'Dy', 'dy', 'Ho', 'ho', 'Er', 'er',
                     'Tm', 'tm', 'Yb', 'yb', 'Lu', 'lu', 'Ac', 'ac', 
                     'Th', 'th', 'Pa', 'pa', 'U',  'u',  'Np', 'np',
                     'Pu', 'pu', 'Am', 'am', 'Cm', 'cm', 'Bk', 'bk',
                     'Cf', 'cf', 'Es', 'es', 'Fm', 'fm', 'Md', 'md',
                     'No', 'no', 'Lr', 'lr']

    if  t_elem in metal_sec:
         return True
    else:
         return False

class Metals:

    """This class is just a list of metal elements """

    def __init__(self):
     
        self.elements = ['Li', 'Na', 'K',  'Rb', 'Cs', 'Fr', 
                 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra', 
                 'Sc', 'Y', 
                 'B',  'Si', 'Ge', 'As', 'Sb', 'Te', 'Po', 
                 'Ti', 'Zr', 'Hf', 'Rf',
                 'V',  'Nb', 'Ta', 'Db', 
                 'Cr', 'Mo', 'W',  'Sg', 
                 'Mn', 'Tc', 'Re', 'Bh',   
                 'Fe', 'Ru', 'Os', 'Hs',    
                 'Co', 'Rh', 'Ir', 'Mt',   
                 'Ni', 'Pd', 'Pt', 'Ds',   
                 'Cu', 'Ag', 'Au', 'Rg',    
                 'Zn', 'Cd', 'Hg',   
                 'Al', 'Ga', 'In', 'Tl',  
                 'Sn', 'Pb', 'Bi',  
                 'La', 'Ce', 'Pr', 'Nd',  
                 'Pm', 'Sm', 'Eu', 'Gd', 
                 'Tb', 'Dy', 'Ho', 'Er', 
                 'Tm', 'Yb', 'Lu', 'Ac',  
                 'Th', 'Pa', 'U',  'Np', 
                 'Pu', 'Am', 'Cm', 'Bk', 
                 'Cf', 'Es', 'Fm', 'Md', 
                 'No', 'Lr' ]

 
def notSingleEl(t_str):

    # This check if atoms in a molecule with atoms are all the same element type such as O, Br, Cl etc 
    sameEleList = ["O", "Br", "F", "H", "Cl", "I", "At"]
      
    notSing = True
    atoms   = getAtomSet(t_str)
    diff_ats = []
    for a_at in atoms:
        if not a_at in diff_ats:
            diff_ats.append(a_at)
    len_ats = len(diff_ats)
    if len_ats > 1:
        notSing = False
    elif len_ats ==0:
        print "No atoms in %s ? error!"%t_str
        sys.exit(1) 
  
    return notSing         

def SingleEl(t_str):

    # This deletes all molecules with atoms of only the same element type such as O, Br, Cl etc
    sameEleList = ["C", "O", "N", "S", "P", "B", "Br", "F", "H", "Cl", "I", "At"]

    notSing = True
    ats = t_str.split("(")
    if len(ats) > 0:
         if ats[0].find("[") != -1:
             a_gran_p = ats[0].strip().split("[")[0]
         else:
             a_gran_p = ats[0].strip()
         if not a_gran_p in sameEleList:
             return False
         for a_grp in ats[1:]:
             a_p = ""
             ch_s = []
             if a_grp.find(")") != -1:
                 # repeat units of neighbors
                 unit_cont, s_na = a_grp.strip().split(")")
                 # each unit
                 a_p, ch_s = getSmallFamily(unit_cont)
                 for a_el in ch_s:
                     if a_el.strip() != a_gran_p:
                         return False
    return notSing

def outTempPDB(t_name, t_geo_name, t_atoms):

         try :
             t_file = open(t_name, "w")
         except IOError :
             print "%s can not be opened for reading "%t_name
             sys.exit(1)
         else :
              t_file.write("%s%s%s\n"%("HEADER".ljust(20), "Coordination Geometry: ".ljust(30), t_geo_name.ljust(30)))
              t_file.write("%s\n"%("CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1".ljust(80)))
              i=1
              str_i = str(i)
              for a_key in sorted(t_atoms.iterkeys()):
                  a_atom= t_atoms[a_key]
                  if i==1:
                      id = "M"
                  else:
                      str_i1 = str(i-1)
                      id = "L"+str_i1
                  x= "%5.4f"%a_atom[0]
                  y= "%5.4f"%a_atom[1]
                  z= "%5.4f"%a_atom[2]
                  t_file.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n"%("ATOM".ljust(6), str_i.rjust(5)," ", id.rjust(3), " ",
                                                     " XXX".ljust(4), " A".ljust(2), "1".rjust(4), "    ".ljust(4), 
                                                     x.rjust(8), y.rjust(8), z.rjust(8), "1.00".rjust(6), 
                                                     "20.00".rjust(6), id.rjust(12), " ".rjust(2)))
                  i=i+1    
                  str_i = str(i)
                  
              t_file.close()
