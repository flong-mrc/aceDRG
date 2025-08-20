#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:50:22 2024

@author: flong
"""
# coding in IPython notebook ;-)
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
# load molvs library
from molvs import tautomer
from molvs import standardize
from molvs import standardize_smiles

t1 = Chem.MolFromSmiles('C1=CC=CNC(=O)1')

