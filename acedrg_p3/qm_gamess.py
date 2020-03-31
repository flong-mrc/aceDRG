# Python script
#
#
#     Copyright (C) 2014 --- 2019 Fei Long,  G. Murshudov
#
#     This code is distributed under the terms and conditions of the
#     CCP4 Program Suite Licence Agreement as a CCP4 Library.
#
#====================================================================
## The date of last modification: 01/08/2017
#

from __future__ import absolute_import
import os,os.path,sys
import platform
import glob,shutil
import re,string
import time
import math

from rdkit      import rdBase

from rdkit      import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdmolops
from rdkit.Chem import Pharm3D 
from rdkit.Chem.Pharm3D import EmbedLib
from rdkit.Geometry import rdGeometry 

if os.name != 'nt':
    import fcntl
import signal

from .exebase       import CExeCode

from .acedrgRDKit   import AcedrgRDKit
from .filetools     import FileTransformer
from .filetools     import Ccp4MmCifObj
from .chem          import ChemCheck

from .utility       import listComp
from .utility       import listComp2
from .utility       import listCompDes
from .utility       import listCompAcd

