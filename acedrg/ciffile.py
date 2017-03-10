
import os, sys
from pdbecif import startools
from collections import OrderedDict

def cifdict_from_file(file):
  st = startools.StarTokeniser()
  st.start_matching(file)
  e0 = OrderedDict()
  e1 = OrderedDict()
  e2 = OrderedDict()
  e3_list = list()
  cou = 0
  for token in st:
    tt = token.type_string
    tv = token.value
    if tt == 'COMMENT':
      pass

    elif tt == 'DATA_BLOCK':
      e1 = OrderedDict()
      e0[tv] = e1

    elif tt == 'LOOP':
      cou = 0

    elif tt == 'DATA_NAME':
      k2, sep, k3 = tv.partition('.')
      if not cou:
        e3_list = list()
        e2 = OrderedDict()
        e1[k2] = e2

      e3 = []
      e3_list.append(e3)
      e2[k3] = e3
      cou += 1

    elif tt in ('STRING', 'SQUOTE_STRING', 'DQUOTE_STRING', 'NULL'):
      if ' ' in tv:
        if tt == 'DQUOTE_STRING':
          tv = '"' + tv + '"'

        elif tt == 'SQUOTE_STRING':
          tv = "'" + tv + "'"

        else:
          assert False

      if cou >= len(e3_list):
        cou = 0

      e3_list[cou].append(tv)
      cou += 1

  for k1, e1 in e0.items():
    for k2, e2 in e1.items():
      cou = 0
      for k3, e3 in e2.items():
        if cou:
          assert cou == len(e3)

        else:
          cou = len(e3)

  return e0

def to_old_version(e0):
  mon = e0['data_comp_list']['_chem_comp']['id'][0]
  e1 = e0['data_comp_' + mon]
  e2 = e1['_chem_comp_bond']
  if 'aromatic' in e2:
    ta = zip(e2['type'], e2.pop('aromatic'))
    e2['type'] = ['aromatic' if a == 'y' else t.lower() for t, a in ta]

  else:
    e2['type'] = [t.lower() for t in e2['type']]

  e2 = e1['_chem_comp_atom']
  if 'charge' in e2:
    ke3 = list(e2.items())
    e2.clear()
    for k3, e3 in ke3:
      if k3 != 'partial_charge':
        e2['partial_charge' if k3 == 'charge' else k3] = e3

def rename_atoms(e0, rename_dict):
  mon = e0['data_comp_list']['_chem_comp']['id'][0]
  e1 = e0['data_comp_' + mon]
  for e2 in e1.values():
    changes = list()
    for k3, e3, in e2.items():
      if k3.startswith('atom_id'):
        changes.append((k3, [rename_dict.get(e4, e4) for e4 in e3]))

    e2.update(changes)

def aliases(e0):
  mon = e0['data_comp_list']['_chem_comp']['id'][0]
  e1 = e0['data_comp_' + mon]
  e2 = e1['_chem_comp_atom']
  id_dict = dict()
  alias_dict = dict()
  chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  i1 = i2 = 0
  for id, type in zip(e2['atom_id'], e2['type_symbol']):
    if i2 >= len(chars):
      i1 += 1
      i2 = 0

    alias = type + chars[i1] + chars[i2]
    id_dict[alias] = id
    alias_dict[id] = alias
    i2 += 1

  return id_dict, alias_dict

def complete_aa(e0):
  e2 = e0['data_comp_list']['_chem_comp']
  grp = e2['group'][0]
  if grp.endswith('peptide'):
    mon = e2['id'][0]
    e1 = e0['data_comp_' + mon]
    e2 = e1['_chem_comp_atom']
    if not 'OXT' in e2['atom_id']:
      for e3 in e2.values():
        e3.append('.')

      e2['comp_id'][-1] = mon
      e2['atom_id'][-1] = 'OXT'
      e2['type_symbol'][-1] = 'O'
      if 'charge' in e2:
        e2['charge'][-1] = '0'

      if 'partial_charge' in e2:
        e2['partial_charge'][-1] = '0'

      e2 = e1['_chem_comp_bond']
      for e3 in e2.values():
        e3.append('.')

      e2['comp_id'][-1] = mon
      e2['atom_id_1'][-1] = 'C'
      e2['atom_id_2'][-1] = 'OXT'
      e2['type'][-1] = 'deloc'

def cifdict_to_file(e0, file):
  with open(file, 'w') as ostream:
    for k1, e1 in e0.items():
      print >>ostream, k1
      for k2, e2 in e1.items():
        print >>ostream, 'loop_'
        e3_list = list()
        for k3, e3 in e2.items():
          print >>ostream, '.'.join((k2, k3))
          e3_list.append(e3)

        for row in zip(*e3_list):
          print >>ostream, ' '.join(row)

if __name__ == '__main__':
  e0 = cifdict_from_file(sys.argv[1])
  cifdict_to_file(e0, sys.argv[2])

