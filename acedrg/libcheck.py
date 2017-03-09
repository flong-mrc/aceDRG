
import os, sys, re, shutil
import subprocess as SP
from collections import namedtuple
import ciffile

def isfile():
  return lambda x: os.path.isfile(x)

def isstr():
  return lambda x: type(x) is str and bool(x)

def eqnone():
  return lambda x: x is None

def eq(y):
  t = type(y)
  if t == str:
    l = y.lower()
    return lambda x: l == x.lower()

  else:
    return lambda x: y == t(x)

default_list = (
  ('FILE_L', ''),
  ('MON', ''),
  ('FILE_PDB', ''),
  ('FILE_SMILE', ''),
  ('FILE_SDF', ''),
  ('FILE_MOL', ''),
  ('FILE_CIF', ''),
  ('FILE_CSD', ''),
  ('HFLAG', 'Y'),
  ('IND', 'N'),
  ('FILE_O', 'libcheck'),
  ('FILE_L2', ''),
  ('ANGLE', 0.0),
  ('LIST', 'M'),
  ('REF', 'Y'),
  ('TEST', 0),
  ('COOR', 'N'),
  ('LCOOR', 'Y'),
  ('NODIST', 'N'),
  ('SRCH', 'N'),
)
mode_list = (
  ('min2full',
    (
      ('NODIST', eq('Y')),
      ('REF', eq('0')),
      ('LCOOR', eq('N')),
      ('FILE_L', isfile()),
      ('MON', isstr()),
    ),
  ),
  ('lib2full',
    (
      ('REF', eq('0')),
      ('MON', isstr()),
    ),
  ),
  ('sml2full',
    (
      ('NODIST', eq('Y')),
      ('REF', eq('0')),
      ('LCOOR', eq('N')),
      ('FILE_SMILE', isfile()),
      ('MON', isstr()),
    ),
  ),
  ('mol2full',
    (
      ('NODIST', eq('Y')),
      ('REF', eq('0')),
      ('LCOOR', eq('N')),
      ('FILE_MOL', isfile()),
    ),
  ),
)

ccp4top = os.path.abspath(os.environ['CCP4'])
monomers = os.path.abspath(os.environ['CLIBD_MON'])
libcheck_exe = os.path.join(ccp4top, 'bin', 'libcheck')
acedrg_sh = os.path.join(ccp4top, 'bin', 'acedrg')

def libcheck():
  if len(sys.argv) > 1 and sys.argv[1] == '-i':
    SP.Popen((libcheck_exe, '-i')).wait()
    return

  kwd, default = zip(*default_list)
  param_dict = dict(zip(kwd, default))
  line_list = sys.stdin.readlines()
  assert line_list.pop(0).strip().lower() in 'yn'
  for line in line_list:
    line = line.strip()
    if not line.startswith('#'):
      key, sep, value = line.strip().partition(' ')
      if key:
        param_dict[key.upper()] = value.lstrip()

  assert not sep
  default_filter = dict(zip(kwd, [eq(value) for value in default]))
  pattern = namedtuple('Pattern', kwd)
  params = pattern(**param_dict)
  inp_mode = None
  for mode, changed in mode_list:
    filter = dict(default_filter)
    filter.update(changed)
    inp_mode = mode
    for x, fun in zip(params, pattern(**filter)):
      if not fun(x):
        inp_mode = None
        break

    if inp_mode:
      break

  test = len(sys.argv) > 1
  if inp_mode:
    prefix = 'acedrg'
    mon = params.MON
    if inp_mode == 'min2full':
      sys.argv[:] = (acedrg_sh, '-c', params.FILE_L, '-r', mon, '-o', prefix)

    elif inp_mode == 'sml2full':
      sys.argv[:] = (acedrg_sh, '-i', params.FILE_SMILE, '-r', mon, '-o', prefix)

    elif inp_mode == 'mol2full':
      with open(params.FILE_MOL) as istream:
        mon = istream.readline().strip()

      if not mon.isalnum():
        mon = 'UNL'

      sys.argv[:] = (acedrg_sh, '-m', params.FILE_MOL, '-r', mon, '-o', prefix)

    elif inp_mode == 'lib2full':
      inp = 'input.lib'
      cifdict = ciffile.cifdict_from_file(os.path.join(monomers, mon[0].lower(), mon + '.cif'))
      ciffile.complete_aa(cifdict)
      ciffile.cifdict_to_file(cifdict, inp)
      sys.argv[:] = (acedrg_sh, '-c', inp, '-r', mon, '-o', prefix)

    with open('acedrg.cmd', 'w') as ostream:
      print >>ostream, ' '.join(sys.argv)

    if not test:
      import acedrg
      acedrg.main()
      cifdict = ciffile.cifdict_from_file(prefix + '.cif')
      ciffile.to_old_version(cifdict)
      ciffile.cifdict_to_file(cifdict, 'libcheck.lib')

      shutil.move(prefix + '.pdb', 'libcheck_' + mon + '.pdb')
      shutil.move(prefix + '_TMP', 'acedrg_tmp')

  else:
    with open('libcheck.cmd', 'w') as ostream:
      for line in line_list:
        ostream.write(line)

    if not test:
      sp = SP.Popen(libcheck_exe, stdin=SP.PIPE)
      for line in line_list:
        sp.stdin.write(line)

      sp.wait()

if __name__ == '__main__':
  libcheck()

