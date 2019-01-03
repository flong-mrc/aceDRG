#!/usr/bin/env ccp4-python

###########################################
#
# Dictionary Accumulator
# 20/12/2018
#
# Authors: Rob Nicholls and Garib Murshudov
#
###########################################

import sys
from gemmi import cif

PRINT_CONTENT_TO_COUT = False
CREATE_NEW_CIF = True
DISPLAY_FINAL_METADATA = False

DUPLICATE_MODIFICATION_FAIL = False
DUPLICATE_LINK_FAIL = False

# Outstanding issues:
# (1) Need option to perform graph matching instead of checking atomic nomenclature
# (2) Need to consider whether we need to deal with the case where comp_id != three_letter_code. In order to deal with this, need to know how to extract row with given three_letter_code instead of indexing by comp_id.
# (3) Doing so will fix another issue, which is the assumption that the first column in the data_comp_XXX loops is always comp_id - need to remove this line: "row = category.find_row(id)"
# (4) Assumption is that the correct data block is always "data_comp_XXX", where XXX is the three_letter_code (not the comp_id). This assumption may or may not be reasonable.
# (5) Need to use ARGPARSE for interpreting arguments.
# (6) Need to remove limitations on column tags. Current assumption is that all CIF dictionaries that are to be combined have the same columns/tags in the data_comp_list. At the minute, column tags are hard coded. Need to test using input CIF files generated from different sources.

def check_identical(block_id1,block_id2,cif1,cif2):
   data_block1 = cif1.find_block(block_id1)
   data_block2 = cif2.find_block(block_id2)
   if data_block1.get_mmcif_category_names() != data_block2.get_mmcif_category_names():
      return False
   for loop_name in data_block1.get_mmcif_category_names():
      category1 = data_block1.find_mmcif_category(loop_name)
      category2 = data_block2.find_mmcif_category(loop_name)
      if len(category1) != len(category2):
         return False
      for i in range(1,len(category1)):
         if tuple(category1[i]) != tuple(category2[i]):
            return False
   return True

try:
   if PRINT_CONTENT_TO_COUT:
      for path in sys.argv[1:]:
         doc = cif.read_file(path)
         for block in doc:
            print("BLOCK: ",block.name)
            for loop in block.get_mmcif_category_names():
               print("LOOP: ",loop)
               category = block.find_mmcif_category(loop)
               print("COLS: ")
               for column_id in category.tags:
                  for row in block.find_loop(column_id):
                     print(column_id,row)
               print("ROWS:")
               for row in category:
                  print(row)
               print("")
            print("--- END OF BLOCK ---")
         print("--- END OF FILE ---")

   if CREATE_NEW_CIF:
      print("#######################################")
      print("###### Accumulating dictionaries ######")
      print("#######################################")
      
      # Read input CIF files
      input_cif_paths = sys.argv[1:]
      data = []
      for path in input_cif_paths:
         data.append({"name":path,"data":cif.read_file(path),"comp":[],"mod":[],"link":[]})

      chemcomp_list = []   # lists comps to be included, as well as reference to cif data blocks

      USE_COMPLIST = 0
      USE_MODLIST = 0
      USE_LINKLIST = 0

      for el in data:
         input_cif = el["data"]
         INVALID_CIF = 1
         if input_cif.find_block("comp_list"):
            if len(input_cif.find_block("comp_list").get_mmcif_category_names()) > 0:
               USE_COMPLIST = 1
               INVALID_CIF = 0
         if input_cif.find_block("mod_list"):
            if len(input_cif.find_block("mod_list").get_mmcif_category_names()) > 0:
               USE_MODLIST = 1
               INVALID_CIF = 0
         if input_cif.find_block("link_list"):
            if len(input_cif.find_block("link_list").get_mmcif_category_names()) > 0:
               USE_LINKLIST = 1
               INVALID_CIF = 0
         if INVALID_CIF:
            raise Exception("Input CIF does not contain any comp_list mod_list or link_list blocks: {}".format(el["name"]))

      new_cif = cif.Document()

      if USE_COMPLIST:
         complist_block = new_cif.add_new_block("comp_list")
         complist_loop = complist_block.init_mmcif_loop("_chem_comp",["id","three_letter_code","name","group","number_atoms_all","number_atoms_nh","desc_level"])
         chemcomp_category = complist_block.find_mmcif_category("_chem_comp")
         # Process comp_list from all input CIFs, and get list of all code-CIF pairs that will be combined
         for el in data:
            input_cif = el["data"]
            block = input_cif.find_block("comp_list")
            if block:
               if len(block.get_mmcif_category_names()) == 0:
                  continue
               category = block.find_mmcif_category("_chem_comp")
               if not (set(category.tags) - set(chemcomp_category.tags)):
                  # in future, we may want to expand table rather than requiring identical tags (columns)
                  for (id,three_letter_code) in block.find(['_chem_comp.id','_chem_comp.three_letter_code']):
                     if id != three_letter_code:
                        raise Exception("Comp ID not equal to three letter code: {} {}".format(id,three_letter_code))
                     IS_DUPLICATE = 0
                     for (new_id,new_three_letter_code) in complist_block.find(['_chem_comp.id','_chem_comp.three_letter_code']):
                        if new_id != new_three_letter_code:
                           raise Exception("Comp ID not equal to three letter code: {} {}".format(new_id,new_three_letter_code))
                        if id==new_id:
                           IS_DUPLICATE = 1
                           break
                     if IS_DUPLICATE:
                        for item in chemcomp_list:
                           if item[0] == id:
                              item.append(input_cif)
                        print("Duplicate code encountered: {}. Only the first instance will be retained.".format(id))
                     else:
                        row = category.find_row(id)
                        complist_loop.add_row(row)
                        chemcomp_list.append([id,input_cif])
                     el["comp"].append(id)
               else:
                  raise Exception("Incompatible comp_list._chem_comp tags")

      if USE_LINKLIST:
         linklist_block = new_cif.add_new_block("link_list")
         linklist_loop = linklist_block.init_mmcif_loop("_chem_link",["id","comp_id_1","mod_id_1","group_comp_1","comp_id_2","mod_id_2","group_comp_2","name"])
         chemlink_category = linklist_block.find_mmcif_category("_chem_link")
         for el in data:
            input_cif = el["data"]
            block = input_cif.find_block("link_list")
            if block:
               if len(block.get_mmcif_category_names()) == 0:
                  continue
               category = block.find_mmcif_category("_chem_link")
               if not (set(category.tags) - set(chemlink_category.tags)):
                  for id in block.find_loop('_chem_link.id'):
                     row = category.find_row(id)
                     el["link"].append({"id":id,"row":row})
               else:
                  print(chemlink_category.tags)
                  print(category.tags)
                  raise Exception("Incompatible link_list._chem_link tags")

      if USE_MODLIST:
         modlist_block = new_cif.add_new_block("mod_list")
         modlist_loop = modlist_block.init_mmcif_loop("_chem_mod",["id","name","comp_id","group_id"])
         chemmod_category = modlist_block.find_mmcif_category("_chem_mod")
         for el in data:
            input_cif = el["data"]
            block = input_cif.find_block("mod_list")
            if block:
               if len(block.get_mmcif_category_names()) == 0:
                  continue
               category = block.find_mmcif_category("_chem_mod")
               if not (set(category.tags) - set(chemmod_category.tags)):
                  for id in block.find_loop('_chem_mod.id'):
                     IS_DUPLICATE = 0
                     for new_id in modlist_block.find_loop('_chem_mod.id'):
                        if id==new_id:
                           IS_DUPLICATE = 1
                           break
                     row = category.find_row(id)
                     if IS_DUPLICATE:
                        if DUPLICATE_MODIFICATION_FAIL:
                           raise Exception("Duplicate modification encountered: {}".format(id))
                     el["mod"].append({"id":id,"row":row})
               else:
                  print(chemmod_category.tags)
                  print(category.tags)
                  raise Exception("Incompatible mod_list._chem_mod tags")

# --- COMPONENTS ---

      if USE_COMPLIST:
         # Perform checks:
         # (1) Atom IDs are unique
         # (2) Duplications are compatible
         #     At present just checks for same atom_ids, but in future there will be an option for graph matching
         for code_cif in chemcomp_list:
            data_block_name = "comp_{}".format(code_cif[0])
            block1 = code_cif[1].find_block(data_block_name)
            atom_ids1 = block1.find_loop('_chem_comp_atom.atom_id')
            atom_set1 = set(atom_ids1)
            if len(atom_ids1) != len(atom_set1):
               raise Exception("Non-unique atom_id found in comp_id: {}".format(code_cif[0]))
            for i in range(2,len(code_cif)):
               block2 = code_cif[i].find_block(data_block_name)
               atom_ids2 = block2.find_loop('_chem_comp_atom.atom_id')
               atom_set2 = set(atom_ids2)
               if atom_set1 != atom_set2:
                  raise Exception("Different atomic compositions in instances of duplicate comp_id: {}".format(code_cif[0]))

         # Add all data_comp_ blocks
         for code_cif in chemcomp_list:
            print("Adding component to dictionary: {}".format(code_cif[0]))
            data_block_name = "comp_{}".format(code_cif[0])
            data_block = code_cif[1].find_block(data_block_name)
            if(data_block):
               new_block = new_cif.add_new_block(data_block_name)
               for loop in data_block.get_mmcif_category_names():
                  category = data_block.find_mmcif_category(loop)
                  new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
                  for row in category:
                     new_loop.add_row(row)
            else:
               raise Exception("Input CIF does not contain {} block".format(data_block_name))

# --- MODIFICATIONS ---

      if USE_MODLIST:

         def is_integer(x):
            try:
               int(x)
               return True
            except ValueError:
               return False

         def get_modified_code(code):
            arr = code.split('-')
            if len(arr) > 1:
               if(is_integer(arr[-1])):
                  arr[-1] = str(int(arr[-1]) + 1)
                  return "-".join(arr)
            return code+"-1"

         chemmod_list_renamed = []
         for el in data:
            mods = el["mod"]
            for i in range(0,len(mods)):
               NEW_ENTRY = True
               IS_UNIQUE = True
               mod_id = mods[i]["id"]
               for entry in chemmod_list_renamed:
                  if entry[2] == mod_id:
                     NEW_ENTRY = False
                     if check_identical("mod_{}".format(entry[0]),"mod_{}".format(mod_id),entry[1],el["data"]):
                        print("Skipping identical modification: {}".format(entry[2]))
                        IS_UNIQUE = False
                        break
               if NEW_ENTRY:
                  chemmod_list_renamed.append([mod_id,el["data"],mod_id,mods[i]["row"]])
               elif IS_UNIQUE:
                  new_name = get_modified_code(mod_id)
                  IS_VALID = False
                  while not IS_VALID:
                     IS_VALID = True
                     for entry in chemmod_list_renamed:
                        if new_name == entry[2]:
                           new_name = get_modified_code(entry[2])
                           IS_VALID = False
                           break
                  print("Renaming duplicated modification name {} to: {}".format(mod_id,new_name))
                  chemmod_list_renamed.append([mod_id,el["data"],new_name,mods[i]["row"]])
                  mods[i]["new_id"] = new_name

         # Add row to data_mod_list
         for code_cif in chemmod_list_renamed:
            code_cif[3][0] = code_cif[2]
            modlist_loop.add_row(code_cif[3])

         # Add all data_mod_ blocks
         for code_cif in chemmod_list_renamed:
            print("Adding modification to dictionary: {}".format(code_cif[2]))
            data_block_name = "mod_{}".format(code_cif[0])
            data_block_name_new = "mod_{}".format(code_cif[2])
            data_block = code_cif[1].find_block(data_block_name)
            if(data_block):
               new_block = new_cif.add_new_block(data_block_name_new)
               for loop in data_block.get_mmcif_category_names():
                  category = data_block.find_mmcif_category(loop)
                  new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
                  if code_cif[0] != code_cif[2]:   # Modify renamed mod_id in table
                     idx = -1
                     for i in range(0,len(category.tags)):
                        if str(category.tags[i].split('.')[-1]) == 'mod_id':
                           idx = i
                           break
                     if idx >= 0:
                        for row in category:
                           row[idx] = code_cif[2]
                  for row in category:
                     new_loop.add_row(row)
            else:
               raise Exception("Input CIF does not contain {} block".format(data_block_name))

# --- LINKS ---

      if USE_LINKLIST:

         if USE_MODLIST:
            # Update links with new modification IDs
            for el in data:
               for link in el["link"]:
                  for mod in el["mod"]:
                     for i in [2,5]:
                        if link["row"][i] == mod["id"]:
                           if "new_id" in mod:
                              link["row"][i] = mod["new_id"]

         chemlink_list = []
         for el in data:
            for link in el["link"]:
               NEW_LINK = True
               for existing_link in chemlink_list:
                  if link["id"] == existing_link["id"]:
                     link_block_id = "link_{}".format(link["id"])
                     if check_identical(link_block_id,link_block_id,el["data"],existing_link["data"]["data"]):
                        if list(link["row"]) == list(existing_link["link"]["row"]):
                           print("Skipping identical link: {} from {}".format(link["id"],el["name"]))
                           NEW_LINK = False
                           break
                        print(" ".join(list(link["row"])))
                        print(" ".join(list(existing_link["link"]["row"])))
                        raise Exception("Identical links encountered with inconsistent renamed modifications.")
                     raise Exception("Non-identical links encountered with the same ID: {}".format(link["id"]))
               if NEW_LINK:
                  chemlink_list.append({"id":link["id"],"data":el,"link":link})

#for link in chemlink_list:
#            print("LINK: {} {} {}".format(link["id"],link["data"]["name"],link["data"]["data"]))
#            print("  {}".format(link["link"]["row"]))

         # Add rows to the link_list block and add all data_mod_ blocks
         for link in chemlink_list:
            print("Adding link to dictionary: {}".format(link["id"]))
            linklist_loop.add_row(link["link"]["row"])
            data_block_name = "link_{}".format(link["id"])
            data_block = link["data"]["data"].find_block(data_block_name)
            if(data_block):
               new_block = new_cif.add_new_block(data_block_name)
               for loop in data_block.get_mmcif_category_names():
                  category = data_block.find_mmcif_category(loop)
                  new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
                  for row in category:
                     new_loop.add_row(row)
            else:
               raise Exception("Input CIF does not contain {} block".format(data_block_name))

# --- OUTPUT ---

      if DISPLAY_FINAL_METADATA:
         for el in data:
            print("DATA:")
            for obj in el:
               print("  {} : {}".format(obj,el[obj]))

      print("Files combined")
      fileout = 'output.cif'
      new_cif.write_file(fileout)
      print("Output written to: %s" % fileout)

except Exception as e:
   print("Error: %s" % e)
   print("Cannot continue - program terminated.")
   sys.exit(1)

print("Completed.")
