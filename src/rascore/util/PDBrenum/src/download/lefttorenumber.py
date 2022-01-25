from src.download.modules import *
from src.download.lookfilesinside import look_what_is_inside


def left_to_renumber_mmCIF(default_input_path_to_mmCIF=current_directory + "/mmCIF",
                           default_output_path_to_mmCIF=current_directory + "/output_mmCIF"):
    without_already_renumbered_mmCIF = list()
    # output_mmCIF_files_were_found_4Char = set()
    # input_mmCIF_files_were_found_4Char = set()
    output_mmCIF_files_were_found_set = set()
    input_mmCIF_files_were_found_set = set()

    mmCIF_files_were_found = look_what_is_inside("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF)
    output_mmCIF_files_were_found = look_what_is_inside("output_mmCIF", default_output_path_to_mmCIF=default_output_path_to_mmCIF)

    for output_mmCIF_file in output_mmCIF_files_were_found:
        output_mmCIF_files_were_found_set.add(output_mmCIF_file)
    for input_mmCIF_file in mmCIF_files_were_found:
        input_mmCIF_files_were_found_set.add(input_mmCIF_file)

    set_difference = input_mmCIF_files_were_found_set - output_mmCIF_files_were_found_set

    for mmCIF_file in mmCIF_files_were_found:
        if mmCIF_file in set_difference:
            without_already_renumbered_mmCIF.append(mmCIF_file)

    return without_already_renumbered_mmCIF


def left_to_renumber_PDB(default_input_path_to_PDB=current_directory + "/PDB",
                         default_output_path_to_PDB=current_directory + "/output_PDB"):
    without_already_renumbered_PDB = list()
    output_PDB_files_were_found_4Char = set()
    input_PDB_files_were_found_4Char = set()

    input_PDB_files_were_found = look_what_is_inside("PDB", default_input_path_to_PDB=default_input_path_to_PDB)
    output_PDB_files_were_found = look_what_is_inside("output_PDB", default_output_path_to_PDB=default_output_path_to_PDB)

    for output_PDB_file in output_PDB_files_were_found:
        output_PDB_files_were_found_4Char.add(output_PDB_file[:4])
    for input_PDB_file in input_PDB_files_were_found:
        input_PDB_files_were_found_4Char.add(input_PDB_file[3:7])

    set_difference = input_PDB_files_were_found_4Char - output_PDB_files_were_found_4Char
    list_difference = list(set_difference)

    for PDB_id in list_difference:
        without_already_renumbered_PDB.append("pdb" + PDB_id + ".ent.gz")

    return without_already_renumbered_PDB


""" usage of left_to_renumber_mmCIF()"""
# left_to_renumber_mmCIF()
# left_to_renumber_PDB()
