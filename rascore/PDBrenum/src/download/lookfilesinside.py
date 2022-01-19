from src.download.modules import *


def look_what_is_inside(format_to_look_at,
                        default_input_path_to_mmCIF=current_directory + "/mmCIF",
                        default_input_path_to_PDB=current_directory + "/PDB",
                        default_input_path_to_SIFTS=current_directory + "/SIFTS",
                        default_output_path_to_mmCIF=current_directory + "/output_mmCIF",
                        default_output_path_to_PDB=current_directory + "/output_PDB",
                        default_input_path_to_mmCIF_assembly=current_directory + "/mmCIF_assembly",
                        default_input_path_to_PDB_assembly=current_directory + "/PDB_assembly",
                        default_output_path_to_mmCIF_assembly=current_directory + "/output_mmCIF_assembly",
                        default_output_path_to_PDB_assembly=current_directory + "/output_PDB_assembly"):
    if format_to_look_at == "SIFTS":
        if not os.path.exists(default_input_path_to_SIFTS):
            os.makedirs(default_input_path_to_SIFTS)
        result = [f for f in listdir(default_input_path_to_SIFTS) if isfile(join(default_input_path_to_SIFTS, f))]
        return result
    if format_to_look_at == "mmCIF":
        if not os.path.exists(default_input_path_to_mmCIF):
            os.makedirs(default_input_path_to_mmCIF)
        result = [f for f in listdir(default_input_path_to_mmCIF) if isfile(join(default_input_path_to_mmCIF, f))]
        return result
    if format_to_look_at == "PDB":
        if not os.path.exists(default_input_path_to_PDB):
            os.makedirs(default_input_path_to_PDB)
        result = [f for f in listdir(default_input_path_to_PDB) if isfile(join(default_input_path_to_PDB, f))]
        return result
    if format_to_look_at == "output_mmCIF":
        if not os.path.exists(default_output_path_to_mmCIF):
            os.makedirs(default_output_path_to_mmCIF)
        result = [f for f in listdir(default_output_path_to_mmCIF) if isfile(join(default_output_path_to_mmCIF, f))]
        return result
    if format_to_look_at == "output_PDB":
        if not os.path.exists(default_output_path_to_PDB):
            os.makedirs(default_output_path_to_PDB)
        result = [f for f in listdir(default_output_path_to_PDB) if isfile(join(default_output_path_to_PDB, f))]
        return result
    if format_to_look_at == "mmCIF_assembly":
        if not os.path.exists(default_input_path_to_mmCIF_assembly):
            os.makedirs(default_input_path_to_mmCIF_assembly)
        result = [f for f in listdir(default_input_path_to_mmCIF_assembly) if isfile(join(default_input_path_to_mmCIF_assembly, f))]
        return result
    if format_to_look_at == "PDB_assembly":
        if not os.path.exists(default_input_path_to_PDB_assembly):
            os.makedirs(default_input_path_to_PDB_assembly)
        result = [f for f in listdir(default_input_path_to_PDB_assembly) if isfile(join(default_input_path_to_PDB_assembly, f))]
        return result
    if format_to_look_at == "output_mmCIF_assembly":
        if not os.path.exists(default_output_path_to_mmCIF_assembly):
            os.makedirs(default_output_path_to_mmCIF_assembly)
        result = [f for f in listdir(default_output_path_to_mmCIF_assembly) if isfile(join(default_output_path_to_mmCIF_assembly, f))]
        return result
    if format_to_look_at == "output_PDB_assembly":
        if not os.path.exists(default_output_path_to_PDB_assembly):
            os.makedirs(default_output_path_to_PDB_assembly)
        result = [f for f in listdir(default_output_path_to_PDB_assembly) if isfile(join(default_output_path_to_PDB_assembly, f))]
        return result


""" usage of look_what_is_inside()"""
# input_mmCIF_files_were_found = look_what_is_inside("mmCIF")
# input_PDB_files_were_found = look_what_is_inside("PDB")
# input_SIFTS_files_were_found = look_what_is_inside("SIFTS")
# output_PDB_files_were_found = look_what_is_inside("output_PDB")
# output_mmCIF_files_were_found = look_what_is_inside("output_mmCIF")

# print(len(input_mmCIF_files_were_found))
# print(len(input_PDB_files_were_found))
# print(len(input_SIFTS_files_were_found))
# print(len(output_mmCIF_files_were_found))
# print(len(output_PDB_files_were_found))
