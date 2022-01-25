from src.download.modules import *


def delete_outdated_files(files_to_refresh, path_to_the_files):
    os.chdir(path_to_the_files)
    for file_for_removal in files_to_refresh:
        try:
            os.remove(file_for_removal)
        except FileNotFoundError:
            pass
    os.chdir(current_directory)


"""usage of delete_outdated_files"""
# delete_outdated_files(refresher_for_mmCIF, format_to_delete = "mmCIF")
# delete_outdated_files(refresher_for_PDB, format_to_delete = "PDB")
# delete_outdated_files(refresher_for_SIFTS, format_to_delete = "SIFTS")
