from src.download.modules import *


def what_is_left_to_download(files_were_found, file_names_with_null_if_files_absent):
    left_to_download = list()
    df_file_names_with_null_if_files_absent = pd.DataFrame(file_names_with_null_if_files_absent, columns=['all'])
    df_files_were_found = pd.DataFrame(files_were_found, columns=['found'])

    merged_all_vs_found = df_file_names_with_null_if_files_absent.merge(df_files_were_found, how="left", left_on="all", right_on="found")
    difference_all_vs_found = merged_all_vs_found.where((merged_all_vs_found["found"] != merged_all_vs_found["all"]))

    for file_for_download in difference_all_vs_found.dropna(how="all")["all"]:
        left_to_download.append(file_for_download)

    return left_to_download


"""usage of what_is_left_to_download()"""
# left_to_download_mmCIF = what_is_left_to_download(input_mmCIF_files_were_found, all_mmCIF_file_names_files_from_latest_catalog)
# left_to_download_PDB = what_is_left_to_download(input_PDB_files_were_found, all_PDB_file_names_files_from_latest_catalog)
# left_to_download_SIFTS = what_is_left_to_download(input_SIFTS_files_were_found, all_SIFTS_file_names_files_from_latest_catalog)

# print(len(left_to_download_mmCIF))
# print(len(left_to_download_PDB))
# print(len(left_to_download_SIFTS))
