from src.download.modules import *


def latest_catalog_reader():
    files_of_current_directory = os.listdir(current_directory)
    paths_to_ls_lR = list()
    paths_to_xml = list()

    for file_in_current_dir in files_of_current_directory:
        if file_in_current_dir.startswith('ls-lR'):
            file_in_current_dir = current_directory + "/" + file_in_current_dir + "/" + 'ls-lR'
            paths_to_ls_lR.append(file_in_current_dir)
        if file_in_current_dir.startswith('xml'):
            file_in_current_dir = current_directory + "/" + file_in_current_dir + "/" + 'xml'
            paths_to_xml.append(file_in_current_dir)

    paths_to_ls_lR_sorted = sorted(paths_to_ls_lR, reverse=True)
    paths_to_xml_sorted = sorted(paths_to_xml, reverse=True)

    try:
        path_to_the_latest_listing = paths_to_ls_lR_sorted[0]

        df_catalog_the_latest_listing = pd.read_csv(Path(path_to_the_latest_listing),
                                                    names=["1", "2", "3", "4", "Data_size", "Month", "Day", "Time",
                                                           "file_name", "10", "file_names_path"], sep="\s+",
                                                    low_memory=False)

        # mmCIF
        df_catalog_the_latest_mmCIF_listing_dropna = df_catalog_the_latest_listing.dropna()
        df_catalog_the_latest_mmCIF_listing_dropna_cif_gz = df_catalog_the_latest_mmCIF_listing_dropna[
            df_catalog_the_latest_mmCIF_listing_dropna['file_name'].str.endswith('cif.gz')]
        df_catalog_the_latest_mmCIF_listing_dropna_cif_gz_34kb = df_catalog_the_latest_mmCIF_listing_dropna_cif_gz[
            df_catalog_the_latest_mmCIF_listing_dropna_cif_gz.Data_size == 34.0]

        all_mmCIF_files = list()
        for mmCIF_file_name in df_catalog_the_latest_mmCIF_listing_dropna_cif_gz_34kb["file_name"]:
            all_mmCIF_files.append(mmCIF_file_name)

        # PDB
        df_catalog_the_latest_PDB_listing_dropna = df_catalog_the_latest_listing.dropna()
        df_catalog_the_latest_PDB_listing_dropna_ent_gz = df_catalog_the_latest_PDB_listing_dropna[
            df_catalog_the_latest_PDB_listing_dropna['file_name'].str.endswith('ent.gz')]
        df_catalog_the_latest_PDB_listing_dropna_ent_gz_34kb = df_catalog_the_latest_PDB_listing_dropna_ent_gz[
            df_catalog_the_latest_PDB_listing_dropna_ent_gz.Data_size == 35.0]

        all_PDB_files = list()
        for PDB_file_name in df_catalog_the_latest_PDB_listing_dropna_ent_gz_34kb["file_name"]:
            all_PDB_files.append(PDB_file_name)

        # SIFTS
        path_to_the_latest_list = paths_to_xml_sorted[0]
        df_catalog_the_latest_listing = pd.read_csv(path_to_the_latest_list, header=None, sep="\s+", low_memory=False)
        all_SIFTS_files = set()
        for key in df_catalog_the_latest_listing:
            try:
                if df_catalog_the_latest_listing[key].str.contains(".xml.gz").any():
                    for val in df_catalog_the_latest_listing[key]:
                        if val.endswith(".xml.gz"):
                            if "/" in val:
                                all_SIFTS_files.add(val.split("/")[-1])
                            else:
                                all_SIFTS_files.add(val)
            except AttributeError:
                pass

    except IndexError:
        print("Sorry, nothing to read from. Try catalog_downloader() command first.")
        all_mmCIF_files = None
        all_PDB_files = None
        all_SIFTS_files = None

    return [all_mmCIF_files, all_PDB_files, all_SIFTS_files]


"""usage of latest_catalog_reader()"""
# all_data_list_reader = latest_catalog_reader()
# all_mmCIF_file_names_files_from_latest_catalog = all_data_list_reader[0]
# all_PDB_file_names_files_from_latest_catalog = all_data_list_reader[1]
# all_SIFTS_file_names_files_from_latest_catalog = all_data_list_reader[2]

# print(len(all_mmCIF_file_names_files_from_latest_catalog))
# print(len(all_PDB_file_names_files_from_latest_catalog))
# print(len(all_SIFTS_file_names_files_from_latest_catalog))
