from src.download.modules import *


def my_hook(t):
    last_b = [0]

    def update_to(b=1, b_size=1, t_size=None):
        if t_size is not None:
            t.total = t_size
        t.update((b - last_b[0]) * b_size)
        last_b[0] = b
    return update_to


def downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes):
    socket.setdefaulttimeout(600)
    for _ in range(10):
        try:
            last_slash = ftp_to_download.rsplit('/', 1)[-1]
            with tqdm.tqdm(unit="B", unit_scale=True, desc="Downloading mmCIF/SIFTS catalogs " + last_slash, position=0, leave=True) as t:
                reporthook = my_hook(t)
                urllib.request.urlretrieve(ftp_to_download, Path(str(where_the_file_goes) + "/" + last_slash), reporthook=reporthook)
            break
        except Exception:
            time.sleep(1)


def catalog_downloader():
    """PDB ls-lR catalog"""
    ftp_for_all_mmCIF_and_PDB = "ftp://ftp.rcsb.org/pub/pdb/data/structures/ls-lR"

    ftp_to_download = ftp_for_all_mmCIF_and_PDB
    last_slash = ftp_to_download.rsplit('/', 1)[-1]
    today_date = date.today()
    today_date_str = today_date.strftime("_%Y_%m_%d")

    where_the_file_goes = current_directory + "/" + last_slash + today_date_str
    if not os.path.exists(where_the_file_goes):
        os.makedirs(where_the_file_goes)
        downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes)
    else:
        if not os.path.isfile(where_the_file_goes + "/ls-lR"):
            downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes)

    # downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes)

    # reading txt file SIFTS parent catalog and creating pandas df out of it
    df_catalog_listing_everything = pd.read_csv(Path(str(where_the_file_goes) + "/" + last_slash),
                                                names=["1", "2", "3", "4", "Data_size", "Month", "Day", "Time",
                                                       "file_name", "10", "file_names_path"], sep="\s+",
                                                low_memory=False)

    # Dropping all unnecessary rows leaving only files with 'cif.gz' endings
    df_mmCIF_catalog_dropna = df_catalog_listing_everything.dropna()
    df_mmCIF_catalog_dropna_cif_gz = df_mmCIF_catalog_dropna[df_mmCIF_catalog_dropna['file_name'].str.endswith('cif.gz')]
    df_mmCIF_catalog_dropna_cif_gz_34kb = df_mmCIF_catalog_dropna_cif_gz[df_mmCIF_catalog_dropna_cif_gz.Data_size == 34.0]

    # Dropping all unnecessary rows leaving only files with 'ent.gz' endings
    df_PDB_catalog_dropna = df_catalog_listing_everything.dropna()
    df_PDB_catalog_dropna_ent_gz = df_PDB_catalog_dropna[df_PDB_catalog_dropna['file_name'].str.endswith('ent.gz')]
    df_PDB_catalog_dropna_ent_gz_35kb = df_PDB_catalog_dropna_ent_gz[df_PDB_catalog_dropna_ent_gz.Data_size == 35.0]

    # creating lists of the mmCIF file_names
    list_of_mmCIF_cif_gz_file_names = list()
    for mmCIF_file_name in df_mmCIF_catalog_dropna_cif_gz_34kb["file_name"]:
        list_of_mmCIF_cif_gz_file_names.append(mmCIF_file_name)

    # creating lists of the PDB file_names
    list_of_PDB_ent_gz_file_names = list()
    for PDB_file_name in df_PDB_catalog_dropna_ent_gz_35kb["file_name"]:
        list_of_PDB_ent_gz_file_names.append(PDB_file_name)

    """SIFTS xml catalog"""
    ftp_all_SIFTS = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml"
    ftp_to_download = ftp_all_SIFTS
    last_slash = ftp_to_download.rsplit('/', 1)[-1]

    where_the_file_goes = current_directory + "/" + last_slash + today_date_str
    if not os.path.exists(where_the_file_goes):
        os.makedirs(where_the_file_goes)
        downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes)
    else:
        if not os.path.isfile(where_the_file_goes + "/xml"):
            downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes)

    # downloader_for_catalog_with_urllib(ftp_to_download, where_the_file_goes)

    # reading txt file SIFTS parent catalog and creating pandas df out of it
    df_catalog_listing_everything = pd.read_csv(Path(str(where_the_file_goes) + "/" + last_slash),
                                                names=["1", "2", "3", "4", "Data_size", "Month", "Day", "Time",
                                                       "file_name", "10", "file_names_path"], sep="\s+",
                                                low_memory=False)

    # Dropping all unnecessary rows leaving only files with 'xml.gz' endings
    df_SIFTS_catalog_dropna = df_catalog_listing_everything.dropna()
    df_SIFTS_catalog_dropna_xml_gz = df_SIFTS_catalog_dropna[
        df_SIFTS_catalog_dropna['file_name'].str.endswith('xml.gz')]
    df_SIFTS_catalog_dropna_xml_gz_27kb = df_SIFTS_catalog_dropna_xml_gz[
        df_SIFTS_catalog_dropna_xml_gz.Data_size == 27.0]

    # creating lists of the SIFTS file_names
    list_of_SIFTS_xml_gz_file_names = list()
    for SIFTS_file_names in df_SIFTS_catalog_dropna_xml_gz_27kb["file_name"]:
        list_of_SIFTS_xml_gz_file_names.append(SIFTS_file_names)

    _4Char_list_of_SIFTS_xml_gz_file_names = list()
    for SIFTS_file_names_4Char in list_of_SIFTS_xml_gz_file_names:
        _4Char_list_of_SIFTS_xml_gz_file_names.append(SIFTS_file_names_4Char[:4])

    _4Char_list_of_PDB_ent_gz_file_names = list()
    for PDB_file_names_4Char in list_of_PDB_ent_gz_file_names:
        _4Char_list_of_PDB_ent_gz_file_names.append(PDB_file_names_4Char[3:7])

    _4Char_list_of_mmCIF_cif_gz_file_names = list()
    for mmCIF_file_names_4Char in list_of_mmCIF_cif_gz_file_names:
        _4Char_list_of_mmCIF_cif_gz_file_names.append(mmCIF_file_names_4Char[:4])

    df_list_of_mmCIF_cif_gz_file_names = pd.DataFrame(
        zip(list_of_mmCIF_cif_gz_file_names, _4Char_list_of_mmCIF_cif_gz_file_names), columns=["mmCIF", "4mmCIF"])
    df_list_of_PDB_ent_gz_file_names = pd.DataFrame(
        zip(list_of_PDB_ent_gz_file_names, _4Char_list_of_PDB_ent_gz_file_names), columns=["PDB", "4PDB"])
    df_list_of_SIFTS_xml_gz_file_names = pd.DataFrame(
        zip(list_of_SIFTS_xml_gz_file_names, _4Char_list_of_SIFTS_xml_gz_file_names), columns=["SIFTS", "4SIFTS"])

    merged_df_mmCIF_PDB_file_names = df_list_of_mmCIF_cif_gz_file_names.merge(df_list_of_PDB_ent_gz_file_names,
                                                                              left_on='4mmCIF', right_on='4PDB',
                                                                              how="left")
    merged_df_mmCIF_PDB_SIFTS_file_names = merged_df_mmCIF_PDB_file_names.merge(df_list_of_SIFTS_xml_gz_file_names,
                                                                                left_on='4mmCIF', right_on='4SIFTS',
                                                                                how="left")

    merged_df_mmCIF_PDB_SIFTS_file_names['SIFTS'] = merged_df_mmCIF_PDB_SIFTS_file_names['SIFTS'].replace(np.nan, "0000")
    merged_df_mmCIF_PDB_SIFTS_file_names['PDB'] = merged_df_mmCIF_PDB_SIFTS_file_names['PDB'].replace(np.nan, "0000")

    SIFTS_file_names_with_null_if_files_absent = list()
    for SIFTS_file_name_null_for_absent in merged_df_mmCIF_PDB_SIFTS_file_names['SIFTS']:
        SIFTS_file_names_with_null_if_files_absent.append(SIFTS_file_name_null_for_absent)

    PDB_file_names_with_null_if_files_absent = list()
    for PDB_file_name_null_for_absent in merged_df_mmCIF_PDB_SIFTS_file_names['PDB']:
        PDB_file_names_with_null_if_files_absent.append(PDB_file_name_null_for_absent)

    mmCIF_file_names_with_null_if_files_absent = list()
    for mmCIF_file_name_null_for_absent in merged_df_mmCIF_PDB_SIFTS_file_names['mmCIF']:
        mmCIF_file_names_with_null_if_files_absent.append(mmCIF_file_name_null_for_absent)

    return (mmCIF_file_names_with_null_if_files_absent,
            PDB_file_names_with_null_if_files_absent,
            SIFTS_file_names_with_null_if_files_absent)


"""usage of catalog_downloader()"""
# all_data_tuple = catalog_downloader()
# mmCIF_file_names_with_null_if_files_absent = all_data_tuple[0]
# PDB_file_names_with_null_if_files_absent = all_data_tuple[1]
# SIFTS_file_names_with_null_if_files_absent = all_data_tuple[2]

# print(len(mmCIF_file_names_with_null_if_files_absent))
# print(len(PDB_file_names_with_null_if_files_absent))
# print(len(SIFTS_file_names_with_null_if_files_absent))
