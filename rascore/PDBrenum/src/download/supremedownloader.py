from src.download.modules import *
from src.download import lefttodownload, catalogdownloader, lookfilesinside, latestcatreader
from src.download.downloadwithThreadPool import run_downloads_with_ThreadPool, url_formation_for_pool, download_pdb_assemblies_list_with_lxml


# default_input_path_to_mmCIF = current_directory + "/mmCIF"
# default_input_path_to_PDB = current_directory + "/PDB"
# default_input_path_to_SIFTS = current_directory + "/SIFTS"
# default_output_path_to_mmCIF = current_directory + "/output_mmCIF"
# default_output_path_to_PDB = current_directory + "/output_PDB"
#
# default_input_path_to_mmCIF_assemblies = current_directory + "/mmCIF_assembly"
# default_input_path_to_PDB_assemblies = current_directory + "/PDB_assembly"
# default_output_path_to_mmCIF_assemblies = current_directory + "/output_mmCIF_assembly"
# default_output_path_to_PDB_assemblies = current_directory + "/output_PDB_assembly"


def supreme_download_master(format_of_db, job_type=None,
                            default_input_path_to_mmCIF=current_directory + "/mmCIF",
                            default_input_path_to_PDB=current_directory + "/PDB",
                            default_input_path_to_SIFTS=current_directory + "/SIFTS",
                            default_input_path_to_mmCIF_assembly=current_directory + "/mmCIF_assembly",
                            default_input_path_to_PDB_assembly=current_directory + "/PDB_assembly",

                            default_output_path_to_mmCIF=current_directory + "/output_mmCIF",
                            default_output_path_to_PDB=current_directory + "/output_PDB",
                            default_output_path_to_mmCIF_assemblies=current_directory + "/output_mmCIF_assembly",
                            default_output_path_to_PDB_assemblies=current_directory + "/output_PDB_assembly"):

    catalogdownloader.catalog_downloader()

    if job_type == "refresh":
        if os.path.exists(default_input_path_to_SIFTS):
            shutil.rmtree(default_input_path_to_SIFTS)
        if format_of_db == "mmCIF":
            if os.path.exists(default_input_path_to_mmCIF):
                shutil.rmtree(default_input_path_to_mmCIF)
            if os.path.exists(default_output_path_to_mmCIF):
                shutil.rmtree(default_output_path_to_mmCIF)

        if format_of_db == "mmCIF_assembly":
            if os.path.exists(default_input_path_to_mmCIF_assembly):
                shutil.rmtree(default_input_path_to_mmCIF_assembly)
            if os.path.exists(default_output_path_to_mmCIF_assemblies):
                shutil.rmtree(default_output_path_to_mmCIF_assemblies)

        if format_of_db == "PDB":
            if os.path.exists(default_input_path_to_PDB):
                shutil.rmtree(default_input_path_to_PDB)
            if os.path.exists(default_output_path_to_PDB):
                shutil.rmtree(default_output_path_to_PDB)

        if format_of_db == "PDB_assembly":
            if os.path.exists(default_input_path_to_PDB_assembly):
                shutil.rmtree(default_input_path_to_PDB_assembly)
            if os.path.exists(default_output_path_to_PDB_assemblies):
                shutil.rmtree(default_output_path_to_PDB_assemblies)

        if format_of_db == "all":
            if os.path.exists(default_input_path_to_PDB):
                shutil.rmtree(default_input_path_to_PDB)
            if os.path.exists(default_input_path_to_mmCIF):
                shutil.rmtree(default_input_path_to_mmCIF)
            if os.path.exists(default_input_path_to_PDB_assembly):
                shutil.rmtree(default_input_path_to_PDB_assembly)
            if os.path.exists(default_input_path_to_mmCIF_assembly):
                shutil.rmtree(default_input_path_to_mmCIF_assembly)

            if os.path.exists(default_output_path_to_mmCIF):
                shutil.rmtree(default_output_path_to_mmCIF)
            if os.path.exists(default_output_path_to_mmCIF_assemblies):
                shutil.rmtree(default_output_path_to_mmCIF_assemblies)
            if os.path.exists(default_output_path_to_PDB):
                shutil.rmtree(default_output_path_to_PDB)
            if os.path.exists(default_output_path_to_PDB_assemblies):
                shutil.rmtree(default_output_path_to_PDB_assemblies)

    if format_of_db == "mmCIF":
        all_data_from_catreader = latestcatreader.latest_catalog_reader()
        all_mmCIF_files_from_latest_catalog = all_data_from_catreader[0]
        all_SIFTS_files_from_latest_catalog = all_data_from_catreader[2]

        input_mmCIF_files_were_found = lookfilesinside.look_what_is_inside("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF)
        left_to_download_mmCIF = lefttodownload.what_is_left_to_download(input_mmCIF_files_were_found, all_mmCIF_files_from_latest_catalog)
        urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF", left_to_download_mmCIF, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
        run_downloads_with_ThreadPool("mmCIF", urls_to_target_mmCIF_files, default_input_path_to_mmCIF=default_input_path_to_mmCIF)

        input_SIFTS_files_were_found = lookfilesinside.look_what_is_inside("SIFTS", default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        left_to_download_SIFTS = lefttodownload.what_is_left_to_download(input_SIFTS_files_were_found, all_SIFTS_files_from_latest_catalog)
        urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", left_to_download_SIFTS, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        return left_to_download_mmCIF

    if format_of_db == "mmCIF_assembly":
        all_data_from_catreader = latestcatreader.latest_catalog_reader()
        all_mmCIF_files = all_data_from_catreader[0]
        all_SIFTS_files_from_latest_catalog = all_data_from_catreader[2]

        lefttodownload_mmCIF_assemblies = list()
        input_mmCIF_assembly_files_were_found = lookfilesinside.look_what_is_inside(
            "mmCIF_assembly", default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)

        all_mmCIF_files_4char = set()
        for mmCIF_file in all_mmCIF_files:
            all_mmCIF_files_4char.add(mmCIF_file[:4])

        input_mmCIF_assembly_files_were_found_4char = set()
        for mmCIF_assembly_file in input_mmCIF_assembly_files_were_found:
            input_mmCIF_assembly_files_were_found_4char.add(mmCIF_assembly_file[:4])

        set_difference = all_mmCIF_files_4char - input_mmCIF_assembly_files_were_found_4char
        list_difference = list(set_difference)

        for mmCIF_id in list_difference:
            lefttodownload_mmCIF_assemblies.append(mmCIF_id + ".cif.gz")

        urls_to_target_mmCIF_assembly_files = url_formation_for_pool("mmCIF_assembly", lefttodownload_mmCIF_assemblies)
        run_downloads_with_ThreadPool("mmCIF_assembly", urls_to_target_mmCIF_assembly_files,
                                      default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)

        input_SIFTS_files_were_found = lookfilesinside.look_what_is_inside("SIFTS", default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        left_to_download_SIFTS = lefttodownload.what_is_left_to_download(input_SIFTS_files_were_found, all_SIFTS_files_from_latest_catalog)
        urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", left_to_download_SIFTS, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        return lefttodownload_mmCIF_assemblies

    if format_of_db == "PDB":
        all_data_from_catreader = latestcatreader.latest_catalog_reader()
        all_PDB_files_from_latest_catalog = all_data_from_catreader[1]
        all_SIFTS_files_from_latest_catalog = all_data_from_catreader[2]

        input_PDB_files_were_found = lookfilesinside.look_what_is_inside("PDB", default_input_path_to_PDB=default_input_path_to_PDB)
        left_to_download_PDB = lefttodownload.what_is_left_to_download(input_PDB_files_were_found, all_PDB_files_from_latest_catalog)
        urls_to_target_PDB_files = url_formation_for_pool("PDB", left_to_download_PDB, default_input_path_to_PDB=default_input_path_to_PDB)
        run_downloads_with_ThreadPool("PDB", urls_to_target_PDB_files, default_input_path_to_PDB=default_input_path_to_PDB)

        input_SIFTS_files_were_found = lookfilesinside.look_what_is_inside("SIFTS", default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        left_to_download_SIFTS = lefttodownload.what_is_left_to_download(input_SIFTS_files_were_found, all_SIFTS_files_from_latest_catalog)
        urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", left_to_download_SIFTS, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        return left_to_download_PDB

    if format_of_db == "PDB_assembly":
        all_data_from_catreader = latestcatreader.latest_catalog_reader()
        all_SIFTS_files_from_latest_catalog = all_data_from_catreader[2]

        download_all_PDB_assemblies = download_pdb_assemblies_list_with_lxml()
        input_PDB_assembly_files_were_found = lookfilesinside.look_what_is_inside(
            "PDB_assembly", default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)
        try:
            len(download_all_PDB_assemblies)
        except TypeError:
            return print("Cannot reach https://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/all/ maybe try again later")
        lefttodownload_PDB_assemblies = [assembly for assembly in download_all_PDB_assemblies
                                         if assembly.rsplit('/', 1)[-1] not in input_PDB_assembly_files_were_found]
        run_downloads_with_ThreadPool("PDB_assembly", lefttodownload_PDB_assemblies,
                                      default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)

        input_SIFTS_files_were_found = lookfilesinside.look_what_is_inside("SIFTS", default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        left_to_download_SIFTS = lefttodownload.what_is_left_to_download(input_SIFTS_files_were_found, all_SIFTS_files_from_latest_catalog)
        urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", left_to_download_SIFTS, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        return lefttodownload_PDB_assemblies

    if format_of_db == "all":
        all_data_from_catreader = latestcatreader.latest_catalog_reader()
        all_mmCIF_files_from_latest_catalog = all_data_from_catreader[0]
        all_PDB_files_from_latest_catalog = all_data_from_catreader[1]
        all_SIFTS_files_from_latest_catalog = all_data_from_catreader[2]

        input_mmCIF_files_were_found = lookfilesinside.look_what_is_inside("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF)
        input_PDB_files_were_found = lookfilesinside.look_what_is_inside("PDB", default_input_path_to_PDB=default_input_path_to_PDB)
        input_SIFTS_files_were_found = lookfilesinside.look_what_is_inside("SIFTS", default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        left_to_download_mmCIF = lefttodownload.what_is_left_to_download(input_mmCIF_files_were_found, all_mmCIF_files_from_latest_catalog)
        left_to_download_PDB = lefttodownload.what_is_left_to_download(input_PDB_files_were_found, all_PDB_files_from_latest_catalog)
        left_to_download_SIFTS = lefttodownload.what_is_left_to_download(input_SIFTS_files_were_found, all_SIFTS_files_from_latest_catalog)

        urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF", left_to_download_mmCIF, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
        urls_to_target_PDB_files = url_formation_for_pool("PDB", left_to_download_PDB, default_input_path_to_PDB=default_input_path_to_PDB)
        urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", left_to_download_SIFTS, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        run_downloads_with_ThreadPool("mmCIF", urls_to_target_mmCIF_files, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
        run_downloads_with_ThreadPool("PDB", urls_to_target_PDB_files, default_input_path_to_PDB=default_input_path_to_PDB)
        run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        # PDB_assembly
        download_all_PDB_assemblies = download_pdb_assemblies_list_with_lxml()
        input_PDB_assembly_files_were_found = lookfilesinside.look_what_is_inside("PDB_assembly")
        try:
            len(download_all_PDB_assemblies)
        except TypeError:
            return print("Cannot reach https://ftp.wwpdb.org/pub/pdb/data/biounit/PDB/all/ maybe try again later")
        lefttodownload_PDB_assemblies = [assembly for assembly in download_all_PDB_assemblies
                                         if assembly.rsplit('/', 1)[-1] not in input_PDB_assembly_files_were_found]
        run_downloads_with_ThreadPool("PDB_assembly", lefttodownload_PDB_assemblies,
                                      default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)

        # mmCIF_assembly
        lefttodownload_mmCIF_assemblies = list()
        input_mmCIF_assembly_files_were_found = lookfilesinside.look_what_is_inside(
            "mmCIF_assembly", default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)

        all_mmCIF_files_4char = set()
        for mmCIF_file in all_mmCIF_files_from_latest_catalog:
            all_mmCIF_files_4char.add(mmCIF_file[:4])

        input_mmCIF_assembly_files_were_found_4char = set()
        for mmCIF_assembly_file in input_mmCIF_assembly_files_were_found:
            input_mmCIF_assembly_files_were_found_4char.add(mmCIF_assembly_file[:4])

        set_difference = all_mmCIF_files_4char - input_mmCIF_assembly_files_were_found_4char
        list_difference = list(set_difference)

        for mmCIF_id in list_difference:
            lefttodownload_mmCIF_assemblies.append(mmCIF_id + ".cif.gz")

        urls_to_target_mmCIF_assembly_files = url_formation_for_pool("mmCIF_assembly", lefttodownload_mmCIF_assemblies,
                                                                     default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
        run_downloads_with_ThreadPool("mmCIF_assembly", urls_to_target_mmCIF_assembly_files,
                                      default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)

        return [left_to_download_mmCIF, left_to_download_PDB, lefttodownload_mmCIF_assemblies, lefttodownload_PDB_assemblies]


"""usage of supreme_download_master"""
# supreme_download_master("mmCIF")
# supreme_download_master("PDB")
