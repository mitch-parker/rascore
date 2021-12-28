#####################################################################################################################################################
# DEPENDENCIES #
#####################################################################################################################################################


# import os
# import argparse

from src.download.modules import *
from src.download.lefttorenumber import left_to_renumber_mmCIF, left_to_renumber_PDB
from src.download.inputtextfileparser import input_text_file_parser
from src.download.shortusagemessage import short_usage_messenger
from src.download.longusagemessage import long_usage_messenger
from src.download.supremedownloader import supreme_download_master
from src.download.lookfilesinside import look_what_is_inside
from src.download.downloadwithThreadPool import run_downloads_with_ThreadPool, url_formation_for_pool, download_pdb_assemblies_list_with_lxml

from src.renum.PDB.new_PDB import ProcessPool_run_renum_PDB
from src.renum.shared.write_log import log_writer
from src.renum.mmCIF.ProcessPool_run_renum_mmCIF import ProcessPool_run_renum_mmCIF
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

current_directory = os.getcwd()
exception_AccessionIDs = ["P42212", "Q17104", "Q27903", "Q93125", "P03069", "D3DLN9", "Q96UT3", "P0ABE7", "P00192",
                          "P76805", "Q8XCE3", "P00720", "Q38170", "Q94N07", "P0AEX9", "P02928", "Q2M6S0"]

# os.environ["HTTPS_PROXY"] = "http://155.247.166.25:8080"
# os.environ["HTTP_PROXY"] = "http://155.247.166.25:8080"
# os.environ["FTP_PROXY"] = "http://155.247.166.25:8080"
# os.environ["https_proxy"] = "http://155.247.166.25:8080"
# os.environ["http_proxy"] = "http://155.247.166.25:8080"
# os.environ["ftp_proxy"] = "http://155.247.166.25:8080"

# IMPORTANT!!! Check your network connection: Wi-Fi or Wired Ethernet

#####################################################################################################################################################
# ARGUMENTS #
#####################################################################################################################################################


argpar = argparse.ArgumentParser(usage=short_usage_messenger(), add_help=False)

argpar.add_argument("-h", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("--help", action="store_true", help=argparse.SUPPRESS)

argpar.add_argument("-rftf", "--renumber_from_text_file", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-rfla", "--renumber_from_list_of_arguments", metavar="6dbp 3v03 2jit", nargs="*", type=str, help=argparse.SUPPRESS)

argpar.add_argument("-dftf", "--download_from_text_file", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-dfla", "--download_from_list_of_arguments", metavar="6dbp 3v03 2jit", nargs="+", type=str, help=argparse.SUPPRESS)

argpar.add_argument("-redb", "--renumber_entire_database", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("-dall", "--download_entire_database", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("-refr", "--refresh_entire_database", action="store_true", help=argparse.SUPPRESS)

argpar.add_argument("-PDB", "--PDB_format_only", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("-mmCIF", "--mmCIF_format_only", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("-PDB_assembly", "--PDB_assembly_format_only", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("-mmCIF_assembly", "--mmCIF_assembly_format_only", action="store_true", help=argparse.SUPPRESS)
argpar.add_argument("-all", "--all_formats", action="store_true", help=argparse.SUPPRESS)

argpar.add_argument("-sipm", "--set_default_input_path_to_mmCIF", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sipma", "--set_default_input_path_to_mmCIF_assembly", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sipp", "--set_default_input_path_to_PDB", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sippa", "--set_default_input_path_to_PDB_assembly", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sips", "--set_default_input_path_to_SIFTS", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sopm", "--set_default_output_path_to_mmCIF", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sopma", "--set_default_output_path_to_mmCIF_assembly", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sopp", "--set_default_output_path_to_PDB", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-soppa", "--set_default_output_path_to_PDB_assembly", type=str, help=argparse.SUPPRESS)

argpar.add_argument("-sdmn", "--set_default_mmCIF_num", type=int, help=argparse.SUPPRESS)
argpar.add_argument("-sdpn", "--set_default_PDB_num", type=int, help=argparse.SUPPRESS)

argpar.add_argument("-nproc", "--set_number_of_processes", type=int, help=argparse.SUPPRESS)
argpar.add_argument("-offz", "--set_to_off_mode_gzip", action="store_true", help=argparse.SUPPRESS)


args = argpar.parse_args()


#####################################################################################################################################################
# FLAGS #
#####################################################################################################################################################


if args.help:
    print(long_usage_messenger())

if args.h:
    print(short_usage_messenger())

if args.set_default_input_path_to_mmCIF:
    default_input_path_to_mmCIF = args.set_default_input_path_to_mmCIF
else:
    default_input_path_to_mmCIF = current_directory + "/mmCIF"

if args.set_default_input_path_to_mmCIF_assembly:
    default_input_path_to_mmCIF_assembly = args.set_default_input_path_to_mmCIF_assembly
else:
    default_input_path_to_mmCIF_assembly = current_directory + "/mmCIF_assembly"

if args.set_default_input_path_to_PDB:
    default_input_path_to_PDB = args.set_default_input_path_to_PDB
else:
    default_input_path_to_PDB = current_directory + "/PDB"

if args.set_default_input_path_to_PDB_assembly:
    default_input_path_to_PDB_assembly = args.set_default_input_path_to_PDB_assembly
else:
    default_input_path_to_PDB_assembly = current_directory + "/PDB_assembly"

if args.set_default_input_path_to_SIFTS:
    default_input_path_to_SIFTS = args.set_default_input_path_to_SIFTS
else:
    default_input_path_to_SIFTS = current_directory + "/SIFTS"

if args.set_default_output_path_to_mmCIF:
    default_output_path_to_mmCIF = args.set_default_output_path_to_mmCIF
else:
    default_output_path_to_mmCIF = current_directory + "/output_mmCIF"

if args.set_default_output_path_to_mmCIF_assembly:
    default_output_path_to_mmCIF_assembly = args.set_default_output_path_to_mmCIF_assembly
else:
    default_output_path_to_mmCIF_assembly = current_directory + "/output_mmCIF_assembly"

if args.set_default_output_path_to_PDB:
    default_output_path_to_PDB = args.set_default_output_path_to_PDB
else:
    default_output_path_to_PDB = current_directory + "/output_PDB"

if args.set_default_output_path_to_PDB_assembly:
    default_output_path_to_PDB_assembly = args.set_default_output_path_to_PDB_assembly
else:
    default_output_path_to_PDB_assembly = current_directory + "/output_PDB_assembly"

if args.set_default_mmCIF_num:
    default_mmCIF_num = args.set_default_mmCIF_num
else:
    default_mmCIF_num = 50000

if args.set_default_PDB_num:
    default_PDB_num = args.set_default_PDB_num
else:
    default_PDB_num = 5000

if args.set_to_off_mode_gzip:
    gzip_mode = "off"
else:
    gzip_mode = "on"

if args.set_number_of_processes:
    nproc = args.set_number_of_processes
else:
    nproc = None


if __name__ == "__main__":

    #################################################################################################################################################
    # PARTIAL DB WORK #
    #################################################################################################################################################

    # RENUMBER
    # RENUMBER FROM TEXT FILE or RENUMBER FROM LIST OF ARGUMENTS
    if args.renumber_from_text_file or args.renumber_from_list_of_arguments:
        if args.renumber_from_text_file:
            parsed_input_text = (input_text_file_parser(args.renumber_from_text_file))
        else:
            parsed_input_text = args.renumber_from_list_of_arguments

        if args.all_formats:
            urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF", parsed_input_text, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            urls_to_target_mmCIF_assembly_files = url_formation_for_pool("mmCIF_assembly", parsed_input_text,
                                                                         default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            urls_to_target_PDB_files = url_formation_for_pool("PDB", parsed_input_text, default_input_path_to_PDB=default_input_path_to_PDB)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            run_downloads_with_ThreadPool("mmCIF", urls_to_target_mmCIF_files, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            run_downloads_with_ThreadPool("mmCIF_assembly", urls_to_target_mmCIF_assembly_files,
                                          default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            run_downloads_with_ThreadPool("PDB", urls_to_target_PDB_files, default_input_path_to_PDB=default_input_path_to_PDB)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            # renum PDB
            passed_as_arg_file_4Char_PDB = list()
            for file_name in parsed_input_text:
                passed_as_arg_file_4Char_PDB.append(file_name[:4])
            input_PDB_files_were_found = look_what_is_inside("PDB", default_input_path_to_PDB=default_input_path_to_PDB)
            target_files_list_PDB = list()
            for file_name in input_PDB_files_were_found:
                if file_name[3:7] in passed_as_arg_file_4Char_PDB:
                    target_files_list_PDB.append(file_name)
            ProcessPool_run_renum_PDB("PDB", target_files_list_PDB, default_input_path_to_PDB, default_input_path_to_SIFTS,
                                      default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)

            # renum mmCIF_assembly
            input_mmCIF_files_were_found = look_what_is_inside("mmCIF_assembly",
                                                               default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            passed_as_arg_file_4Char_mmCIF = list()
            for file_name in parsed_input_text:
                passed_as_arg_file_4Char_mmCIF.append(file_name[:4])
            target_files_list_mmCIF = list()
            for file_name in input_mmCIF_files_were_found:
                if file_name[:4] in passed_as_arg_file_4Char_mmCIF:
                    target_files_list_mmCIF.append(file_name)

            ProcessPool_run_renum_mmCIF("mmCIF_assembly", target_files_list_mmCIF, default_input_path_to_mmCIF_assembly,
                                        default_input_path_to_SIFTS, default_output_path_to_mmCIF_assembly, default_mmCIF_num,
                                        gzip_mode, exception_AccessionIDs, nproc)

            # renum mmCIF
            input_mmCIF_files_were_found = look_what_is_inside("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            passed_as_arg_file_4Char_mmCIF = list()
            for file_name in parsed_input_text:
                passed_as_arg_file_4Char_mmCIF.append(file_name[:4])
            target_files_list_mmCIF = list()
            for file_name in input_mmCIF_files_were_found:
                if file_name[:4] in passed_as_arg_file_4Char_mmCIF:
                    target_files_list_mmCIF.append(file_name)
            res = ProcessPool_run_renum_mmCIF("mmCIF", target_files_list_mmCIF, default_input_path_to_mmCIF, default_input_path_to_SIFTS,
                                              default_output_path_to_mmCIF, default_mmCIF_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

            # renum PDB_assembly
            if not os.path.exists(default_input_path_to_PDB_assembly):
                os.makedirs(default_input_path_to_PDB_assembly)
            urls_to_target_PDB_files = list()
            all_urls_set = download_pdb_assemblies_list_with_lxml()
            for url in all_urls_set:
                for PDB_id in parsed_input_text:
                    if PDB_id in url:
                        urls_to_target_PDB_files.append(url)
            run_downloads_with_ThreadPool("PDB_assembly", urls_to_target_PDB_files,
                                          default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)
            input_PDB_files_were_found = look_what_is_inside("PDB_assembly", default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)

            passed_as_arg_file_4Char_PDB = list()
            for file_name in parsed_input_text:
                passed_as_arg_file_4Char_PDB.append(file_name[:4])
            target_files_list_PDB = list()
            for file_name in input_PDB_files_were_found:
                if file_name[:4] in passed_as_arg_file_4Char_PDB:
                    target_files_list_PDB.append(file_name)
            ProcessPool_run_renum_PDB("PDB_assembly", target_files_list_PDB, default_input_path_to_PDB_assembly, default_input_path_to_SIFTS,
                                      default_output_path_to_PDB_assembly, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)

        elif args.PDB_assembly_format_only:
            if not os.path.exists(default_input_path_to_PDB_assembly):
                os.makedirs(default_input_path_to_PDB_assembly)
            urls_to_target_PDB_files = list()
            all_urls_set = download_pdb_assemblies_list_with_lxml()
            for url in all_urls_set:
                for PDB_id in parsed_input_text:
                    if PDB_id in url:
                        urls_to_target_PDB_files.append(url)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("PDB_assembly", urls_to_target_PDB_files,
                                          default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            input_PDB_files_were_found = look_what_is_inside("PDB_assembly", default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)
            passed_as_arg_file_4Char_PDB = list()
            for file_name in parsed_input_text:
                passed_as_arg_file_4Char_PDB.append(file_name[:4])
            target_files_list_PDB = list()
            for file_name in input_PDB_files_were_found:
                if file_name[:4] in passed_as_arg_file_4Char_PDB:
                    target_files_list_PDB.append(file_name)

            res = ProcessPool_run_renum_PDB("PDB_assembly", target_files_list_PDB, default_input_path_to_PDB_assembly, default_input_path_to_SIFTS,
                                            default_output_path_to_PDB_assembly, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.PDB_format_only:
            urls_to_target_PDB_files = url_formation_for_pool("PDB", parsed_input_text, default_input_path_to_PDB=default_input_path_to_PDB)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("PDB", urls_to_target_PDB_files, default_input_path_to_PDB=default_input_path_to_PDB)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            input_PDB_files_were_found = look_what_is_inside("PDB", default_input_path_to_PDB=default_input_path_to_PDB)
            passed_as_arg_file_4Char_PDB = list()
            for file_name in parsed_input_text:
                passed_as_arg_file_4Char_PDB.append(file_name[:4])
            target_files_list_PDB = list()
            for file_name in input_PDB_files_were_found:
                if file_name[3:7] in passed_as_arg_file_4Char_PDB:
                    target_files_list_PDB.append(file_name)

            res = ProcessPool_run_renum_PDB("PDB", target_files_list_PDB, default_input_path_to_PDB, default_input_path_to_SIFTS,
                                            default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.mmCIF_assembly_format_only:
            urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF_assembly", parsed_input_text,
                                                                default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("mmCIF_assembly", urls_to_target_mmCIF_files,
                                          default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            input_mmCIF_files_were_found = look_what_is_inside("mmCIF_assembly",
                                                               default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            passed_as_args_files_list_4Char = list()
            for file_name in parsed_input_text:
                passed_as_args_files_list_4Char.append(file_name[:4])

            target_files_list = list()
            for file_name in input_mmCIF_files_were_found:
                if file_name[:4] in passed_as_args_files_list_4Char:
                    target_files_list.append(file_name)

            if not os.path.exists(default_output_path_to_mmCIF_assembly):
                os.makedirs(default_output_path_to_mmCIF_assembly)

            res = ProcessPool_run_renum_mmCIF("mmCIF_assembly", target_files_list, default_input_path_to_mmCIF_assembly,
                                              default_input_path_to_SIFTS, default_output_path_to_mmCIF_assembly, default_mmCIF_num,
                                              gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        else:
            urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF", parsed_input_text, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("mmCIF", urls_to_target_mmCIF_files, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            input_mmCIF_files_were_found = look_what_is_inside("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            passed_as_args_files_list_4Char = list()
            for file_name in parsed_input_text:
                passed_as_args_files_list_4Char.append(file_name[:4])

            target_files_list = list()
            for file_name in input_mmCIF_files_were_found:
                if file_name[:4] in passed_as_args_files_list_4Char:
                    target_files_list.append(file_name)

            if not os.path.exists(default_output_path_to_mmCIF):
                os.makedirs(default_output_path_to_mmCIF)
            res = ProcessPool_run_renum_mmCIF("mmCIF", target_files_list, default_input_path_to_mmCIF, default_input_path_to_SIFTS,
                                              default_output_path_to_mmCIF, default_mmCIF_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

    # DOWNLOAD
    # DOWNLOAD FROM TEXT FILE or DOWNLOAD FROM LIST
    if args.download_from_text_file or args.download_from_list_of_arguments:
        if args.download_from_text_file:
            parsed_input_text = (input_text_file_parser(args.download_from_text_file))
        else:
            parsed_input_text = args.download_from_list_of_arguments

        if args.all_formats:
            urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF", parsed_input_text, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            urls_to_target_PDB_files = url_formation_for_pool("PDB", parsed_input_text, default_input_path_to_PDB=default_input_path_to_PDB)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("mmCIF", urls_to_target_mmCIF_files, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            run_downloads_with_ThreadPool("PDB", urls_to_target_PDB_files, default_input_path_to_PDB=default_input_path_to_PDB)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("mmCIF_assembly", urls_to_target_mmCIF_files,
                                          default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            # PDB_assembly download
            urls_to_target_PDB_files = list()
            all_urls_set = download_pdb_assemblies_list_with_lxml()
            for url in all_urls_set:
                for PDB_id in parsed_input_text:
                    if PDB_id in url:
                        urls_to_target_PDB_files.append(url)
            run_downloads_with_ThreadPool("PDB_assembly", urls_to_target_PDB_files,
                                          default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)

        elif args.PDB_assembly_format_only:
            urls_to_target_PDB_files = list()
            all_urls_set = download_pdb_assemblies_list_with_lxml()
            for url in all_urls_set:
                for PDB_id in parsed_input_text:
                    if PDB_id in url:
                        urls_to_target_PDB_files.append(url)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("PDB_assembly", urls_to_target_PDB_files,
                                          default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        elif args.PDB_format_only:
            urls_to_target_PDB_files = url_formation_for_pool("PDB", parsed_input_text, default_input_path_to_PDB=default_input_path_to_PDB)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("PDB", urls_to_target_PDB_files, default_input_path_to_PDB=default_input_path_to_PDB)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        elif args.mmCIF_assembly_format_only:
            urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF_assembly", parsed_input_text,
                                                                default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("mmCIF_assembly", urls_to_target_mmCIF_files,
                                          default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        else:
            urls_to_target_mmCIF_files = url_formation_for_pool("mmCIF", parsed_input_text, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            urls_to_target_SIFTS_files = url_formation_for_pool("SIFTS", parsed_input_text, default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            run_downloads_with_ThreadPool("mmCIF", urls_to_target_mmCIF_files, default_input_path_to_mmCIF=default_input_path_to_mmCIF)
            run_downloads_with_ThreadPool("SIFTS", urls_to_target_SIFTS_files, default_input_path_to_SIFTS=default_input_path_to_SIFTS)

    #################################################################################################################################################
    # WHOLE DB WORK #
    #################################################################################################################################################

    # RENUMBER ENTIRE DB
    if args.renumber_entire_database:
        if args.all_formats:
            print("Starting to renumber all databases...")
            print("Please, be patient...")
            # renumber mmCIF
            supreme_download_master("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            mmCIF_files_left_to_renumber = left_to_renumber_mmCIF(default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                                                  default_output_path_to_mmCIF=default_output_path_to_mmCIF)
            res = ProcessPool_run_renum_mmCIF("mmCIF", mmCIF_files_left_to_renumber, default_input_path_to_mmCIF, default_input_path_to_SIFTS,
                                              default_output_path_to_mmCIF, default_mmCIF_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

            # renumber mmCIF_assembly
            supreme_download_master("mmCIF_assembly", default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            mmCIF_assembly_left_to_renumber = left_to_renumber_mmCIF(default_input_path_to_mmCIF=default_input_path_to_mmCIF_assembly,
                                                                     default_output_path_to_mmCIF=default_output_path_to_mmCIF_assembly)
            ProcessPool_run_renum_mmCIF("mmCIF_assembly", mmCIF_assembly_left_to_renumber, default_input_path_to_mmCIF_assembly,
                                        default_input_path_to_SIFTS, default_output_path_to_mmCIF_assembly, default_mmCIF_num,
                                        gzip_mode, exception_AccessionIDs, nproc)

            # renumber PDB_assembly
            supreme_download_master("PDB_assembly", default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            PDB_assembly_files_left_to_renumber = left_to_renumber_mmCIF(default_input_path_to_mmCIF=default_input_path_to_PDB_assembly,
                                                                         default_output_path_to_mmCIF=default_output_path_to_PDB_assembly)
            ProcessPool_run_renum_PDB("PDB_assembly", PDB_assembly_files_left_to_renumber, default_input_path_to_PDB_assembly,
                                      default_input_path_to_SIFTS, default_output_path_to_PDB_assembly, default_PDB_num,
                                      gzip_mode, exception_AccessionIDs, nproc)

            # renumber PDB
            supreme_download_master("PDB", default_input_path_to_PDB=default_input_path_to_PDB,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            PDB_files_left_to_renumber = left_to_renumber_PDB(default_input_path_to_PDB=default_input_path_to_PDB,
                                                              default_output_path_to_PDB=default_output_path_to_PDB)
            ProcessPool_run_renum_PDB("PDB", PDB_files_left_to_renumber, default_input_path_to_PDB, default_input_path_to_SIFTS,
                                      default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)

        elif args.PDB_assembly_format_only:
            print("Starting to renumber entire PDB_assembly database...")
            print("Please, be patient...")
            supreme_download_master("PDB_assembly", default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            PDB_assembly_files_left_to_renumber = left_to_renumber_mmCIF(default_input_path_to_mmCIF=default_input_path_to_PDB_assembly,
                                                                         default_output_path_to_mmCIF=default_output_path_to_PDB_assembly)
            res = ProcessPool_run_renum_PDB("PDB_assembly", PDB_assembly_files_left_to_renumber, default_input_path_to_PDB_assembly,
                                            default_input_path_to_SIFTS, default_output_path_to_PDB_assembly, default_PDB_num,
                                            gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.PDB_format_only:
            print("Starting to renumber entire PDB database...")
            print("Please, be patient...")
            supreme_download_master("PDB", default_input_path_to_PDB=default_input_path_to_PDB,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            PDB_files_left_to_renumber = left_to_renumber_PDB(default_input_path_to_PDB=default_input_path_to_PDB,
                                                              default_output_path_to_PDB=default_output_path_to_PDB)
            res = ProcessPool_run_renum_PDB("PDB", PDB_files_left_to_renumber, default_input_path_to_PDB, default_input_path_to_SIFTS,
                                            default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.mmCIF_assembly_format_only:
            print("Starting to renumber entire mmCIF_assembly database...")
            print("Please, be patient...")
            supreme_download_master("mmCIF_assembly", default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            mmCIF_assembly_left_to_renumber = left_to_renumber_mmCIF(default_input_path_to_mmCIF=default_input_path_to_mmCIF_assembly,
                                                                     default_output_path_to_mmCIF=default_output_path_to_mmCIF_assembly)
            res = ProcessPool_run_renum_mmCIF("mmCIF_assembly", mmCIF_assembly_left_to_renumber, default_input_path_to_mmCIF_assembly,
                                              default_input_path_to_SIFTS, default_output_path_to_mmCIF_assembly, default_mmCIF_num,
                                              gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        else:
            print("Starting to renumber entire mmCIF database...")
            print("Please, be patient...")
            supreme_download_master("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            mmCIF_files_left_to_renumber = left_to_renumber_mmCIF(default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                                                  default_output_path_to_mmCIF=default_output_path_to_mmCIF)
            res = ProcessPool_run_renum_mmCIF("mmCIF", mmCIF_files_left_to_renumber, default_input_path_to_mmCIF, default_input_path_to_SIFTS,
                                              default_output_path_to_mmCIF, default_mmCIF_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

    # DOWNLOAD ENTIRE DB
    if args.download_entire_database:
        if args.all_formats:
            print("Starting to download all databases...")
            print("Please, be patient...")
            supreme_download_master("all", default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                    default_input_path_to_PDB=default_input_path_to_PDB,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS,
                                    default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly,
                                    default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly)

        elif args.PDB_format_only:
            print("Starting to download entire PDB database...")
            print("Please, be patient...")
            supreme_download_master("PDB", default_input_path_to_PDB=default_input_path_to_PDB,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        elif args.PDB_assembly_format_only:
            print("Starting to download entire PDB_assembly database...")
            print("Please, be patient...")
            supreme_download_master("PDB_assembly", default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        elif args.mmCIF_assembly_format_only:
            print("Starting to download entire mmCIF_assembly database...")
            print("Please, be patient...")
            supreme_download_master("mmCIF_assembly", default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)

        else:
            print("Starting to download entire mmCIF database...")
            print("Please, be patient...")
            supreme_download_master("mmCIF", default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                    default_input_path_to_SIFTS=default_input_path_to_SIFTS)

    # REFRESH ENTIRE DB
    if args.refresh_entire_database:
        if args.all_formats:
            print("Starting to refresh all databases...")
            print("Please, be patient...")
            left_to_refresh = supreme_download_master("all", "refresh", default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                                      default_input_path_to_PDB=default_input_path_to_PDB,
                                                      default_input_path_to_SIFTS=default_input_path_to_SIFTS,
                                                      default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly,
                                                      default_input_path_to_PDB_assembly=default_input_path_to_mmCIF_assembly)
            left_to_refresh_mmCIF = left_to_refresh[0]
            left_to_refresh_PDB = left_to_refresh[1]
            lefttodownload_mmCIF_assemblies = left_to_refresh[2]
            lefttodownload_PDB_assemblies = left_to_refresh[3]
            ProcessPool_run_renum_PDB("PDB_assembly", lefttodownload_PDB_assemblies, default_input_path_to_PDB_assembly,
                                      default_input_path_to_SIFTS, default_output_path_to_PDB_assembly,
                                      default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)
            ProcessPool_run_renum_PDB("PDB", left_to_refresh_PDB, default_input_path_to_PDB, default_input_path_to_SIFTS, default_output_path_to_PDB,
                                      default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)

            ProcessPool_run_renum_mmCIF("mmCIF_assembly", lefttodownload_mmCIF_assemblies, default_input_path_to_mmCIF_assembly,
                                        default_input_path_to_SIFTS, default_output_path_to_mmCIF_assembly, default_mmCIF_num,
                                        gzip_mode, exception_AccessionIDs, nproc)

            res = ProcessPool_run_renum_mmCIF("mmCIF", left_to_refresh_mmCIF, default_input_path_to_mmCIF, default_input_path_to_SIFTS,
                                              default_output_path_to_mmCIF, default_mmCIF_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.PDB_assembly_format_only:
            print("Starting to refresh entire PDB_assembly database...")
            print("Please, be patient...")
            left_to_refresh_PDB_assembly = supreme_download_master("PDB_assembly", "refresh",
                                                                   default_input_path_to_PDB_assembly=default_input_path_to_PDB_assembly,
                                                                   default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            res = ProcessPool_run_renum_PDB("PDB_assembly", left_to_refresh_PDB_assembly, default_input_path_to_PDB_assembly,
                                            default_input_path_to_SIFTS, default_output_path_to_PDB_assembly,
                                            default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.PDB_format_only:
            print("Starting to refresh entire PDB database...")
            print("Please, be patient...")
            left_to_refresh_PDB = supreme_download_master("PDB", "refresh", default_input_path_to_PDB=default_input_path_to_PDB,
                                                          default_input_path_to_SIFTS=default_input_path_to_SIFTS)
            res = ProcessPool_run_renum_PDB("PDB", left_to_refresh_PDB, default_input_path_to_PDB, default_input_path_to_SIFTS,
                                            default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        elif args.mmCIF_assembly_format_only:
            print("Starting to refresh entire mmCIF_assembly database...")
            print("Please, be patient...")
            left_to_refresh_mmCIF_assembly = supreme_download_master("mmCIF_assembly", "refresh",
                                                                     default_input_path_to_mmCIF_assembly=default_input_path_to_mmCIF_assembly,
                                                                     default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            res = ProcessPool_run_renum_mmCIF("mmCIF_assembly", left_to_refresh_mmCIF_assembly, default_input_path_to_mmCIF_assembly,
                                              default_input_path_to_SIFTS, default_output_path_to_mmCIF_assembly, default_mmCIF_num,
                                              gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)

        else:
            print("Starting to refresh entire mmCIF database...")
            print("Please, be patient...")
            left_to_refresh_mmCIF = supreme_download_master("mmCIF", "refresh", default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                                            default_input_path_to_SIFTS=default_input_path_to_SIFTS)

            res = ProcessPool_run_renum_mmCIF("mmCIF", left_to_refresh_mmCIF, default_input_path_to_mmCIF, default_input_path_to_SIFTS,
                                              default_output_path_to_mmCIF, default_mmCIF_num, gzip_mode, exception_AccessionIDs, nproc)
            log_writer(res)
