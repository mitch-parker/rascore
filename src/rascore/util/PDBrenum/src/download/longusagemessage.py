def long_usage_messenger():
    return """
PDB.py
optional arguments:
-h, --help            show this help message and exit

-rftf text_file_with_PDB.tsv, --renumber_from_text_file text_file_with_PDB.tsv
This option will download and renumber specified files
usage $ python3 PDB.py -rftf text_file_with_PDB_in_it.tsv -mmCIF
usage $ python3 PDB.py -rftf text_file_with_PDB_in_it.tsv -PDB
usage $ python3 PDB.py -rftf text_file_with_PDB_in_it.tsv -mmCIF_assembly
usage $ python3 PDB.py -rftf text_file_with_PDB_in_it.tsv -PDB_assembly
usage $ python3 PDB.py -rftf text_file_with_PDB_in_it.tsv -all

-rfla [6dbp 3v03 2jit ...], --renumber_from_list_of_arguments [6dbp 3v03 2jit ...]
This option will download and renumber specified files
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -mmCIF
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -PDB
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -mmCIF_assembly
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -PDB_assembly
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -all

-dftf text_file_with_PDB.tsv, --download_from_text_file text_file_with_PDB.tsv
This option will read given input file parse by space
or tab or comma or new line and download it example 
usage $ python3 PDB.py -dftf text_file_with_PDB_in_it.tsv -mmCIF
usage $ python3 PDB.py -dftf text_file_with_PDB_in_it.tsv -PDB
usage $ python3 PDB.py -dftf text_file_with_PDB_in_it.tsv -mmCIF_assembly
usage $ python3 PDB.py -dftf text_file_with_PDB_in_it.tsv -PDB_assembly
usage $ python3 PDB.py -dftf text_file_with_PDB_in_it.tsv -all

-dfla [6dbp 3v03 2jit ...], --download_from_list_of_arguments 6dbp 3v03 2jit ...]
This option will read given list of arguments separated by space. 
Format of the list should be without any commas or quotation marks
usage $ python3 PDB.py -dfla 6dbp 3v03 2jit -mmCIF
usage $ python3 PDB.py -dfla 6dbp 3v03 2jit -PDB
usage $ python3 PDB.py -dfla 6dbp 3v03 2jit -mmCIF_assembly
usage $ python3 PDB.py -dfla 6dbp 3v03 2jit -PDB_assembly
usage $ python3 PDB.py -dfla 6dbp 3v03 2jit -all

-redb, --renumber_entire_database
This option will download and renumber entire PDB database in PDB or/and mmCIF format
usage $ python3 PDB.py -redb -mmCIF
usage $ python3 PDB.py -redb -PDB
usage $ python3 PDB.py -redb -mmCIF_assembly
usage $ python3 PDB.py -redb -PDB_assembly
usage $ python3 PDB.py -redb -all 

-dall, --download_entire_database
This option will download entire mmCIF database
usage $ python3 PDB.py -dall -mmCIF
usage $ python3 PDB.py -dall -PDB
usage $ python3 PDB.py -dall -mmCIF_assembly
usage $ python3 PDB.py -dall -PDB_assembly
usage $ python3 PDB.py -dall -all

-refr, --refresh_entire_database
This option will delete outdated files and download
fresh ones. This option makes sense and only works if
you work with entire database
usage $ python3 PDB.py -refr -mmCIF
usage $ python3 PDB.py -refr -PDB
usage $ python3 PDB.py -refr -mmCIF_assembly
usage $ python3 PDB.py -refr -PDB_assembly
usage $ python3 PDB.py -refr -all    

-PDB, --PDB_format_only
This option will specify working format to pdb format
-mmCIF, --mmCIF_format_only
This option will specify working format to mmCIF format (default)
-PDB_assembly, --PDB_assembly_format_only
This option will specify working format to pdb format
-mmCIF_assembly, --mmCIF_assembly_format_only
This option will specify working format to mmCIF format
-all, --all_formats   This option will work with both formats

argpar.add_argument("-sipm", "--set_default_input_path_to_mmCIF", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sipma", "--set_default_input_path_to_mmCIF_assembly", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sipp", "--set_default_input_path_to_PDB", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sippa", "--set_default_input_path_to_PDB_assembly", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sips", "--set_default_input_path_to_SIFTS", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sopm", "--set_default_output_path_to_mmCIF", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sopma", "--set_default_output_path_to_mmCIF_assembly", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-sopp", "--set_default_output_path_to_PDB", type=str, help=argparse.SUPPRESS)
argpar.add_argument("-soppa", "--set_default_output_path_to_PDB_assembly", type=str, help=argparse.SUPPRESS)

-sipm, --set_default_input_path_to_mmCIF
This option will set default input path to mmCIF files (default: <./mmCIF>)
usage $ python3 PDB.py -sipm /Users/bulatfaezov/PycharmProjects/renum/venv/mmCIF
-sipp, --set_default_input_path_to_PDB
This option will set default input path to PDB files (default: <./PDB>)
usage $ python3 PDB.py -sipp /Users/bulatfaezov/PycharmProjects/renum/venv/PDB
-sipma, --set_default_input_path_to_mmCIF_assembly
This option will set default input path to mmCIF_assembly files (default: <./mmCIF_assembly>)
usage $ python3 PDB.py -sipm /Users/bulatfaezov/PycharmProjects/renum/venv/mmCIF_assembly
-sippa, --set_default_input_path_to_PDB_assembly
This option will set default input path to PDB_assembly files (default: <./PDB_assembly>)
usage $ python3 PDB.py -sipp /Users/bulatfaezov/PycharmProjects/renum/venv/PDB_assembly
-sips, --set_default_input_path_to_SIFTS
This option will set default input path to SIFTS files (default: <./SIFTS>)
usage $ python3 PDB.py -sips /Users/bulatfaezov/PycharmProjects/renum/venv/SIFTS
-sopm, --set_default_output_path_to_mmCIF
This option will set default output path to mmCIF files (default: <./output_mmCIF>)
usage $ python3 PDB.py -sopm /Users/bulatfaezov/PycharmProjects/renum/venv/output_mmCIF
-sopp, --set_default_output_path_to_PDB
This option will set default output path to PDB files (default: <./output_PDB>)
usage $ python3 PDB.py -sopp /Users/bulatfaezov/PycharmProjects/renum/venv/output_PDB
-sopma, --set_default_output_path_to_mmCIF_assembly
This option will set default output path to mmCIF_assembly files (default: <./output_mmCIF_assembly>)
usage $ python3 PDB.py -sopm /Users/bulatfaezov/PycharmProjects/renum/venv/output_mmCIF_assembly
-soppa, --set_default_output_path_to_PDB_assembly
This option will set default output path to PDB_assembly files (default: <./output_PDB_assembly>)
usage $ python3 PDB.py -sopp /Users/bulatfaezov/PycharmProjects/renum/venv/output_PDB_assembly

-sdmn, --set_default_mmCIF_num
This option will set default mmCIF number which will be added to 1 to end numbering in cases 
when there are no UniProt numbering (default: 50000)
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -mmCIF -sdmn 50000

-sdpn, --set_default_PDB_num
This option will set default PDB number which will be added to 1 to end numbering in cases 
when there are no UniProt numbering (default: 5000)
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -mmCIF -sdpn 5000

"-offz", "--set_to_off_mode_gzip"
By default program will compress files with gzip this option will turn that off
(default: gzip is on)
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -mmCIF -offz

"-nproc", "--set_number_of_processes"
By default program will use all available CPUs. User can reduce number of CPUs for PDBrenum.
In this example: only 4 CPUs will be used by the PDBrenum even if more CPUs available
(default: nproc = None)
usage $ python3 PDB.py -rfla 6dbp 3v03 2jit -mmCIF -nproc 4


Roland Dunbrack's Lab
Fox Chase Cancer Center
Philadelphia, PA
2020
        """
