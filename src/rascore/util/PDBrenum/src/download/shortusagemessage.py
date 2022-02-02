def short_usage_messenger():
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

-redb, --renumber_entire_database
This option will download and renumber entire PDB database in PDB or/and mmCIF format
usage $ python3 PDB.py -redb -mmCIF
usage $ python3 PDB.py -redb -PDB
usage $ python3 PDB.py -redb -mmCIF_assembly
usage $ python3 PDB.py -redb -PDB_assembly
usage $ python3 PDB.py -redb -all

for more info --help

Roland Dunbrack's Lab
Fox Chase Cancer Center
Philadelphia, PA
2020
        """
