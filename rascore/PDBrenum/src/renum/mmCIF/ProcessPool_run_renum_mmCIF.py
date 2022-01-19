from src.download.modules import *
from src.download.lookfilesinside import look_what_is_inside
from src.renum.mmCIF.new_mmCIF import master_mmCIF_renumber_function


REMARK_mmCIF = ["#\n",
                "loop_\n",
                "_database_PDB_remark.id       1\n",
                "_database_PDB_remark.text\n",
                ";File processed by PDBrenum: http://dunbrack3.fccc.edu/PDBrenum\n",
                "Author sequence numbering is replaced with UniProt numbering according to\n",
                "alignment by SIFTS (https://www.ebi.ac.uk/pdbe/docs/sifts/).\n",
                "Only chains with UniProt sequences in SIFTS are renumbered.\n",
                "Residues in UniProt chains without UniProt residue numbers in SIFTS\n",
                "(e.g., sequence tags) are given residue numbers 50000+label_seq_id\n",
                "(where label_seq_id is the 1-to-N residue numbering of each chain.\n",
                "Ligands are numbered 50000+their residue number in the original file.\n",
                "The _poly_seq_scheme table contains a correspondence between the\n",
                "1-to-N sequence (seq_id), the new numbering based on UniProt (pdb_seq_num =\n",
                "auth_seq_id in the _atom_site records), and the author numbering\n",
                "in the original mmCIF file from the PDB (auth_seq_num).\n",
                ";\n",
                "#\n"]


def check_assemblies(mmCIF_assembly, default_output_path_to_mmCIF_assembly):
    output_mmCIF_assembly_files_were_found_list = list()
    output_mmCIF_assembly_files_were_found_list.append(mmCIF_assembly)
    new_order_with_remark = 0

    for name in output_mmCIF_assembly_files_were_found_list:
        not_gzip = 1
        try:
            list_of_lines_from_assembly_file = gzip.open(
                Path(str(default_output_path_to_mmCIF_assembly) + "/" + name), 'rt').readlines()
        except OSError:
            # maybe not archived
            try:
                list_of_lines_from_assembly_file = open(
                    Path(str(default_output_path_to_mmCIF_assembly) + "/" + name), 'rt').readlines()
                not_gzip = 0
            except Exception:
                # broken archive
                os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))
                continue
        except Exception:
            # broken archive
            os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))
            continue

        # check if file startswith "_atom_site" table at the beginning
        try:
            if "_atom_site" in list_of_lines_from_assembly_file[3] and "loop_" in list_of_lines_from_assembly_file[2]:
                pass
            else:
                if "renum" in name:
                    for line in list_of_lines_from_assembly_file:
                        if line.startswith("_entry.id"):
                            new_order_with_remark = (list_of_lines_from_assembly_file[:(list_of_lines_from_assembly_file.index(line)) + 1]
                                                     + REMARK_mmCIF
                                                     + list_of_lines_from_assembly_file[(list_of_lines_from_assembly_file.index(line)) + 2:])
                    if new_order_with_remark == 0:
                        new_order_with_remark = list_of_lines_from_assembly_file
                        print(name)

                    if not_gzip != 0:
                        with gzip.open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name), "wt") as gzip_out:
                            for listitem in new_order_with_remark:
                                gzip_out.write(listitem)

                    else:
                        with open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name), "wt") as file_out:
                            for listitem in new_order_with_remark:
                                file_out.write(listitem)

                else:
                    new_order_with_remark = list_of_lines_from_assembly_file
                    if not_gzip != 0:
                        with gzip.open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name.split(".")[0] + "_renum.cif.gz"),
                                       "wt") as gzip_out:
                            for listitem in new_order_with_remark:
                                gzip_out.write(listitem)
                        os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))
                    else:
                        with open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name.split(".")[0] + "_renum.cif"), "wt") as file_out:
                            for listitem in new_order_with_remark:
                                file_out.write(listitem)
                        os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))

                return name

        except IndexError:
            # empty file
            os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))
            continue

        try:
            new_order_for_assembly_file = (list_of_lines_from_assembly_file[:1]
                                           + list_of_lines_from_assembly_file[list_of_lines_from_assembly_file.index("#\n", 2):]
                                           + list_of_lines_from_assembly_file[2:list_of_lines_from_assembly_file.index("#\n", 2)]
                                           + ["#\n"])

            if "renum" in name:
                for line in new_order_for_assembly_file:
                    if line.startswith("_entry.id"):
                        new_order_with_remark = (new_order_for_assembly_file[:(new_order_for_assembly_file.index(line)) + 1]
                                                 + REMARK_mmCIF
                                                 + new_order_for_assembly_file[(new_order_for_assembly_file.index(line)) + 2:])

                if new_order_with_remark != 0:
                    if not_gzip != 0:
                        with gzip.open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name), "wt") as gzip_out:
                            for listitem in new_order_with_remark:
                                gzip_out.write(listitem)

                    else:
                        with open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name), "wt") as file_out:
                            for listitem in new_order_with_remark:
                                file_out.write(listitem)
            else:
                new_order_with_remark = new_order_for_assembly_file
                if not_gzip != 0:
                    with gzip.open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name.split(".")[0] + "_renum.cif.gz"), "wt") as gzip_out:
                        for listitem in new_order_with_remark:
                            gzip_out.write(listitem)
                    os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))
                else:
                    with open(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name.split(".")[0] + "_renum.cif"), "wt") as file_out:
                        for listitem in new_order_with_remark:
                            file_out.write(listitem)
                    os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))
            return name

        except ValueError:
            # file isn't complete
            os.remove(Path(str(default_output_path_to_mmCIF_assembly) + "/" + name))


def ProcessPool_run_renum_mmCIF(format_mmCIF, mmCIF_to_renumber, default_input_path_to_mmCIF,
                                default_input_path_to_SIFTS, default_output_path_to_mmCIF, default_mmCIF_num,
                                gzip_mode, exception_AccessionIDs, nproc):
    first_res = 0

    for i in range(3):
        if not os.path.exists(default_output_path_to_mmCIF):
            os.makedirs(default_output_path_to_mmCIF)

        # renumber loop
        resulting = list()
        executor = ProcessPoolExecutor(max_workers=nproc)
        partial_master_mmCIF_renumber_function = partial(master_mmCIF_renumber_function,
                                                         default_input_path_to_mmCIF=default_input_path_to_mmCIF,
                                                         default_input_path_to_SIFTS=default_input_path_to_SIFTS,
                                                         default_output_path_to_mmCIF=default_output_path_to_mmCIF,
                                                         default_mmCIF_num=default_mmCIF_num, gzip_mode=gzip_mode,
                                                         exception_AccessionIDs=exception_AccessionIDs)
        jobs = [executor.submit(partial_master_mmCIF_renumber_function, mmCIF_files) for mmCIF_files in mmCIF_to_renumber]
        for job in tqdm.tqdm(as_completed(jobs), total=len(jobs), miniters=1, position=0,
                             leave=True, desc="Renumbering " + format_mmCIF + " files"):
            result = job.result()
            if result is None:
                continue
            resulting.append(result)

        if i == 0:
            first_res = resulting

        if format_mmCIF == "mmCIF_assembly":
            output_mmCIF = look_what_is_inside('output_mmCIF_assembly', default_output_path_to_mmCIF_assembly=default_output_path_to_mmCIF)
        else:
            output_mmCIF = look_what_is_inside('output_mmCIF', default_output_path_to_mmCIF=default_output_path_to_mmCIF)

        # checker loop
        check_list = list()
        executor = ProcessPoolExecutor(max_workers=nproc)
        partial_reform_assembly = partial(check_assemblies, default_output_path_to_mmCIF_assembly=default_output_path_to_mmCIF)
        jobs = [executor.submit(partial_reform_assembly, assembly_files) for assembly_files in output_mmCIF]
        for job in tqdm.tqdm(as_completed(jobs), total=len(jobs), miniters=1, position=0,
                             leave=True, desc="Checking " + format_mmCIF + " files"):
            resultus = job.result()
            check_list.append(resultus)

        if format_mmCIF == "mmCIF_assembly":
            output_mmCIF = look_what_is_inside('output_mmCIF_assembly', default_output_path_to_mmCIF_assembly=default_output_path_to_mmCIF)
        else:
            output_mmCIF = look_what_is_inside('output_mmCIF', default_output_path_to_mmCIF=default_output_path_to_mmCIF)

        output_mmCIF_4char = set()
        for n in output_mmCIF:
            output_mmCIF_4char.add(n[:4])

        if len(check_list) <= len(output_mmCIF):
            break
        else:
            new_round_mmCIF_to_renumber = set()
            for n in mmCIF_to_renumber:
                if n[:4] in output_mmCIF_4char:
                    continue
                else:
                    new_round_mmCIF_to_renumber.add(n)
            mmCIF_to_renumber = new_round_mmCIF_to_renumber

    return first_res
