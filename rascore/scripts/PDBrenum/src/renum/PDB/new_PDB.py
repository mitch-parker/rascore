from src.download.modules import *
from src.renum.shared.handling_chain_numbering_clashes import handling_chain_numbering_clashes
from src.renum.shared.SIFTS_tree_parser import SIFTS_tree_parser
from src.renum.shared.renumbered_count_in_chains import renumbered_count_in_chains
from src.download.downloadwithThreadPool import download_with_pool, url_formation_for_pool
PDBrenum_REMARK = [
    "REMARK   0  File processed by PDBrenum: http://dunbrack3.fccc.edu/PDBrenum      ",
    "REMARK   0  Author sequence numbering is replaced with UniProt numbering        ",
    "REMARK   0  according to alignment by SIFTS                                     ",
    "REMARK   0  (https://www.ebi.ac.uk/pdbe/docs/sifts/).                           ",
    "REMARK   0  Only chains with UniProt sequences in SIFTS are renumbered.         ",
    "REMARK   0  Residues in UniProt chains without UniProt residue numbers in SIFTS ",
    "REMARK   0  (e.g., sequence tags) are given residue numbers 5000+label_seq_id   ",
    "REMARK   0  (where label_seq_id is the 1-to-N residue numbering of each chain.  ",
    "REMARK   0  Ligands are numbered 5000+their residue number in the original      ",
    "REMARK   0  file. The _poly_seq_scheme table contains a correspondence between  ",
    "REMARK   0  the 1-to-N sequence (seq_id), the new numbering based on UniProt    ",
    "REMARK   0  (pdb_seq_num = auth_seq_id in the _atom_site records), and          ",
    "REMARK   0  the author numbering in the original mmCIF file from the PDB        ",
    "REMARK   0  (auth_seq_num).                                                     "]


def SIFTS_data_parser_for_PDB(tuple_PDBe_for_PDB_and_tuple_PDB, tuple_PDBe_for_UniProt_and_tuple_UniProt,
                              default_PDB_num, chains_to_change="all"):
    df_PDBe_UniProt = pd.DataFrame(tuple_PDBe_for_UniProt_and_tuple_UniProt, columns=['PDBe', 'UniProt', "AccessionID"])
    df_PDBe_UniProt = df_PDBe_UniProt.drop_duplicates(subset="PDBe", keep='first')
    df_PDBe_PDB = pd.DataFrame(tuple_PDBe_for_PDB_and_tuple_PDB, columns=['PDBe', 'PDB'])
    df_PDBe_PDB = df_PDBe_PDB.drop_duplicates(subset="PDBe", keep='first')

    df_PDBe_PDB_UniProt = df_PDBe_PDB.merge(df_PDBe_UniProt, left_on="PDBe", right_on="PDBe", how='left')
    df_PDBe_PDB_UniProt['UniProt'] = df_PDBe_PDB_UniProt['UniProt'].replace(np.nan, "5000")
    df_PDBe_PDB_UniProt["Uni_moD"] = np.where(df_PDBe_PDB_UniProt['UniProt'] != "5000", df_PDBe_PDB_UniProt['UniProt'], df_PDBe_PDB_UniProt["PDBe"])
    df_PDBe_PDB_UniProt.loc[:, 'new_col_Uni'] = df_PDBe_PDB_UniProt.Uni_moD.map(lambda x: x[0])
    df_PDBe_PDB_UniProt["UniProt_5k"] = df_PDBe_PDB_UniProt.new_col_Uni.apply(lambda x: (int(x) + default_PDB_num if type(x) == str else x))
    df_PDBe_PDB_UniProt.loc[df_PDBe_PDB_UniProt['UniProt'] != '5000', 'UniProt_5k'] = df_PDBe_PDB_UniProt['new_col_Uni']

    Three_Rows_CIF_Num_Uni = []
    if chains_to_change == "all":
        for index, rows in df_PDBe_PDB_UniProt.iterrows():
            intermediate_list = [rows.PDBe, rows.UniProt_5k, rows.Uni_moD, rows.PDB, rows.AccessionID]
            Three_Rows_CIF_Num_Uni.append(intermediate_list)

    else:
        for index, rows in df_PDBe_PDB_UniProt.iterrows():
            if rows.PDB[2].strip() in chains_to_change:
                intermediate_list = [rows.PDBe, rows.UniProt_5k, rows.Uni_moD, rows.PDB, rows.AccessionID]
            else:
                intermediate_list = [rows.PDBe, rows.PDB[0], rows.Uni_moD, rows.PDB, rows.AccessionID]
            Three_Rows_CIF_Num_Uni.append(intermediate_list)

    df_PDBe_PDB_UniProt["Three_Rows_CIF_Num_Uni"] = Three_Rows_CIF_Num_Uni
    df_PDBe_PDB_UniProt_without_null = df_PDBe_PDB_UniProt[df_PDBe_PDB_UniProt.PDB.map(lambda x: x[0]) != "null"]
    df_PDBe_PDB_UniProt_without_null_index_PDBe = df_PDBe_PDB_UniProt_without_null.set_index("PDBe")

    return [df_PDBe_PDB_UniProt_without_null_index_PDBe, df_PDBe_PDB_UniProt]


def try_SIFTS_tree_parser(default_input_path_to_SIFTS, SIFTS_name):
    product_tree_SIFTS = 0
    for _ in range(3):
        try:
            product_tree_SIFTS = SIFTS_tree_parser(
                gzip.open(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name), 'rt'))
            break
        except EOFError:
            os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name))
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name], default_input_path_to_SIFTS=default_input_path_to_SIFTS)[0],
                               default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        except ValueError:
            os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name))
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name], default_input_path_to_SIFTS=default_input_path_to_SIFTS)[0],
                               default_input_path_to_SIFTS=default_input_path_to_SIFTS)
        except OSError:
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name], default_input_path_to_SIFTS=default_input_path_to_SIFTS)[0],
                               default_input_path_to_SIFTS=default_input_path_to_SIFTS)
    return product_tree_SIFTS


def try_PDB(default_input_path_to_PDB, PDB):
    split = 0

    for _ in range(3):
        try:
            split = gzip.open(Path(str(default_input_path_to_PDB) + "/" + PDB), 'rt').read().splitlines()
            break
        except EOFError:
            try:
                re.search('\.pdb(.*).gz', PDB).group(1)
                os.remove(Path(str(default_input_path_to_PDB) + "/" + PDB))
                download_with_pool(url_formation_for_pool("PDB_assembly", [PDB], default_input_path_to_PDB_assembly=default_input_path_to_PDB)[0],
                                   default_input_path_to_PDB_assembly=default_input_path_to_PDB)
            except AttributeError:
                os.remove(Path(str(default_input_path_to_PDB) + "/" + PDB))
                download_with_pool(url_formation_for_pool("PDB", [PDB], default_input_path_to_PDB=default_input_path_to_PDB)[0],
                                   default_input_path_to_PDB=default_input_path_to_PDB)

        except ValueError:
            try:
                re.search('\.pdb(.*).gz', PDB).group(1)
                os.remove(Path(str(default_input_path_to_PDB) + "/" + PDB))
                download_with_pool(url_formation_for_pool("PDB_assembly", [PDB], default_input_path_to_PDB_assembly=default_input_path_to_PDB)[0],
                                   default_input_path_to_PDB_assembly=default_input_path_to_PDB)
            except AttributeError:
                os.remove(Path(str(default_input_path_to_PDB) + "/" + PDB))
                download_with_pool(url_formation_for_pool("PDB", [PDB], default_input_path_to_PDB=default_input_path_to_PDB)[0],
                                   default_input_path_to_PDB=default_input_path_to_PDB)
        except OSError:
            try:
                re.search('\.pdb(.*).gz', PDB).group(1)
                download_with_pool(url_formation_for_pool("PDB_assembly", [PDB], default_input_path_to_PDB_assembly=default_input_path_to_PDB)[0],
                                   default_input_path_to_PDB_assembly=default_input_path_to_PDB)
            except AttributeError:
                download_with_pool(url_formation_for_pool("PDB", [PDB], default_input_path_to_PDB=default_input_path_to_PDB)[0],
                                   default_input_path_to_PDB=default_input_path_to_PDB)
    return split


def if_no_SIFTS_data_log_for_PDB(default_input_path_to_PDB, PDB_id, PDB):
    split = try_PDB(default_input_path_to_PDB, PDB)
    res_number_name_chainID_from_PDB_tuple = list()
    chains_set = set()
    log_message = list()

    for n in split:
        if n.startswith("ATOM") or n.startswith("TER") or n.startswith("ANISOU") or n.startswith("ANISOU") or n.startswith("SIGUIJ"):
            res_number_name_chainID_from_PDB_tuple.append((n[22:27].strip(" "), n[17:20], n[21]))
            chains_set.add(n[21])

    if len(res_number_name_chainID_from_PDB_tuple) == 0:
        log_message.append([PDB_id, "-", "-", "-", "-", "-", "-", "-", "-", "-"])
        return log_message

    df_mmCIF = pd.DataFrame(list(zip(res_number_name_chainID_from_PDB_tuple, res_number_name_chainID_from_PDB_tuple)))
    df_mmCIF = df_mmCIF.rename(columns={0: "PDB_old", 1: "PDB_old_copy"})
    df_mmCIF = df_mmCIF.set_index("PDB_old")
    df_mmCIF = df_mmCIF.drop_duplicates()

    for chain in chains_set:
        count_res_in_chain = 0
        for resnum_resname_chain in df_mmCIF.PDB_old_copy:
            if chain == resnum_resname_chain[2]:
                count_res_in_chain += 1
        log_message.append([PDB_id, "-", chain, "-", "-", len(df_mmCIF), "-", count_res_in_chain, "0", "0"])
    return log_message


def copy_file(inpath, file_name, outpath, postfix, gzip_mode):
    if file_name.endswith(".ent.gz") and file_name.startswith("pdb"):
        PDB_id = file_name[3:file_name.rfind(".ent.gz")]
    else:
        PDB_id = file_name[:4]

    absolute_path_in = inpath + "/" + file_name
    absolute_path_out = outpath + "/" + PDB_id + postfix
    if gzip_mode == "off":
        with gzip.open(absolute_path_in, 'rb') as f_in:
            with open(absolute_path_out, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(absolute_path_in, absolute_path_out + ".gz")


def PDB_parser(split, df_PDBe_PDB_UniProt_without_null_index_PDBe, default_PDB_num):
    res_number_name_chainID_from_PDB_tuple = list()
    num_ins_code_name_chain_HETATM = list()
    missing_res_remark_465 = list()
    skipper_for_remark_465 = True
    Num_is_Too_Big = True

    for n in split:
        if n.startswith("ATOM") or n.startswith("TER") or n.startswith("ANISOU") or n.startswith("ANISOU") or n.startswith("SIGUIJ"):
            res_number_name_chainID_from_PDB_tuple.append((n[22:27].strip(" "), n[17:20], n[21]))
        if n.startswith("HETATM"):
            num_ins_code_name_chain_HETATM.append((n[22:27].strip(" "), n[17:20], n[21]))

        if n.startswith("REMARK 465"):
            if not skipper_for_remark_465:
                missing_res_remark_465.append((n[20:27].strip(" "), n[15:18], n[19]))
            if n[15:27] == "RES C SSSEQI":
                skipper_for_remark_465 = False

    df_mmCIF = pd.DataFrame(list(zip(res_number_name_chainID_from_PDB_tuple, res_number_name_chainID_from_PDB_tuple)))
    df_mmCIF = df_mmCIF.rename(columns={0: "PDB_old", 1: "PDB_old_copy"})
    df_mmCIF = df_mmCIF.set_index("PDB_old")
    df_mmCIF = df_mmCIF.drop_duplicates()

    df_final = df_mmCIF.merge(df_PDBe_PDB_UniProt_without_null_index_PDBe, left_on="PDB_old_copy", right_on="PDB", how='left')
    df_final['Uni_moD'] = df_final['Uni_moD'].replace(np.nan, "5000")
    df_final["Uni_moD"] = np.where(df_final['Uni_moD'] != "5000", df_final['Uni_moD'], df_final["PDB_old_copy"])
    df_final.loc[:, 'new_col_Uni'] = df_final.Uni_moD.map(lambda x: x[0])
    df_final["UniProt_5k"] = df_final.new_col_Uni.apply(lambda x: (int(x) + default_PDB_num if x.isdigit() else x))
    df_final.loc[df_final['UniProt'] != '5000', 'UniProt_5k'] = df_final['new_col_Uni']

    Three_Rows_CIF_Num_Uni = list()
    for index, rows in df_final.iterrows():
        intermediate_list = [rows.PDB_old_copy, rows.UniProt_5k, rows.Uni_moD]
        if type(rows.UniProt_5k) == int:
            if len(str(rows.UniProt_5k)) > 4:
                Num_is_Too_Big = False
        else:
            if len(rows.UniProt_5k) > 4:
                Num_is_Too_Big = False
        Three_Rows_CIF_Num_Uni.append(intermediate_list)

    df_final["Three_Rows_CIF_Num_Uni"] = Three_Rows_CIF_Num_Uni
    df_final_index_PDBe = df_final.set_index("PDB_old_copy")
    df_final_index_PDBe_drop_NAN = df_final_index_PDBe.dropna(subset=['PDB'])
    pd_series_index_PDBe = df_final_index_PDBe["Three_Rows_CIF_Num_Uni"]

    PDB_str = df_final_index_PDBe_drop_NAN.PDB.map(
        lambda x: x[1] + "{:>2}".format(x[2]) + "{:>4}".format(x[0]) + " " if x[0].isdigit() else x[1] + "{:>2}".format(x[2]) + "{:>5}".format(x[0]))
    df_final_index_PDBe_drop_NAN = df_final_index_PDBe_drop_NAN.merge(PDB_str.rename('PDB_str'), left_index=True, right_index=True)
    df_final_poly_corrected = df_final_index_PDBe_drop_NAN.drop(columns=['UniProt', 'AccessionID', "new_col_Uni", "UniProt_5k", "Uni_moD"])
    renum_str = df_final_poly_corrected.PDB.map(lambda x: x[1] + "{:>2}".format(x[2])) + df_final_poly_corrected.Three_Rows_CIF_Num_Uni.map(
        lambda x: "{:>4}".format(str(int(x[1])))) + " "
    df_final_poly_corrected = df_final_poly_corrected.merge(renum_str.rename('renum_str'), left_index=True, right_index=True)
    df_final_poly_corrected = df_final_poly_corrected.reset_index(drop=True)

    return [pd_series_index_PDBe, num_ins_code_name_chain_HETATM, df_final_poly_corrected, missing_res_remark_465, Num_is_Too_Big]


def non_poly_num(pd_series_index_PDBe, num_ins_code_name_chain_HETATM):
    working_range_list = list()
    for n in range(1, 10000):
        working_range_list.append(n)

    chain_and_number = list()
    for n in pd_series_index_PDBe:
        chain_and_number.append((n[0][2], n[1]))

    chain_label = chain_and_number[0][0]
    numbers_per_chain = list()
    chain_label_with_numbers_per_chain = list()
    for n in chain_and_number:
        if chain_label == n[0]:
            d = str(n[1])
            n_numeric = ''.join(d for d in d if d.isdigit())
            numbers_per_chain.append(int(n_numeric))
        else:
            numbers_per_chain = list(np.unique(numbers_per_chain))
            chain_label_with_numbers_per_chain.append((chain_label, numbers_per_chain))
            numbers_per_chain = list()
            chain_label = n[0]
            if chain_label == n[0]:
                d = str(n[1])
                n_numeric = ''.join(d for d in d if d.isdigit())
                numbers_per_chain.append(int(n_numeric))

    numbers_per_chain = list(np.unique(numbers_per_chain))
    chain_label_with_numbers_per_chain.append((chain_label, numbers_per_chain))

    available_numbers_for_chains = list()
    for n in chain_label_with_numbers_per_chain:
        for num in n[1]:
            if num in working_range_list:
                working_range_list.remove(num)

        available_numbers_for_chains.append((n[0], working_range_list))
        working_range_list = list()
        for n_ in range(1, 10000):
            working_range_list.append(n_)

    chain_and_num_available = list()
    for n in available_numbers_for_chains:
        for num in n[1]:
            chain_and_num_available.append((n[0], num))

    numbers_from_num_ins_code_name_chain_HETATM = list()
    chain_from_num_ins_code_name_chain_HETATM = list()
    for n in num_ins_code_name_chain_HETATM:
        numbers_from_num_ins_code_name_chain_HETATM.append(n[0])
        chain_from_num_ins_code_name_chain_HETATM.append(n[2])

    df_nonpoly = pd.DataFrame(
        list(zip(num_ins_code_name_chain_HETATM, chain_from_num_ins_code_name_chain_HETATM, numbers_from_num_ins_code_name_chain_HETATM)), columns=[
            'PDB', "PDB_chain", "numbers"])
    df_nonpoly_dropped_dup = df_nonpoly.drop_duplicates(subset="PDB", keep='first')
    df_nonpoly_dropped_dup_sorted = df_nonpoly_dropped_dup.sort_values(["PDB_chain", "numbers"], ascending=(True, True)).reset_index(drop=True)
    small_ref_table = df_nonpoly_dropped_dup_sorted.set_index(["PDB_chain", "PDB"]).count(level="PDB_chain")

    all_nonpoly_chains = list()
    for n in small_ref_table.index:
        all_nonpoly_chains.append(n)

    checked_chains_list = list()
    for n in chain_and_num_available:
        checked_chains_list.append(n[0])
    checked_chains_list_uniq = list(np.unique(checked_chains_list))

    available_numbers_to_chains = list()
    for n in all_nonpoly_chains:
        if n not in checked_chains_list_uniq:
            for z in range(1, 10000):
                available_numbers_to_chains.append((n, z))

    chain_and_num_available.extend(available_numbers_to_chains)

    df_chain_and_num_available = pd.DataFrame(chain_and_num_available, columns=['available_chain', "available_number"])
    df_chain_and_num_available_sorted = df_chain_and_num_available.drop_duplicates(
        subset=["available_chain", "available_number"], keep="first").sort_values(
        ["available_chain", "available_number"], ascending=(True, False)).reset_index(drop=True)

    df_for_nonpoly_replace = pd.DataFrame(list(), columns=['available_chain', "available_number"])
    for n in df_nonpoly_dropped_dup_sorted.set_index(["PDB_chain", "PDB"]).count(level="PDB_chain").index:
        temporal_df_for_addition_of_available_num = df_chain_and_num_available_sorted.where(
            df_chain_and_num_available_sorted['available_chain'] == n).dropna()[0:(small_ref_table["numbers"][n])]
        df_for_nonpoly_replace = df_for_nonpoly_replace.append(temporal_df_for_addition_of_available_num, ignore_index=True)

    df_final_nonpoly = pd.merge(left=df_nonpoly_dropped_dup_sorted, right=df_for_nonpoly_replace, left_index=True, right_index=True)

    df_final_nonpoly.loc[:, 'PDB_str'] = df_final_nonpoly.PDB.map(
        lambda x: x[1] + "{:>2}".format(x[2]) + "{:>4}".format(x[0]) + " " if x[0].isdigit() else x[1] + "{:>2}".format(x[2]) + "{:>5}".format(x[0]))
    df_final_nonpoly_corrected = df_final_nonpoly.drop(columns=['PDB_chain', 'numbers', "available_chain"])
    df_final_nonpoly_corrected.loc[:, 'renum_str'] = df_final_nonpoly_corrected.PDB.map(
        lambda x: x[1] + "{:>2}".format(x[2])) + df_final_nonpoly_corrected.available_number.map(lambda x: "{:>4}".format(str(int(x)))) + " "

    return df_final_nonpoly_corrected


def remark_465(missing_res_remark_465, df_PDBe_PDB_UniProt):
    df_PDBe_PDB_UniProt_nulls = df_PDBe_PDB_UniProt.loc[df_PDBe_PDB_UniProt['PDB'].apply(lambda x: x[0] == "null")]
    df_PDBe_PDB_UniProt_nulls = df_PDBe_PDB_UniProt_nulls.reset_index(drop=True)

    df_mmCIF_remark_465 = pd.DataFrame(list(zip(missing_res_remark_465)))
    df_mmCIF_remark_465 = df_mmCIF_remark_465.rename(columns={0: "PDB_old"})
    df_mmCIF_remark_465 = df_mmCIF_remark_465.drop_duplicates()

    df_remark_465_final = df_mmCIF_remark_465.merge(df_PDBe_PDB_UniProt_nulls, left_index=True, right_index=True)
    df_remark_465_final.loc[:, 'PDB_str'] = df_remark_465_final.PDB_old.map(
        lambda x: x[1] + "{:>2}".format(x[2]) + "{:>6}".format(x[0]) + " " if x[0].isdigit() else x[1] + "{:>2}".format(x[2]) + "{:>5}".format(x[0]))

    df_final_poly_remark_465_corrected = df_remark_465_final.drop(columns=['UniProt', 'AccessionID', "new_col_Uni", "UniProt_5k", "Uni_moD"])
    df_final_poly_remark_465_corrected.loc[:, 'renum_str'] = df_final_poly_remark_465_corrected.PDB_old.map(
        lambda x: x[1] + "{:>2}".format(x[2])) + df_final_poly_remark_465_corrected.Three_Rows_CIF_Num_Uni.map(
        lambda x: "{:>6}".format(str(int(x[1])))) + " "
    df_final_poly_remark_465_corrected = df_final_poly_remark_465_corrected.reset_index(drop=True)

    return df_final_poly_remark_465_corrected


def final_dict_formation(df_final_poly_corrected, df_final_nonpoly_corrected, final_remark_465, chains_to_change):
    all_data_df = df_final_poly_corrected.append(df_final_nonpoly_corrected, ignore_index=True, sort=False)
    all_data_df = all_data_df.append(final_remark_465, ignore_index=True, sort=False)
    all_data_df_drop_dup = all_data_df.drop_duplicates(subset="PDB_str", keep='first')
    not_in_chains_to_change = all_data_df_drop_dup.PDB.map(lambda x: x if x[2] in chains_to_change else "not_in_chains_to_change")
    all_data_merged_not_in_chain_to_change = all_data_df_drop_dup.merge(not_in_chains_to_change.rename('not_in_chains_to_change'),
                                                                        left_index=True, right_index=True)
    all_data_df_drop_dup_drop_chains = all_data_merged_not_in_chain_to_change[
        all_data_merged_not_in_chain_to_change.not_in_chains_to_change != 'not_in_chains_to_change']

    different_indent_PDB_str = list()
    for n in all_data_df_drop_dup_drop_chains["PDB_str"]:
        if "-" in n:
            n = n[:5] + n[6:] + " "
        different_indent_PDB_str.append(n)
        different_indent_PDB_str.append(n[:3] + " " + n[3:])  # present HET
        different_indent_PDB_str.append(n[:5] + " " + n[5:])  # present most common
        different_indent_PDB_str.append(n[:5] + "  " + n[5:])  # present REMARK 500

    different_indent_renum_str = list()
    for n in all_data_df_drop_dup_drop_chains["renum_str"]:
        different_indent_renum_str.append(n)
        different_indent_renum_str.append(n[:3] + " " + n[3:])  # present HET
        different_indent_renum_str.append(n[:5] + " " + n[5:])  # present most common
        different_indent_renum_str.append(n[:5] + "  " + n[5:])  # present REMARK 500

    dict_for_replacement = dict(zip(different_indent_PDB_str, different_indent_renum_str))
    return dict_for_replacement


def replace_all(lines, dict_for_replacement):
    location_of_the_value = 0
    for key, value in dict_for_replacement.items():
        if key in lines:
            if location_of_the_value == lines.find(key):
                continue
            lines = lines.replace(key, value)
            location_of_the_value = lines.find(value)
    return lines


def master_PDB_renumber_function(input_PDB_files_were_found, default_input_path_to_PDB, default_input_path_to_SIFTS,
                                 default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs):
    if not os.path.exists(default_output_path_to_PDB):
        os.makedirs(default_output_path_to_PDB)

    input_PDB_files_were_found_list = list()
    input_PDB_files_were_found_list.append(input_PDB_files_were_found)

    for PDB in input_PDB_files_were_found_list:

        try:
            assembly_num = re.search('\.pdb(.*).gz', PDB).group(1)
            SIFTS_name = PDB[:4] + ".xml.gz"
            PDB_id = PDB[:4]
        except AttributeError:
            assembly_num = ""
            SIFTS_name = PDB[3:7] + ".xml.gz"
            PDB_id = PDB[3:7]

        # for no corresponding SIFTS
        try:
            gzip.open(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name), 'rt')
        except FileNotFoundError:
            copy_file(default_input_path_to_PDB, PDB, default_output_path_to_PDB, "_renum.pdb" + assembly_num, gzip_mode)
            log_message = if_no_SIFTS_data_log_for_PDB(default_input_path_to_PDB, PDB_id, PDB)
            return log_message

        # for zero byte SIFTS
        if os.path.getsize(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name)) == 0:
            copy_file(default_input_path_to_PDB, PDB, default_output_path_to_PDB, "_renum.pdb" + assembly_num, gzip_mode)
            log_message = if_no_SIFTS_data_log_for_PDB(default_input_path_to_PDB, PDB_id, PDB)
            return log_message

        product_tree_SIFTS = try_SIFTS_tree_parser(default_input_path_to_SIFTS, SIFTS_name)
        if product_tree_SIFTS == 0:
            continue

        tuple_PDBe_for_PDB_and_tuple_PDB = product_tree_SIFTS[0]
        tuple_PDBe_for_UniProt_and_tuple_UniProt = product_tree_SIFTS[1]
        UniProt_conversion_dict = product_tree_SIFTS[2]

        # for no UniProt in SIFTS
        if len(tuple_PDBe_for_UniProt_and_tuple_UniProt) == 0:
            copy_file(default_input_path_to_PDB, PDB, default_output_path_to_PDB, "_renum.pdb" + assembly_num, gzip_mode)
            log_message = if_no_SIFTS_data_log_for_PDB(default_input_path_to_PDB, PDB_id, PDB)
            return log_message

        split = try_PDB(default_input_path_to_PDB, PDB)
        if split == 0:
            continue

        product_of_SIFTS_data_parser = SIFTS_data_parser_for_PDB(tuple_PDBe_for_PDB_and_tuple_PDB,
                                                                 tuple_PDBe_for_UniProt_and_tuple_UniProt,
                                                                 default_PDB_num, 'all')
        df_PDBe_PDB_UniProt = product_of_SIFTS_data_parser[1]

        handling_chain_numbering = handling_chain_numbering_clashes(df_PDBe_PDB_UniProt, exception_AccessionIDs)
        chains_to_change = handling_chain_numbering[0]
        combined_tuple_PDBe_UniProt_AccessionID = handling_chain_numbering[1]
        longest_AccessionID_list = handling_chain_numbering[3]
        chains_to_change_one_to_end = handling_chain_numbering[4]

        product_of_SIFTS_data_parser = SIFTS_data_parser_for_PDB(tuple_PDBe_for_PDB_and_tuple_PDB, combined_tuple_PDBe_UniProt_AccessionID,
                                                                 default_PDB_num, chains_to_change)
        df_PDBe_PDB_UniProt_without_null_index_PDBe = product_of_SIFTS_data_parser[0]
        df_PDBe_PDB_UniProt = product_of_SIFTS_data_parser[1]

        renumbered_count = renumbered_count_in_chains(chains_to_change_one_to_end, df_PDBe_PDB_UniProt_without_null_index_PDBe,
                                                      PDB_id, UniProt_conversion_dict, longest_AccessionID_list)
        chain_total_renum = renumbered_count[0]
        nothing_changed = renumbered_count[1]

        chain_total_renum.append(nothing_changed)
        mod_log_message = chain_total_renum

        if nothing_changed == 0:
            copy_file(default_input_path_to_PDB, PDB, default_output_path_to_PDB, "_renum.pdb" + assembly_num, gzip_mode)
            return mod_log_message

        parsed_PDB = PDB_parser(split, df_PDBe_PDB_UniProt_without_null_index_PDBe, default_PDB_num)
        pd_series_index_PDBe = parsed_PDB[0]
        num_ins_code_name_chain_HETATM = parsed_PDB[1]
        df_final_poly_corrected = parsed_PDB[2]
        missing_res_remark_465 = parsed_PDB[3]
        Num_is_Too_Big = parsed_PDB[4]

        # when numbers get too big
        if not Num_is_Too_Big:
            copy_file(default_input_path_to_PDB, PDB, default_output_path_to_PDB, "_renum.pdb" + assembly_num, gzip_mode)
            return mod_log_message

        df_final_nonpoly_corrected = non_poly_num(pd_series_index_PDBe, num_ins_code_name_chain_HETATM)
        if len(missing_res_remark_465) != 0:
            final_remark_465 = remark_465(missing_res_remark_465, df_PDBe_PDB_UniProt)
        else:
            final_remark_465 = None

        dict_for_replacement = final_dict_formation(df_final_poly_corrected, df_final_nonpoly_corrected, final_remark_465, chains_to_change)

        outF = open(Path(str(default_output_path_to_PDB) + "/" + PDB_id + ".pdb" + assembly_num), "w")
        # PDBrenum REMARK 0 insert
        start_remark_0 = 0
        for lines in split:
            if lines.startswith("HEADER"):
                start_remark_0 += 1
            else:
                split = split[:start_remark_0] + PDBrenum_REMARK + split[start_remark_0:]
                break

        # actual renumbering
        for lines in split:
            lines = replace_all(lines, dict_for_replacement)
            outF.write(lines)
            outF.write("\n")
        outF.close()

        if gzip_mode == "on":
            with open(Path(str(default_output_path_to_PDB) + "/" + PDB_id + ".pdb" + assembly_num), 'rb') as f_in:
                with gzip.open(Path(str(default_output_path_to_PDB) + "/" + PDB_id + "_renum.pdb" + assembly_num + ".gz"), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(Path(str(default_output_path_to_PDB) + "/" + PDB_id + ".pdb" + assembly_num))

        return mod_log_message


def ProcessPool_run_renum_PDB(format_to_download, input_PDB_files_were_found, default_input_path_to_PDB, default_input_path_to_SIFTS,
                              default_output_path_to_PDB, default_PDB_num, gzip_mode, exception_AccessionIDs, nproc):
    if not os.path.exists(default_output_path_to_PDB):
        os.makedirs(default_output_path_to_PDB)

    resulting = list()
    executor = ProcessPoolExecutor(max_workers=nproc)
    partial_master_PDB_renumber_function = partial(master_PDB_renumber_function,
                                                   default_input_path_to_PDB=default_input_path_to_PDB,
                                                   default_input_path_to_SIFTS=default_input_path_to_SIFTS,
                                                   default_output_path_to_PDB=default_output_path_to_PDB,
                                                   default_PDB_num=default_PDB_num, gzip_mode=gzip_mode,
                                                   exception_AccessionIDs=exception_AccessionIDs)

    jobs = [executor.submit(partial_master_PDB_renumber_function, pdb_files) for pdb_files in input_PDB_files_were_found]

    for job in tqdm.tqdm(as_completed(jobs), total=len(jobs), position=0, miniters=1,
                         leave=True, desc="Renumbering " + format_to_download + " files"):
        result = job.result()
        resulting.append(result)

    return resulting


# from src.download.lookfilesinside import look_what_is_inside
# current_directory = "/home/bulat/Desktop/main/PDB_fix_project/PDBarnum"
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
#
# #
# # current_directory = os.getcwd()
# # default_mmCIF_num = 50000
# default_PDB_num = 5000
# gzip_mode = "on"
# exception_AccessionIDs = ["P42212", "Q17104", "Q27903", "Q93125", "P03069", "D3DLN9", "Q96UT3", "P0ABE7", "P00192",
#                           "P76805", "Q8XCE3", "P00720", "Q38170", "Q94N07", "P0AEX9", "P02928", "Q2M6S0"]

# input_PDB_files_were_found = look_what_is_inside("PDB", default_input_path_to_PDB=default_input_path_to_PDB)

# if __name__ == '__main__':
#     resulting2 = ProcessPool_run_renum_PDB("PDB", input_PDB_files_were_found,
#                                            default_input_path_to_PDB,
#                                            default_input_path_to_SIFTS,
#                                            default_output_path_to_PDB,
#                                            default_PDB_num, gzip_mode, exception_AccessionIDs)
