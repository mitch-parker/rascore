from src.download.modules import *

from src.renum.shared.handling_chain_numbering_clashes import handling_chain_numbering_clashes
from src.renum.shared.renumbered_count_in_chains import renumbered_count_in_chains
from src.download.downloadwithThreadPool import download_with_pool, url_formation_for_pool
from src.download import compressor


def try_MMCIF2Dict(default_input_path_to_mmCIF, mmCIF_name):
    mmcif_dict = 0
    for _ in range(3):
        try:
            mmcif_dict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(gzip.open(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_name), 'rt'))
            break
        except EOFError:
            os.remove(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_name))
            if "assembly" in mmCIF_name:
                download_with_pool(url_formation_for_pool("mmCIF_assembly", [mmCIF_name])[0])
            else:
                download_with_pool(url_formation_for_pool("mmCIF", [mmCIF_name])[0])
        except ValueError:
            os.remove(Path(str(default_input_path_to_mmCIF) + "/" + mmCIF_name))
            if "assembly" in mmCIF_name:
                download_with_pool(url_formation_for_pool("mmCIF_assembly", [mmCIF_name])[0])
            else:
                download_with_pool(url_formation_for_pool("mmCIF", [mmCIF_name])[0])
        except OSError:
            if "assembly" in mmCIF_name:
                download_with_pool(url_formation_for_pool("mmCIF_assembly", [mmCIF_name])[0])
            else:
                download_with_pool(url_formation_for_pool("mmCIF", [mmCIF_name])[0])
    return mmcif_dict


def try_SIFTS_tree_parser(default_input_path_to_SIFTS, SIFTS_name):
    product_tree_SIFTS = 0
    for _ in range(3):
        try:
            handle_SIFTS = gzip.open(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name), 'rt')
            product_tree_SIFTS = SIFTS_tree_parser(handle_SIFTS)
            break
        except EOFError:
            os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name))
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
        except ValueError:
            os.remove(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name))
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
        except OSError:
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
        except Exception:
            download_with_pool(url_formation_for_pool("SIFTS", [SIFTS_name])[0])
    return product_tree_SIFTS


def output_with_this_name_ending(name_ending, path, mmcif_dict, mmCIF_name, gzip_mode, current_directory):
    mmCIF_name = mmCIF_name[:mmCIF_name.rfind(".cif.gz")]
    os.chdir(path)
    io = MMCIFIO()
    io.set_dict(mmcif_dict)
    io.save(mmCIF_name + name_ending)
    if gzip_mode == "on":
        compressor.compress_output_files(mmCIF_name + name_ending, gzip_mode)
        os.remove(mmCIF_name + name_ending)
    os.chdir(current_directory)


def copy_file(inpath, file_name, outpath, postfix, gzip_mode):
    mmCIF_name = file_name[:file_name.rfind(".cif.gz")]
    absolute_path_in = inpath + "/" + file_name
    absolute_path_out = outpath + "/" + mmCIF_name + postfix
    if gzip_mode == "off":
        with gzip.open(absolute_path_in, 'rb') as f_in:
            with open(absolute_path_out[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(absolute_path_in, absolute_path_out)


def if_no_SIFTS_data_log(mmCIF_name, mmcif_dict, log_message):
    strand_id_set = set()
    try:
        pull_chains_for_chains_count = mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"]
    except KeyError:
        try:
            pull_chains_for_chains_count = mmcif_dict["_pdbe_orig_poly_seq_scheme.pdb_strand_id"]
        except KeyError:
            pull_chains_for_chains_count = mmcif_dict["_atom_site.auth_asym_id"]

    for strand in pull_chains_for_chains_count:
        strand_id_set.add(strand)
    strand_id_set = list(strand_id_set)
    strand_id_set.sort()
    for strand in strand_id_set:
        count_elements_in_strand = 0
        for chain_id in pull_chains_for_chains_count:
            if chain_id == strand:
                count_elements_in_strand += 1
        log_message.append([mmCIF_name[:4], strand, "-", "-", "-", "-", count_elements_in_strand, "0", "0"])
    return log_message


def renum_struct_ref_seq_pdbx_auth_seq_align(mmcif_dict):
    try:
        _struct_ref_seq_pdbx_strand_id = mmcif_dict["_struct_ref_seq.pdbx_strand_id"]

        _struct_ref_seq_pdbx_seq_align_beg_ins_code = mmcif_dict["_struct_ref_seq.pdbx_seq_align_beg_ins_code"]
        _struct_ref_seq_pdbx_auth_seq_align_beg = mmcif_dict["_struct_ref_seq.pdbx_auth_seq_align_beg"]
        _struct_ref_seq_db_align_beg = mmcif_dict["_struct_ref_seq.db_align_beg"]
        mmcif_dict["_struct_ref_seq.pdbx_auth_seq_align_beg"] = mmcif_dict["_struct_ref_seq.db_align_beg"]

        _struct_ref_seq_pdbx_seq_align_end_ins_code = mmcif_dict["_struct_ref_seq.pdbx_seq_align_end_ins_code"]
        _struct_ref_seq_pdbx_auth_seq_align_end = mmcif_dict["_struct_ref_seq.pdbx_auth_seq_align_end"]
        _struct_ref_seq_db_align_end = mmcif_dict["_struct_ref_seq.db_align_end"]
        mmcif_dict["_struct_ref_seq.pdbx_auth_seq_align_end"] = mmcif_dict["_struct_ref_seq.db_align_end"]

        if type(_struct_ref_seq_pdbx_seq_align_beg_ins_code) == str:
            if "." in _struct_ref_seq_pdbx_seq_align_beg_ins_code:
                mmcif_dict["_struct_ref_seq.pdbx_seq_align_beg_ins_code"] = "."
            else:
                mmcif_dict["_struct_ref_seq.pdbx_seq_align_beg_ins_code"] = "?"
        if type(_struct_ref_seq_pdbx_seq_align_end_ins_code) == str:
            if "." in _struct_ref_seq_pdbx_seq_align_end_ins_code:
                mmcif_dict["_struct_ref_seq.pdbx_seq_align_end_ins_code"] = "."
            else:
                mmcif_dict["_struct_ref_seq.pdbx_seq_align_end_ins_code"] = "?"

        PDB_ins_code_list = list()
        if type(_struct_ref_seq_pdbx_seq_align_beg_ins_code) != str:
            if "." in _struct_ref_seq_pdbx_seq_align_beg_ins_code:
                for _ in range(len(_struct_ref_seq_pdbx_seq_align_beg_ins_code)):
                    PDB_ins_code_list.append(".")
            else:
                for _ in range(len(_struct_ref_seq_pdbx_seq_align_beg_ins_code)):
                    PDB_ins_code_list.append("?")
            mmcif_dict["_struct_ref_seq.pdbx_seq_align_beg_ins_code"] = PDB_ins_code_list
            mmcif_dict["_struct_ref_seq.pdbx_seq_align_end_ins_code"] = PDB_ins_code_list

    except KeyError:
        pass


def poly_nonpoly_renum(mmcif_dict, df_PDBe_PDB_UniProt, chains_to_change, default_mmCIF_num):
    try:
        _pdbx_poly_seq_scheme_seq_id = mmcif_dict["_pdbx_poly_seq_scheme.seq_id"]
        _pdbx_poly_seq_scheme_asym_id = mmcif_dict["_pdbx_poly_seq_scheme.asym_id"]
        _pdbx_poly_seq_scheme_mon_id = mmcif_dict["_pdbx_poly_seq_scheme.mon_id"]

        _pdbx_poly_seq_scheme_pdb_seq_num = mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"]
        _pdbx_poly_seq_scheme_auth_seq_num = mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"]
        _pdbx_poly_seq_scheme_pdb_mon_id = mmcif_dict["_pdbx_poly_seq_scheme.pdb_mon_id"]
        _pdbx_poly_seq_scheme_auth_mon_id = mmcif_dict["_pdbx_poly_seq_scheme.auth_mon_id"]
        _pdbx_poly_seq_scheme_pdb_strand_id = mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"]
        _pdbx_poly_seq_scheme_pdb_ins_code = mmcif_dict["_pdbx_poly_seq_scheme.pdb_ins_code"]
    except KeyError:
        try:
            _pdbx_poly_seq_scheme_seq_id = mmcif_dict["_pdbe_orig_poly_seq_scheme.seq_id"]
            _pdbx_poly_seq_scheme_asym_id = mmcif_dict["_pdbe_orig_poly_seq_scheme.asym_id"]
            _pdbx_poly_seq_scheme_mon_id = mmcif_dict["_pdbe_orig_poly_seq_scheme.mon_id"]

            _pdbx_poly_seq_scheme_pdb_seq_num = mmcif_dict["_pdbe_orig_poly_seq_scheme.pdb_seq_num"]
            _pdbx_poly_seq_scheme_auth_seq_num = mmcif_dict["_pdbe_orig_poly_seq_scheme.auth_seq_num"]
            _pdbx_poly_seq_scheme_pdb_mon_id = mmcif_dict["_pdbe_orig_poly_seq_scheme.pdb_mon_id"]
            _pdbx_poly_seq_scheme_auth_mon_id = mmcif_dict["_pdbe_orig_poly_seq_scheme.auth_mon_id"]
            _pdbx_poly_seq_scheme_pdb_strand_id = mmcif_dict["_pdbe_orig_poly_seq_scheme.pdb_strand_id"]
            _pdbx_poly_seq_scheme_pdb_ins_code = mmcif_dict["_pdbe_orig_poly_seq_scheme.pdb_ins_code"]

        except KeyError:
            # continue
            return 0

    if type(_pdbx_poly_seq_scheme_pdb_strand_id) == str:
        _pdbx_poly_seq_scheme_pdb_seq_num = [_pdbx_poly_seq_scheme_pdb_seq_num]
        _pdbx_poly_seq_scheme_auth_seq_num = [_pdbx_poly_seq_scheme_auth_seq_num]
        _pdbx_poly_seq_scheme_pdb_mon_id = [_pdbx_poly_seq_scheme_pdb_mon_id]
        _pdbx_poly_seq_scheme_auth_mon_id = [_pdbx_poly_seq_scheme_auth_mon_id]
        _pdbx_poly_seq_scheme_pdb_strand_id = [_pdbx_poly_seq_scheme_pdb_strand_id]
        _pdbx_poly_seq_scheme_pdb_ins_code = [_pdbx_poly_seq_scheme_pdb_ins_code]

    mmCIF_pdbx_poly_seq_scheme_label = list(zip(_pdbx_poly_seq_scheme_seq_id,
                                                _pdbx_poly_seq_scheme_mon_id,
                                                _pdbx_poly_seq_scheme_asym_id))
    mmCIF_pdbx_poly_seq_scheme_pdb = list(zip(_pdbx_poly_seq_scheme_pdb_seq_num,
                                              _pdbx_poly_seq_scheme_pdb_mon_id,
                                              _pdbx_poly_seq_scheme_pdb_strand_id))
    mmCIF_pdbx_poly_seq_scheme_auth = list(zip(_pdbx_poly_seq_scheme_auth_seq_num,
                                               _pdbx_poly_seq_scheme_auth_mon_id,
                                               _pdbx_poly_seq_scheme_pdb_strand_id))

    df_mmCIF_pdbx_poly_seq_scheme = pd.DataFrame(zip(mmCIF_pdbx_poly_seq_scheme_label,
                                                     mmCIF_pdbx_poly_seq_scheme_pdb,
                                                     mmCIF_pdbx_poly_seq_scheme_auth,
                                                     _pdbx_poly_seq_scheme_pdb_ins_code))

    df_mmCIF_pdbx_poly_seq_scheme = df_mmCIF_pdbx_poly_seq_scheme.rename(
        columns={0: "_pdbx_poly_seq_scheme_label", 1: "pdbx_poly_seq_scheme_pdb",
                 2: "pdbx_poly_seq_scheme_auth", 3: "pdbx_poly_seq_scheme_pdb_ins_code"})

    df_pdbx_poly_seq_scheme_pdb_final = df_mmCIF_pdbx_poly_seq_scheme.merge(
        df_PDBe_PDB_UniProt, left_on="_pdbx_poly_seq_scheme_label", right_on="PDBe", how='left')
    df_pdbx_poly_seq_scheme_pdb_final["PDBe_num_and_chain"] = df_pdbx_poly_seq_scheme_pdb_final[
        "_pdbx_poly_seq_scheme_label"].apply(lambda x: (x[0], x[2]))

    df_pdbx_poly_seq_scheme_pdb_final["PDB_num_and_chain"] = np.where(
        df_pdbx_poly_seq_scheme_pdb_final["pdbx_poly_seq_scheme_pdb_ins_code"].apply(lambda x: x == "."),
        df_pdbx_poly_seq_scheme_pdb_final["pdbx_poly_seq_scheme_pdb"].apply(lambda x: (x[0], x[2])),
        df_pdbx_poly_seq_scheme_pdb_final["pdbx_poly_seq_scheme_pdb"].apply(lambda x: x[0]) +
        df_pdbx_poly_seq_scheme_pdb_final["pdbx_poly_seq_scheme_pdb_ins_code"].apply(lambda x: x) + "," +
        df_pdbx_poly_seq_scheme_pdb_final["pdbx_poly_seq_scheme_pdb"].apply(lambda x: x[2]))
    df_pdbx_poly_seq_scheme_pdb_final["PDB_num_and_chain"] = df_pdbx_poly_seq_scheme_pdb_final["PDB_num_and_chain"].apply(
        lambda x: tuple(x.split(",")) if type(x) == str else x)

    df_pdbx_poly_seq_scheme_pdb_final["Uni_or_50k"] = np.where(
        df_pdbx_poly_seq_scheme_pdb_final["PDB_num_and_chain"].apply(lambda x: x[1] in chains_to_change),
        df_pdbx_poly_seq_scheme_pdb_final["UniProt_50k"].apply(lambda x: x),
        df_pdbx_poly_seq_scheme_pdb_final["PDB_num_and_chain"].apply(lambda x: x[0].strip(re.sub('[0-9\-\?\.]+', '', x[0]))))

    try:
        mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"]  # check if key exists
        mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"] = list(df_pdbx_poly_seq_scheme_pdb_final["Uni_or_50k"].values)
        mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"] = _pdbx_poly_seq_scheme_pdb_seq_num
    except KeyError:
        mmcif_dict["_pdbe_orig_poly_seq_scheme.pdb_seq_num"] = list(df_pdbx_poly_seq_scheme_pdb_final["Uni_or_50k"].values)
        mmcif_dict["_pdbe_orig_poly_seq_scheme.auth_seq_num"] = _pdbx_poly_seq_scheme_pdb_seq_num

    nonpoly_present = False

    try:
        _pdbx_nonpoly_scheme_pdb_seq_num = mmcif_dict["_pdbx_nonpoly_scheme.pdb_seq_num"]
        _pdbx_nonpoly_scheme_auth_seq_num = mmcif_dict["_pdbx_nonpoly_scheme.auth_seq_num"]
        _pdbx_nonpoly_scheme_pdb_mon_id = mmcif_dict["_pdbx_nonpoly_scheme.pdb_mon_id"]
        _pdbx_nonpoly_scheme_auth_mon_id = mmcif_dict["_pdbx_nonpoly_scheme.auth_mon_id"]
        _pdbx_nonpoly_scheme_pdb_strand_id = mmcif_dict["_pdbx_nonpoly_scheme.pdb_strand_id"]
        _pdbx_nonpoly_scheme_asym_id = mmcif_dict["_pdbx_nonpoly_scheme.asym_id"]
        dots_for_label = ["." for _ in range(len(_pdbx_nonpoly_scheme_asym_id)) if type(_pdbx_nonpoly_scheme_asym_id) == list]
        nonpoly_present = True
    except KeyError:
        try:
            _pdbx_nonpoly_scheme_pdb_seq_num = mmcif_dict["_pdbe_orig_nonpoly_scheme.pdb_seq_num"]
            _pdbx_nonpoly_scheme_auth_seq_num = mmcif_dict["_pdbe_orig_nonpoly_scheme.auth_seq_num"]
            _pdbx_nonpoly_scheme_pdb_mon_id = mmcif_dict["_pdbe_orig_nonpoly_scheme.pdb_mon_id"]
            _pdbx_nonpoly_scheme_auth_mon_id = mmcif_dict["_pdbe_orig_nonpoly_scheme.auth_mon_id"]
            _pdbx_nonpoly_scheme_pdb_strand_id = mmcif_dict["_pdbe_orig_nonpoly_scheme.pdb_strand_id"]
            _pdbx_nonpoly_scheme_asym_id = mmcif_dict["_pdbe_orig_nonpoly_scheme.asym_id"]
            dots_for_label = ["." for _ in range(len(_pdbx_nonpoly_scheme_asym_id)) if type(_pdbx_nonpoly_scheme_asym_id) == list]
            nonpoly_present = True
        except KeyError:
            pass

    if nonpoly_present:
        if type(_pdbx_nonpoly_scheme_pdb_strand_id) == str:
            _pdbx_nonpoly_scheme_pdb_seq_num = [_pdbx_nonpoly_scheme_pdb_seq_num]
            _pdbx_nonpoly_scheme_auth_seq_num = [_pdbx_nonpoly_scheme_auth_seq_num]
            _pdbx_nonpoly_scheme_pdb_mon_id = [_pdbx_nonpoly_scheme_pdb_mon_id]
            _pdbx_nonpoly_scheme_auth_mon_id = [_pdbx_nonpoly_scheme_auth_mon_id]
            _pdbx_nonpoly_scheme_pdb_strand_id = [_pdbx_nonpoly_scheme_pdb_strand_id]
            _pdbx_nonpoly_scheme_asym_id = [_pdbx_nonpoly_scheme_asym_id]
            dots_for_label = ["."]

        mmCIF_pdbx_nonpoly_scheme_pdb = list(zip(_pdbx_nonpoly_scheme_pdb_seq_num,
                                                 _pdbx_nonpoly_scheme_pdb_mon_id,
                                                 _pdbx_nonpoly_scheme_pdb_strand_id))
        mmCIF_pdbx_nonpoly_scheme_auth = list(zip(_pdbx_nonpoly_scheme_auth_seq_num,
                                                  _pdbx_nonpoly_scheme_auth_mon_id,
                                                  _pdbx_nonpoly_scheme_pdb_strand_id))
        mmCIF_pdbx_nonpoly_scheme_label = list(zip(dots_for_label,
                                                   _pdbx_nonpoly_scheme_pdb_mon_id,
                                                   _pdbx_nonpoly_scheme_asym_id))

        df_mmCIF_pdbx_nonpoly_scheme = pd.DataFrame(zip(mmCIF_pdbx_nonpoly_scheme_pdb,
                                                        mmCIF_pdbx_nonpoly_scheme_auth,
                                                        mmCIF_pdbx_nonpoly_scheme_label))
        df_mmCIF_pdbx_nonpoly_scheme = df_mmCIF_pdbx_nonpoly_scheme.rename(columns={0: "pdbx_nonpoly_scheme_pdb",
                                                                                    1: "pdbx_nonpoly_scheme_auth",
                                                                                    2: "pdbx_nonpoly_scheme_label"})

        df_mmCIF_pdbx_nonpoly_scheme["PDB"] = df_mmCIF_pdbx_nonpoly_scheme["pdbx_nonpoly_scheme_pdb"]
        df_mmCIF_pdbx_nonpoly_scheme["PDB_num_and_chain"] = df_mmCIF_pdbx_nonpoly_scheme["pdbx_nonpoly_scheme_pdb"].apply(lambda x: (x[0], x[2]))
        df_mmCIF_pdbx_nonpoly_scheme["Uni_or_50k"] = df_mmCIF_pdbx_nonpoly_scheme["pdbx_nonpoly_scheme_pdb"].apply(
            lambda x: str(int(x[0]) + default_mmCIF_num + 10000) if x[2] in chains_to_change else x[0])

        try:
            mmcif_dict["_pdbx_nonpoly_scheme.pdb_seq_num"]  # check if key exists
            mmcif_dict["_pdbx_nonpoly_scheme.pdb_seq_num"] = list(df_mmCIF_pdbx_nonpoly_scheme["Uni_or_50k"].values)
            mmcif_dict["_pdbx_nonpoly_scheme.auth_seq_num"] = _pdbx_nonpoly_scheme_pdb_seq_num
        except KeyError:
            try:
                mmcif_dict["_pdbe_orig_nonpoly_scheme.pdb_seq_num"] = list(df_mmCIF_pdbx_nonpoly_scheme["Uni_or_50k"].values)
                mmcif_dict["_pdbe_orig_nonpoly_scheme.auth_seq_num"] = _pdbx_nonpoly_scheme_pdb_seq_num
            except KeyError:
                pass

        poly_nonpoly_append = df_pdbx_poly_seq_scheme_pdb_final.append(df_mmCIF_pdbx_nonpoly_scheme)
        poly_nonpoly_append = poly_nonpoly_append[["PDBe", "PDB", "UniProt", "PDBe_num_and_chain", "PDB_num_and_chain", "AccessionID", "Uni_or_50k"]]
    else:
        poly_nonpoly_append = df_pdbx_poly_seq_scheme_pdb_final[
            ["PDBe", "PDB", "UniProt", "PDBe_num_and_chain", "PDB_num_and_chain", "AccessionID", "Uni_or_50k"]]

    return poly_nonpoly_append


def renumber_tables(formed_columns, mmcif_dict, poly_nonpoly_atom_site, chains_to_change, default_mmCIF_num):
    dot_or_question_tuple = (".", "?")
    for n in formed_columns:
        auth_comp_id = 0
        auth_seq_id = n[0]
        auth_asym_id = n[1]
        try:
            PDB_ins_code = n[2]
            if "ins_code" not in PDB_ins_code:
                auth_comp_id = PDB_ins_code
                PDB_ins_code = 0
        except IndexError:
            PDB_ins_code = 0
        try:
            if auth_comp_id == 0:
                auth_comp_id = n[3]
        except IndexError:
            auth_comp_id = 0

        if "_pdbx_branch_scheme" in auth_seq_id:
            auth_seq_id = "_pdbx_branch_scheme.pdb_seq_num"
            auth_asym_id = "_pdbx_branch_scheme.pdb_asym_id"

        PDB_ins_code_list = list()
        # auth_comp_id_list = mmcif_dict[auth_comp_id] #for debug only
        auth_seq_id_list = mmcif_dict[auth_seq_id]
        auth_asym_id_list = mmcif_dict[auth_asym_id]

        if PDB_ins_code == 0:
            for _ in range(len(auth_seq_id_list)):
                PDB_ins_code_list.append("?")
        else:
            PDB_ins_code_list = mmcif_dict[PDB_ins_code]

        if type(auth_asym_id_list) == str:
            # auth_comp_id_list = [auth_comp_id_list] for debug only
            auth_seq_id_list = [auth_seq_id_list]
            auth_asym_id_list = [auth_asym_id_list]

            if PDB_ins_code == 0:
                PDB_ins_code_list = ["?"]
            else:
                PDB_ins_code_list = [PDB_ins_code]

        if PDB_ins_code != 0:
            dot_to_question = list()
            for ins_code in mmcif_dict[PDB_ins_code]:
                if ins_code == ".":
                    dot_to_question.append("?")
                else:
                    dot_to_question.append(ins_code)
            PDB_ins_code_list = dot_to_question

        auth_seq_id_list_zip = list(zip(auth_seq_id_list, auth_asym_id_list))
        df_mmCIF_auth_seq_id_list_zip = pd.DataFrame(zip(auth_seq_id_list_zip, PDB_ins_code_list))
        df_mmCIF_auth_seq_id_list_zip = df_mmCIF_auth_seq_id_list_zip.rename(columns={0: "auth_seq_id_list_zip", 1: "ins_code"})

        df_mmCIF_auth_seq_id_list_zip["PDB_with_ins_code"] = np.where(df_mmCIF_auth_seq_id_list_zip['ins_code'] != "?",
                                                                      (df_mmCIF_auth_seq_id_list_zip['auth_seq_id_list_zip'].apply(lambda x: x[0])
                                                                       + df_mmCIF_auth_seq_id_list_zip['ins_code'].apply(lambda y: y[0]) + ","
                                                                       + df_mmCIF_auth_seq_id_list_zip['auth_seq_id_list_zip'].apply(lambda x: x[1])),
                                                                      df_mmCIF_auth_seq_id_list_zip['ins_code'])

        df_mmCIF_auth_seq_id_list_zip["PDB_with_ins_code_cor"] = np.where(df_mmCIF_auth_seq_id_list_zip['PDB_with_ins_code'] != "?",
                                                                          df_mmCIF_auth_seq_id_list_zip["PDB_with_ins_code"].apply(
                                                                              lambda x: tuple(x.split(","))),
                                                                          df_mmCIF_auth_seq_id_list_zip["auth_seq_id_list_zip"])

        df_mmCIF_auth_seq_id_list_zip["auth_seq_id_list_zip"] = df_mmCIF_auth_seq_id_list_zip["PDB_with_ins_code_cor"]
        df_mmCIF_auth_seq_id_list_zip = df_mmCIF_auth_seq_id_list_zip.drop(columns=["PDB_with_ins_code_cor", "ins_code", "PDB_with_ins_code"])

        df_auth_seq_id_list_zip_final = df_mmCIF_auth_seq_id_list_zip.merge(poly_nonpoly_atom_site, left_on="auth_seq_id_list_zip",
                                                                            right_on="PDB_num_and_chain", how='left')

        df_auth_seq_id_list_zip_final["question_mark"] = np.where(
            df_auth_seq_id_list_zip_final["auth_seq_id_list_zip"].apply(lambda x: x[0] in dot_or_question_tuple),
            df_auth_seq_id_list_zip_final["auth_seq_id_list_zip"].apply(lambda x: x[0]),
            df_auth_seq_id_list_zip_final["Uni_or_50k"].apply(lambda x: x))
        try:
            df_auth_seq_id_list_zip_final["final"] = np.where(df_auth_seq_id_list_zip_final["question_mark"].apply(lambda x: type(x) == float),
                                                              df_auth_seq_id_list_zip_final["auth_seq_id_list_zip"].apply(
                                                                  lambda x: "." if x[0] == "." else
                                                                  "?" if x[0] == "?" else str(
                                                                      int(''.join(filter(str.isdigit, str(x[0])))) + default_mmCIF_num)
                                                                  if x[1] in chains_to_change else str(int(''.join(filter(str.isdigit, str(x[0])))))),
                                                              df_auth_seq_id_list_zip_final["question_mark"].apply(lambda x: x))
        except ValueError:
            # print("ValueError in table " + auth_seq_id + " has non-numeric value point in file " + mmcif_dict["data_"])
            return print("ValueError in table " + auth_seq_id + " has non-numeric value point in file " + mmcif_dict["data_"])

        df_auth_seq_id_list_zip_final["ins_code"] = df_auth_seq_id_list_zip_final["final"].apply(lambda x: "?"
        if re.sub('[0-9]+', '', x).strip("-").strip(".").strip('?') == ""
        else re.sub('[0-9]+', '', x).strip("-").strip(".").strip('?'))
        df_auth_seq_id_list_zip_final["final"] = df_auth_seq_id_list_zip_final["final"].apply(lambda x: x.strip(re.sub('[0-9\-\?\.]+', '', x)))

        for num in df_auth_seq_id_list_zip_final["final"]:
            if num == "":
                print("Empty str")
            if type(num) == float:
                print("Float or npNAN")

        # actual replacing auth_num with UniProt_num and of ins_code with '?'

        PDB_ins_code_list = list()
        if PDB_ins_code != 0:
            if "." in mmcif_dict[PDB_ins_code]:
                for ins in df_auth_seq_id_list_zip_final["ins_code"].values:
                    if "?" == ins:
                        PDB_ins_code_list.append(".")
                    else:
                        PDB_ins_code_list.append(ins)
                mmcif_dict[PDB_ins_code] = PDB_ins_code_list
            else:
                mmcif_dict[PDB_ins_code] = list(df_auth_seq_id_list_zip_final["ins_code"].values)

        if "_pdbx_branch_scheme" in auth_seq_id:
            mmcif_dict["_pdbx_branch_scheme.auth_seq_num"] = list(df_auth_seq_id_list_zip_final["final"].values)
        else:
            mmcif_dict[auth_seq_id] = list(df_auth_seq_id_list_zip_final["final"].values)

    return mmcif_dict


def column_formation(mmcif_dict):
    mmcif_dict_keys = mmcif_dict.keys()
    aut_seq_all_splitted = list()
    for key in mmcif_dict_keys:
        key_dot_splitted = key.split(".")
        for tab_name_col_name in key_dot_splitted:
            if "auth_seq" in tab_name_col_name:
                if "auth_seq_id" in key:
                    aut_seq_all_splitted.append(key_dot_splitted[:1] + key_dot_splitted[1].split("auth_seq_id"))
                if "auth_seq_num" in key:
                    aut_seq_all_splitted.append(key_dot_splitted[:1] + key_dot_splitted[1].split("auth_seq_num"))

    totaling_combinations = list()
    for table_name_prefix_suffix in aut_seq_all_splitted:
        combinations = list()
        for key in mmcif_dict_keys:
            if table_name_prefix_suffix[0] == key.split(".")[0]:
                # res_num auth_seq_id or auth_seq_num
                if table_name_prefix_suffix[1] in key and table_name_prefix_suffix[2] in key \
                        and "auth_seq_id" in key or "auth_seq_num" in key:
                    combinations.append(key)
                # chain auth_asym_id or strand_id
                if "assembly" in mmcif_dict["data_"]:
                    if table_name_prefix_suffix[1] in key and table_name_prefix_suffix[2] in key \
                            and "orig_auth_asym_id" in key:
                        combinations.append(key)
                else:
                    if table_name_prefix_suffix[1] in key and table_name_prefix_suffix[2] in key \
                            and "auth_asym_id" in key or "strand_id" in key:
                        combinations.append(key)
                # ins_code
                if table_name_prefix_suffix[1] in key and table_name_prefix_suffix[2] in key \
                        and "ins_code" in key:
                    combinations.append(key)
                # monomer_type or auth_comp_id or auth_mon_id or mon_id for _struct_ref_seq_dif
                if table_name_prefix_suffix[1] in key and table_name_prefix_suffix[2] in key \
                        and "auth_comp_id" in key or "auth_mon_id" in key:
                    combinations.append(key)
                elif table_name_prefix_suffix[0] == "_struct_ref_seq_dif" \
                        and "mon_id" in key and "db_mon_id" not in key:
                    combinations.append(key)

        # work assuming all the elements in right order
        # and they are not crossing each other
        if len(combinations) > 4:
            combinations = combinations[:4]

        ordered_combination = list()
        for name in combinations:
            if "auth_seq" in name:
                ordered_combination.insert(0, name)
        for name in combinations:
            if "auth_asym_id" in name or "strand_id" in name:
                ordered_combination.insert(1, name)
        for name in combinations:
            if "ins_code" in name:
                ordered_combination.insert(2, name)
        for name in combinations:
            if "auth_comp_id" in name or "mon_id" in name:
                ordered_combination.insert(3, name)

        # exceptions
        if (  # "pdbx_unobs_or_zero_occ_residues" not in ordered_combination[0]
                "nonpoly_scheme" not in ordered_combination[0]
                and "poly_seq_scheme" not in ordered_combination[0]
                and "ndb_struct_na_base" not in ordered_combination[0]):
            totaling_combinations.append(ordered_combination)

    return totaling_combinations


def mmCIF_parser(mmCIF_name, default_input_path_to_mmCIF, df_PDBe_PDB_UniProt_without_null_index_PDBe, default_mmCIF_num, chains_to_change,
                 chains_to_change_one_to_end):
    mmcif_dict = try_MMCIF2Dict(default_input_path_to_mmCIF, mmCIF_name)
    if mmcif_dict == 0:
        return None

    try:
        _pdbx_poly_seq_scheme_auth_seq_num_before_change = mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"]
    except KeyError:
        _pdbx_poly_seq_scheme_auth_seq_num_before_change = mmcif_dict["_pdbe_orig_poly_seq_scheme.auth_seq_num"]
        pass

    _atom_site_label_comp_id_list = mmcif_dict["_atom_site.label_comp_id"]
    _atom_site_label_seq_id_list = mmcif_dict["_atom_site.label_seq_id"]
    _atom_site_label_asym_id = mmcif_dict["_atom_site.label_asym_id"]
    _atom_site_pdbx_PDB_ins_code = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]

    _atom_site_auth_comp_id = mmcif_dict["_atom_site.auth_comp_id"]
    _atom_site_auth_seq_id = mmcif_dict["_atom_site.auth_seq_id"]
    _atom_site_auth_asym_id = mmcif_dict["_atom_site.auth_asym_id"]
    _atom_site_pdbx_formal_charge = mmcif_dict["_atom_site.pdbx_formal_charge"]

    final_mmCIF_data_list_of_tuples_just_pdb = list(zip(_atom_site_label_seq_id_list, _atom_site_label_comp_id_list, _atom_site_label_asym_id))
    final_mmCIF_data_list_of_tuples_with_auth = list(zip(_atom_site_auth_seq_id, _atom_site_auth_comp_id, _atom_site_auth_asym_id))
    final_mmCIF_data_list_of_tuples_for_df = list(
        zip(final_mmCIF_data_list_of_tuples_just_pdb, final_mmCIF_data_list_of_tuples_with_auth, _atom_site_pdbx_PDB_ins_code))

    df_mmCIF = pd.DataFrame(final_mmCIF_data_list_of_tuples_for_df)
    df_mmCIF = df_mmCIF.rename(columns={0: "One_to_N_mmCIF", 1: "auth_mmCIF", 2: "ins_code"})

    df_mmCIF["One_to_N_mmCIF"]

    df_mmCIF["PDBnum_inc_code"] = np.where(df_mmCIF['ins_code'] != "?",
                                           (df_mmCIF['auth_mmCIF'].apply(lambda x: x[0]) + df_mmCIF["ins_code"].apply(lambda y: y[0]) + ","
                                            + df_mmCIF['auth_mmCIF'].apply(lambda x: x[1]) + "," + df_mmCIF['auth_mmCIF'].apply(lambda x: x[2])),
                                           df_mmCIF["ins_code"])
    df_mmCIF["PDBnum_inc_code_cor"] = np.where(df_mmCIF["PDBnum_inc_code"] != "?", df_mmCIF["PDBnum_inc_code"].apply(lambda x: tuple(x.split(","))),
                                               df_mmCIF["auth_mmCIF"])

    df_mmCIF["auth_mmCIF"] = df_mmCIF["PDBnum_inc_code_cor"]
    df_mmCIF = df_mmCIF.drop(columns=["PDBnum_inc_code_cor", "ins_code", "PDBnum_inc_code"])

    df_PDBe_PDB_UniProt_without_null_index_PDBe = df_PDBe_PDB_UniProt_without_null_index_PDBe.reset_index()
    df_final = df_mmCIF.merge(df_PDBe_PDB_UniProt_without_null_index_PDBe, left_on="One_to_N_mmCIF", right_on="PDBe", how='left')
    df_final = df_final.rename(columns={"PDBe_copy": "PDBe"})
    df_final = df_final.drop_duplicates(subset="auth_mmCIF", keep='first')
    df_final["PDB_num_and_chain"] = df_final["auth_mmCIF"].apply(lambda x: (x[0], x[2]))
    df_final["PDBe_num_and_chain"] = df_final["One_to_N_mmCIF"].apply(lambda x: (x[0], x[2]))

    df_final["Uni_or_50k_NAN"] = np.where(df_final["One_to_N_mmCIF"].apply(lambda x: x[0] != "."),
                                          df_final["UniProt_50k"].apply(lambda x: x),
                                          df_final["PDB_num_and_chain"].apply(
                                              lambda x: str(int(''.join(filter(str.isdigit, x[0]))) + default_mmCIF_num + 10000)
                                              if x[1] in chains_to_change else str(int(''.join(filter(str.isdigit, x[0]))))))
    df_final["Uni_or_50k"] = np.where(df_final["Uni_or_50k_NAN"].apply(lambda x: type(x) == float),
                                      df_final["PDBe_num_and_chain"].apply(
                                          lambda x: "." if x[0] == "." else str(int(''.join(filter(str.isdigit, x[0]))) + default_mmCIF_num)
                                          if x[1] in chains_to_change_one_to_end else str(int(''.join(filter(str.isdigit, x[0]))))),
                                      df_final["Uni_or_50k_NAN"].apply(lambda x: x))

    df_final_atom_site = df_final[["PDBe", "PDB", "UniProt", "PDBe_num_and_chain", "PDB_num_and_chain", "AccessionID", "Uni_or_50k"]]

    return [df_final_atom_site, mmcif_dict]


def SIFTS_tree_parser(handle_SIFTS):
    tree = ET.parse(handle_SIFTS)
    root = tree.getroot()

    crossRefDb_list = list()
    PDBe_val_tuples_in_list = list()
    PDBe_val_tuples_in_list_for_Uni = list()
    PDBe_val_tuples_in_list_for_PDB = list()
    PDB_val_tuples_in_list = list()
    UniProt_val_tuple_in_list = list()
    UniProtdbAccessionId_list = list()
    UniProt_conversion_dict = dict()
    details_list = list()
    # Human_readable_AccessionID_list = list()

    for entity in root:
        if entity.tag.endswith("entity"):
            entity_chainID_list = list(entity.attrib.items())
            if entity_chainID_list[0][0] == "type" and entity_chainID_list[0][1] == "protein":
                for segment in entity:
                    for listResidue in segment:
                        if listResidue.tag.endswith("listMapRegion"):
                            for mapRegion in listResidue:
                                for db in mapRegion:
                                    dbSource_UniProt = list(db.attrib.items())
                                    if "dbSource" == dbSource_UniProt[0][0] and "UniProt" == dbSource_UniProt[0][1]:
                                        if db.text is None:
                                            UniProt = dbSource_UniProt[2][1]
                                        else:
                                            Human_readable = db.text
                                            UniProt_conversion_dict[UniProt] = Human_readable

                        for residue in listResidue:
                            key_val_tuples_in_list_parent = list(residue.attrib.items())
                            if key_val_tuples_in_list_parent[0][0] == "dbSource" and key_val_tuples_in_list_parent[0][1] == "PDBe":
                                PDBe_val_tuples_in_list.append(
                                    (key_val_tuples_in_list_parent[2][1], key_val_tuples_in_list_parent[3][1], entity_chainID_list[1][1]))

                                for crossRefDb in residue:
                                    if crossRefDb.tag.endswith("residueDetail") and crossRefDb.text != "Not_Observed":
                                        details_list.append((("PDBid", root.get("dbAccessionId")), ("Annotation:", crossRefDb.text), (
                                            key_val_tuples_in_list_parent[2][1], key_val_tuples_in_list_parent[3][1], entity_chainID_list[1][1])))

                                    crossRefDb_list.append(crossRefDb.attrib)
                                    key_val_tuples_in_list_child = list(crossRefDb.attrib.items())

                                    if key_val_tuples_in_list_child[0][0] == "dbSource" and key_val_tuples_in_list_child[0][1] == "PDB":
                                        PDB_val_tuples_in_list.append((key_val_tuples_in_list_child[3][1], key_val_tuples_in_list_child[4][1],
                                                                       key_val_tuples_in_list_child[5][1]))
                                        PDBe_val_tuples_in_list_for_PDB.append(
                                            (key_val_tuples_in_list_parent[2][1], key_val_tuples_in_list_parent[3][1], entity_chainID_list[1][1]))

                                    if key_val_tuples_in_list_child[0][0] == "dbSource" and key_val_tuples_in_list_child[0][1] == "UniProt":
                                        UniProt_val_tuple_in_list.append(
                                            (key_val_tuples_in_list_child[3][1], key_val_tuples_in_list_child[4][1], entity_chainID_list[1][1]))
                                        PDBe_val_tuples_in_list_for_Uni.append(
                                            (key_val_tuples_in_list_parent[2][1], key_val_tuples_in_list_parent[3][1], entity_chainID_list[1][1]))
                                        UniProtdbAccessionId_list.append(key_val_tuples_in_list_child[2][1])

    tuple_PDBe_for_PDB_and_tuple_PDB = list(zip(PDBe_val_tuples_in_list_for_PDB, PDB_val_tuples_in_list))
    tuple_PDBe_for_UniProt_and_tuple_UniProt = list(zip(PDBe_val_tuples_in_list_for_Uni, UniProt_val_tuple_in_list, UniProtdbAccessionId_list))

    return [tuple_PDBe_for_PDB_and_tuple_PDB, tuple_PDBe_for_UniProt_and_tuple_UniProt, UniProt_conversion_dict, details_list]


def SIFTS_data_parser_for_mmCIF(tuple_PDBe_for_PDB_and_tuple_PDB, tuple_PDBe_for_UniProt_and_tuple_UniProt, default_mmCIF_num,
                                chains_to_change="all"):
    df_PDBe_UniProt = pd.DataFrame(tuple_PDBe_for_UniProt_and_tuple_UniProt, columns=['PDBe', 'UniProt', "AccessionID"])
    df_PDBe_UniProt = df_PDBe_UniProt.drop_duplicates(subset="PDBe", keep='first')
    df_PDBe_PDB = pd.DataFrame(tuple_PDBe_for_PDB_and_tuple_PDB, columns=['PDBe', 'PDB'])
    df_PDBe_PDB = df_PDBe_PDB.drop_duplicates(subset="PDBe", keep='first')

    df_PDBe_PDB_UniProt = df_PDBe_PDB.merge(df_PDBe_UniProt, left_on="PDBe", right_on="PDBe", how='left')
    df_PDBe_PDB_UniProt['UniProt'] = df_PDBe_PDB_UniProt['UniProt'].replace(np.nan, "50000")
    df_PDBe_PDB_UniProt["Uni_moD"] = np.where(df_PDBe_PDB_UniProt['UniProt'] != "50000", df_PDBe_PDB_UniProt['UniProt'], df_PDBe_PDB_UniProt["PDBe"])
    df_PDBe_PDB_UniProt.loc[:, 'new_col_Uni'] = df_PDBe_PDB_UniProt.Uni_moD.map(lambda x: x[0])
    df_PDBe_PDB_UniProt["UniProt_50k"] = df_PDBe_PDB_UniProt.new_col_Uni.apply(lambda x: str(int(x) + default_mmCIF_num if type(x) == str else x))
    df_PDBe_PDB_UniProt.loc[df_PDBe_PDB_UniProt['UniProt'] != '50000', 'UniProt_50k'] = df_PDBe_PDB_UniProt['new_col_Uni']

    Three_Rows_CIF_Num_Uni = []
    if chains_to_change == "all":
        for index, rows in df_PDBe_PDB_UniProt.iterrows():
            intermediate_list = [rows.PDBe, rows.UniProt_50k, rows.Uni_moD, rows.PDB, rows.AccessionID]
            Three_Rows_CIF_Num_Uni.append(intermediate_list)

    else:
        for index, rows in df_PDBe_PDB_UniProt.iterrows():
            if rows.PDB[2].strip() in chains_to_change:
                intermediate_list = [rows.PDBe, rows.UniProt_50k, rows.Uni_moD, rows.PDB, rows.AccessionID]
            else:
                intermediate_list = [rows.PDBe, rows.PDB[0], rows.Uni_moD, rows.PDB, rows.AccessionID]
            Three_Rows_CIF_Num_Uni.append(intermediate_list)

    df_PDBe_PDB_UniProt["Three_Rows_CIF_Num_Uni"] = Three_Rows_CIF_Num_Uni
    df_PDBe_PDB_UniProt_without_null = df_PDBe_PDB_UniProt[df_PDBe_PDB_UniProt.PDB.map(lambda x: x[0]) != "null"]
    df_PDBe_PDB_UniProt_without_null_index_PDBe = df_PDBe_PDB_UniProt_without_null.set_index("PDBe")

    return [df_PDBe_PDB_UniProt_without_null_index_PDBe, df_PDBe_PDB_UniProt]


def master_mmCIF_renumber_function(input_mmCIF_file_were_found, default_input_path_to_mmCIF,
                                   default_input_path_to_SIFTS, default_output_path_to_mmCIF,
                                   default_mmCIF_num, gzip_mode, exception_AccessionIDs):
    input_mmCIF_assembly_files_were_found_list = list()
    input_mmCIF_assembly_files_were_found_list.append(input_mmCIF_file_were_found)

    for mmCIF_name in input_mmCIF_assembly_files_were_found_list:
        log_message = list()
        SIFTS_name = mmCIF_name[:4] + ".xml.gz"

        # for no SIFTS _no_SIFTS_out.cif.gz
        try:
            gzip.open(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name), 'rt')
        except FileNotFoundError:
            mmcif_dict = try_MMCIF2Dict(default_input_path_to_mmCIF, mmCIF_name)
            if mmcif_dict == 0:
                continue
            copy_file(default_input_path_to_mmCIF, mmCIF_name, default_output_path_to_mmCIF, ".cif.gz", gzip_mode)
            log_message = if_no_SIFTS_data_log(mmCIF_name, mmcif_dict, log_message)
            return log_message

        # for zerobyte SIFTS _zerobyte_SIFTS_out.cif.gz
        if os.path.getsize(Path(str(default_input_path_to_SIFTS) + "/" + SIFTS_name)) == 0:
            mmcif_dict = try_MMCIF2Dict(default_input_path_to_mmCIF, mmCIF_name)
            if mmcif_dict == 0:
                continue
            copy_file(default_input_path_to_mmCIF, mmCIF_name, default_output_path_to_mmCIF, ".cif.gz", gzip_mode)
            log_message = if_no_SIFTS_data_log(mmCIF_name, mmcif_dict, log_message)
            return log_message

        product_tree_SIFTS = try_SIFTS_tree_parser(default_input_path_to_SIFTS, SIFTS_name)
        if product_tree_SIFTS == 0:
            continue

        tuple_PDBe_for_PDB_and_tuple_PDB = product_tree_SIFTS[0]
        tuple_PDBe_for_UniProt_and_tuple_UniProt = product_tree_SIFTS[1]
        UniProt_conversion_dict = product_tree_SIFTS[2]

        # _no UniProt in SIFTS _no_UniProt_in_SIFTS_out.cif.gz
        if tuple_PDBe_for_UniProt_and_tuple_UniProt == list():
            mmcif_dict = try_MMCIF2Dict(default_input_path_to_mmCIF, mmCIF_name)
            if mmcif_dict == 0:
                continue
            copy_file(default_input_path_to_mmCIF, mmCIF_name, default_output_path_to_mmCIF, ".cif.gz", gzip_mode)
            log_message = if_no_SIFTS_data_log(mmCIF_name, mmcif_dict, log_message)
            return log_message

        product_of_SIFTS_data_parser = SIFTS_data_parser_for_mmCIF(tuple_PDBe_for_PDB_and_tuple_PDB, tuple_PDBe_for_UniProt_and_tuple_UniProt,
                                                                   default_mmCIF_num, 'all')
        df_PDBe_PDB_UniProt = product_of_SIFTS_data_parser[1]

        # all good till here
        handling_chain_numbering = handling_chain_numbering_clashes(df_PDBe_PDB_UniProt, exception_AccessionIDs)
        chains_to_change = handling_chain_numbering[0]
        combined_tuple_PDBe_UniProt_AccessionID = handling_chain_numbering[1]
        longest_AccessionID_list = handling_chain_numbering[3]
        chains_to_change_one_to_end = handling_chain_numbering[4]

        product_of_SIFTS_data_parser = SIFTS_data_parser_for_mmCIF(tuple_PDBe_for_PDB_and_tuple_PDB, combined_tuple_PDBe_UniProt_AccessionID,
                                                                   default_mmCIF_num, chains_to_change)
        df_PDBe_PDB_UniProt_without_null_index_PDBe = product_of_SIFTS_data_parser[0]
        df_PDBe_PDB_UniProt = product_of_SIFTS_data_parser[1]

        renumbered_count = renumbered_count_in_chains(chains_to_change_one_to_end, df_PDBe_PDB_UniProt_without_null_index_PDBe,
                                                      mmCIF_name, UniProt_conversion_dict, longest_AccessionID_list)
        chain_total_renum = renumbered_count[0]
        nothing_changed = renumbered_count[1]

        chain_total_renum.append(nothing_changed)
        mod_log_message = chain_total_renum

        # for no change needed _no_change_out.cif.gz
        if nothing_changed == 0:
            copy_file(default_input_path_to_mmCIF, mmCIF_name, default_output_path_to_mmCIF, ".cif.gz", gzip_mode)
            return mod_log_message

        product_of_mmCIF_parser = mmCIF_parser(mmCIF_name, default_input_path_to_mmCIF, df_PDBe_PDB_UniProt_without_null_index_PDBe,
                                               default_mmCIF_num, chains_to_change, chains_to_change_one_to_end)
        df_final_atom_site = product_of_mmCIF_parser[0]
        mmcif_dict = product_of_mmCIF_parser[1]

        poly_nonpoly_append = poly_nonpoly_renum(mmcif_dict, df_PDBe_PDB_UniProt, chains_to_change, default_mmCIF_num)
        poly_nonpoly_atom_site = poly_nonpoly_append.append(df_final_atom_site).drop_duplicates(subset="PDB_num_and_chain", keep='first')

        formed_columns = column_formation(mmcif_dict)
        renumber_tables(formed_columns, mmcif_dict, poly_nonpoly_atom_site, chains_to_change, default_mmCIF_num)

        try:
            output_with_this_name_ending("_renum.cif", default_output_path_to_mmCIF, mmcif_dict, mmCIF_name=mmCIF_name,
                                         gzip_mode=gzip_mode, current_directory=current_directory)
            return mod_log_message
        except IndexError:
            # 5olg data swapped columns
            print("IndexError Warning this file is not renumbered:", mmCIF_name)
            copy_file(default_input_path_to_mmCIF, mmCIF_name, default_output_path_to_mmCIF, ".cif.gz", gzip_mode)
