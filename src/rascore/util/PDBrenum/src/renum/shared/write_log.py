def log_writer(resulting):
    with open("log_corrected.tsv", "w") as f:
        comp_uni_human_uni_PDBid = list()
        pdb_id_set = set()
        formated_item = (
            format("SP", "<3")
            + format("PDB_id", "<7")
            + format("chain_PDB", "<12")
            + format("chain_auth", "<12")
            + format("UniProt", "<20")
            + format("SwissProt", "<20")
            + format("uni_len", ">10")
            + format("chain_len", ">10")
            + format("renum", ">10")
            + format("5k_or_50k", ">10")
        )
        f.write("%s\n" % formated_item)

        for n in resulting:
            if type(n) == list:
                for z in n:
                    if type(z) == list:
                        if z[0][-1] == "*":
                            try:
                                formated_item = (
                                    format("*", "<3")
                                    + format(z[0][:4], "<7")
                                    + format(z[1], "<12")
                                    + format(z[2], "<12")
                                    + format(z[3], "<20")
                                    + format(z[4], "<20")
                                    + format(z[5], ">10")
                                    + format(z[6], ">10")
                                    + format(z[7], ">10")
                                    + format(z[8], ">10")
                                )
                                pdb_id_set.add(z[0][:4])
                                comp_uni_human_uni_PDBid.append((z[3], z[4], z[0][:4]))
                            except Exception:
                                print(z)
                        else:
                            try:
                                formated_item = (
                                    format("+", "<3")
                                    + format(z[0], "<7")
                                    + format(z[1], "<12")
                                    + format(z[2], "<12")
                                    + format(z[3], "<20")
                                    + format(z[4], "<20")
                                    + format(z[5], ">10")
                                    + format(z[6], ">10")
                                    + format(z[7], ">10")
                                    + format(z[8], ">10")
                                )
                                pdb_id_set.add(z[0])
                                comp_uni_human_uni_PDBid.append((z[3], z[4], z[0][:4]))
                            except Exception:
                                print(z)
                        f.write("%s\n" % formated_item)

    uniq_comp_uni_human_uni_PDBid_translation = set()
    for n in comp_uni_human_uni_PDBid:
        if n[0] == "-":
            continue
        uniq_comp_uni_human_uni_PDBid_translation.add(n)

    with open("log_translator.tsv", "w") as file_handle:
        for listitem in uniq_comp_uni_human_uni_PDBid_translation:
            file_handle.write(
                listitem[0] + " " + listitem[1] + " " + listitem[2] + "\n"
            )
