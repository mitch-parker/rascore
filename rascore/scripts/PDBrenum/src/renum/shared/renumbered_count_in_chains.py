def renumbered_count_in_chains(chains_to_change_one_to_end, df_PDBe_PDB_UniProt_without_null_index_PDBe, mmCIF_name,
                               UniProt_conversion_dict, longest_AccessionID_list):
    nothing_changed = 1
    chain_total_renum = list()
    renum_for_all_chains = 0
    total_renum5000 = 0
    chains_to_change = sorted(chains_to_change_one_to_end)
    chain_PDBe_PDB = dict()
    # prot_len = len(df_PDBe_PDB_UniProt_without_null_index_PDBe["Three_Rows_CIF_Num_Uni"])
    # UniProt_total_renum = list()

    for chain in chains_to_change:
        total_count_per_chain = 0
        renum_for_the_chains = 0
        renum5000 = 0
        UniProt_set = set()

        for PDBe_num_Uni_PDB in df_PDBe_PDB_UniProt_without_null_index_PDBe["Three_Rows_CIF_Num_Uni"]:
            if chain == PDBe_num_Uni_PDB[2][2]:
                chain_PDBe_PDB[chain] = PDBe_num_Uni_PDB[3][2]
                if type(PDBe_num_Uni_PDB[4]) != float:
                    UniProt_set.add(PDBe_num_Uni_PDB[4])
                total_count_per_chain += 1
                if int(PDBe_num_Uni_PDB[1]) > 5000:
                    renum5000 += 1
                    total_renum5000 += 1
                elif PDBe_num_Uni_PDB[1] != PDBe_num_Uni_PDB[3][0]:
                    renum_for_all_chains += 1
                    renum_for_the_chains += 1

        for accession in UniProt_set:
            renum_for_accession = 0
            count_accession_len = 0
            for PDBe_num_Uni_PDB in df_PDBe_PDB_UniProt_without_null_index_PDBe["Three_Rows_CIF_Num_Uni"]:
                if accession == PDBe_num_Uni_PDB[4]:
                    if chain == PDBe_num_Uni_PDB[2][2]:
                        count_accession_len += 1
                if chain == PDBe_num_Uni_PDB[2][2] and accession == PDBe_num_Uni_PDB[4]:
                    if PDBe_num_Uni_PDB[1] != PDBe_num_Uni_PDB[3][0]:
                        renum_for_accession += 1

            if len(longest_AccessionID_list) != 0:
                if accession in longest_AccessionID_list:
                    AccessionID_human_read_longest = UniProt_conversion_dict.get(accession)
                    chain_total_renum.append(
                        [mmCIF_name[:4] + "*", chain, chain_PDBe_PDB[chain], accession, AccessionID_human_read_longest, count_accession_len,
                         total_count_per_chain, renum_for_accession, renum5000])
                else:
                    AccessionID_human_read = UniProt_conversion_dict.get(accession)
                    chain_total_renum.append(
                        [mmCIF_name[:4], chain, chain_PDBe_PDB[chain], accession, AccessionID_human_read, count_accession_len, total_count_per_chain,
                         renum_for_accession, renum5000])
            else:
                if type(accession) != float:
                    AccessionID_human_read = UniProt_conversion_dict.get(accession)
                    chain_total_renum.append(
                        [mmCIF_name[:4], chain, chain_PDBe_PDB[chain], accession, AccessionID_human_read, count_accession_len, total_count_per_chain,
                         renum_for_accession, renum5000])

    if renum_for_all_chains == 0 and total_renum5000 == 0:
        nothing_changed = 0

    return [chain_total_renum, nothing_changed]
