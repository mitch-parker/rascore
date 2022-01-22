from src.download.modules import *


def handling_chain_numbering_clashes(df_PDBe_PDB_UniProt, exception_AccessionIDs):
    chains_to_change = set()
    chains_to_change_one_to_end = set()
    AccessionIDs = set()
    chain_AccessionID_dict = dict()

    for PDBe_num_UniProt_PDB_accession in df_PDBe_PDB_UniProt["Three_Rows_CIF_Num_Uni"]:
        if type(PDBe_num_UniProt_PDB_accession[4]) == float:
            continue
        chains_to_change.add(PDBe_num_UniProt_PDB_accession[3][2])
        chains_to_change_one_to_end.add(PDBe_num_UniProt_PDB_accession[2][2])
        AccessionIDs.add(PDBe_num_UniProt_PDB_accession[4])

    for chains in chains_to_change:
        accessions_in_chain = set()
        for PDBe_num_UniProt_PDB_accession in df_PDBe_PDB_UniProt["Three_Rows_CIF_Num_Uni"]:
            if chains == PDBe_num_UniProt_PDB_accession[3][2]:
                if PDBe_num_UniProt_PDB_accession[4] is not np.nan:
                    accessions_in_chain.add(PDBe_num_UniProt_PDB_accession[4])
        chain_AccessionID_dict[chains] = accessions_in_chain

    tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = list()
    combined_tuple_PDBe_UniProt_AccessionID = list()
    longest_AccessionID_list = list()
    clash = 0

    for chain_accession in chain_AccessionID_dict.items():
        chains_to_change_for_AccessionID = list()
        longest_AccessionID = None
        longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = list()
        if len(chain_accession[1]) > 1:
            for accessions in chain_accession[1]:
                tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = list()
                target_UniProt_numbers_in_chain = list()
                diff_another_UniProt_numbers_in_same_chain = list()
                diff_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = list()

                for PDBe_num_UniProt_PDB_accession in df_PDBe_PDB_UniProt["Three_Rows_CIF_Num_Uni"]:
                    if PDBe_num_UniProt_PDB_accession[4] == accessions and PDBe_num_UniProt_PDB_accession[3][2] == chain_accession[0] and \
                            PDBe_num_UniProt_PDB_accession[4] is not np.nan:
                        tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID.append(
                            (PDBe_num_UniProt_PDB_accession[0], PDBe_num_UniProt_PDB_accession[2], PDBe_num_UniProt_PDB_accession[4]))
                        target_UniProt_numbers_in_chain.append(PDBe_num_UniProt_PDB_accession[2])
                        chains_to_change_for_AccessionID.append(PDBe_num_UniProt_PDB_accession[3][2])
                    if PDBe_num_UniProt_PDB_accession[4] != accessions and PDBe_num_UniProt_PDB_accession[3][2] == chain_accession[0] and \
                            PDBe_num_UniProt_PDB_accession[4] is not np.nan:
                        diff_another_UniProt_numbers_in_same_chain.append(PDBe_num_UniProt_PDB_accession[2])
                        diff_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID.append(
                            (PDBe_num_UniProt_PDB_accession[0], PDBe_num_UniProt_PDB_accession[2], PDBe_num_UniProt_PDB_accession[4]))

                for target_Uni in target_UniProt_numbers_in_chain:
                    for diff_Uni in diff_another_UniProt_numbers_in_same_chain:
                        if target_Uni[0] == diff_Uni[0]:
                            clash = 1

                if accessions not in exception_AccessionIDs:
                    if len(longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID) < len(
                            tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID):
                        longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID
                        longest_AccessionID = accessions

                if longest_AccessionID is None:
                    if len(longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID) < len(
                            tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID):
                        longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID
                        longest_AccessionID = accessions

            if clash == 1:
                combined_tuple_PDBe_UniProt_AccessionID.extend(longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID)
                longest_AccessionID_list.append(longest_AccessionID)
            else:
                combined_tuple_PDBe_UniProt_AccessionID.extend(longest_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID)
                combined_tuple_PDBe_UniProt_AccessionID.extend(diff_tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID)
        else:
            for accessions in chain_accession[1]:
                tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID = list()
                target_UniProt_numbers_in_chain = list()

                for PDBe_num_UniProt_PDB_accession in df_PDBe_PDB_UniProt["Three_Rows_CIF_Num_Uni"]:
                    if PDBe_num_UniProt_PDB_accession[4] == accessions and PDBe_num_UniProt_PDB_accession[3][2] == chain_accession[0] and \
                            PDBe_num_UniProt_PDB_accession[4] is not np.nan:
                        tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID.append((PDBe_num_UniProt_PDB_accession[0],
                                                                                         PDBe_num_UniProt_PDB_accession[2],
                                                                                         PDBe_num_UniProt_PDB_accession[4]))
                        target_UniProt_numbers_in_chain.append(PDBe_num_UniProt_PDB_accession[2])
                        chains_to_change_for_AccessionID.append(PDBe_num_UniProt_PDB_accession[3][2])
            combined_tuple_PDBe_UniProt_AccessionID.extend(tuple_PDBe_for_UniProt_and_tuple_UniProt_for_AccessionID)

    return [chains_to_change, combined_tuple_PDBe_UniProt_AccessionID, AccessionIDs, longest_AccessionID_list, chains_to_change_one_to_end]
