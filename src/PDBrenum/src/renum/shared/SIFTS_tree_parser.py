from src.download.modules import *


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
                                            try:
                                                UniProt_conversion_dict[UniProt] = Human_readable
                                            except NameError:
                                                UniProt_conversion_dict["UniProt"] = "NameError"

                        for residue in listResidue:
                            key_val_tuples_in_list_parent = list(residue.attrib.items())
                            if key_val_tuples_in_list_parent[0][0] == "dbSource" and key_val_tuples_in_list_parent[0][1] == "PDBe":
                                PDBe_val_tuples_in_list.append((key_val_tuples_in_list_parent[2][1],
                                                                key_val_tuples_in_list_parent[3][1],
                                                                entity_chainID_list[1][1]))

                                for crossRefDb in residue:
                                    crossRefDb_list.append(crossRefDb.attrib)
                                    key_val_tuples_in_list_child = list(crossRefDb.attrib.items())

                                    if key_val_tuples_in_list_child[0][0] == "dbSource" and key_val_tuples_in_list_child[0][1] == "PDB":
                                        PDB_val_tuples_in_list.append((key_val_tuples_in_list_child[3][1],
                                                                       key_val_tuples_in_list_child[4][1],
                                                                       key_val_tuples_in_list_child[5][1]))
                                        PDBe_val_tuples_in_list_for_PDB.append((key_val_tuples_in_list_parent[2][1],
                                                                                key_val_tuples_in_list_parent[3][1],
                                                                                entity_chainID_list[1][1]))

                                    if key_val_tuples_in_list_child[0][0] == "dbSource" and key_val_tuples_in_list_child[0][1] == "UniProt":
                                        UniProt_val_tuple_in_list.append((key_val_tuples_in_list_child[3][1],
                                                                          key_val_tuples_in_list_child[4][1],
                                                                          entity_chainID_list[1][1]))
                                        PDBe_val_tuples_in_list_for_Uni.append((key_val_tuples_in_list_parent[2][1],
                                                                                key_val_tuples_in_list_parent[3][1],
                                                                                entity_chainID_list[1][1]))
                                        UniProtdbAccessionId_list.append(key_val_tuples_in_list_child[2][1])

    tuple_PDBe_for_PDB_and_tuple_PDB = list(zip(PDBe_val_tuples_in_list_for_PDB, PDB_val_tuples_in_list))
    tuple_PDBe_for_UniProt_and_tuple_UniProt = list(zip(PDBe_val_tuples_in_list_for_Uni, UniProt_val_tuple_in_list, UniProtdbAccessionId_list))

    return [tuple_PDBe_for_PDB_and_tuple_PDB, tuple_PDBe_for_UniProt_and_tuple_UniProt, UniProt_conversion_dict]
