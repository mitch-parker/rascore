# -*- coding: utf-8 -*-
"""
  Copyright 2022 Mitchell Isaac Parker

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""

y32_col = "Y32"
y71_col = "Y71"
sw1_col = "SW1"
sw2_col = "SW2"

id_col = "id"

core_path_col = "core_path"
rcsb_path_col = "rcsb_path"
renum_path_col = "renum_path"
rcsb_assembly_path_col = "rcsb_assembly_path"
renum_assembly_path_col = "renum_assembly_path"
sifts_path_col = "sifts_path"
edia_path_col = "edia_path"
interf_path_col = "interface_path"
pocket_path_col = "pocket_path"

modelid_col = "modelid"
chainid_col = "chainid"
resid_col = "resid"
resname_col = "resname"
atomid_col = "atomid"

pdb_id_col = "pdb_id"
pdb_code_col = "pdb_code"
swiss_id_col = "swiss_id"
uniprot_id_col = "uniprot_id"
pfam_col = "pfam"
prot_col = "prot_name"
method_col = "method"
resolution_col = "resolution"
r_factor_col = "r_factor"
seq_col = "sequence"
len_col = "length"
range_col = "range"
date_col = "deposit_data"
pmid_col = "pmid"

len_a_col = "length_a"
len_b_col = "length_b"
len_c_col = "length_c"
ang_a_col = "angle_alpha"
ang_b_col = "angle_beat"
ang_g_col = "angle_gamma"
space_col = "space_group"
cf_col = "crystal_form"

mut_status_col = "mutation_status"
mut_pos_col = "mutation_position"

interf_id_col = "interface_id"
interf_col = "interface"
interf_area_col = "interface_area"
cb_cont_col = "cb_contacts"
atomid_cont_col = "atomid_contacts"
total_cb_cont_col = "total_cb_contacts"
total_atomid_cont_col = "total_atomid_contacts"
interf_cont_col = "interface_contacts"
cb_dist_col = "cb_distance"
iso_col = "isologous"

pocket_id_col = "pocket_id"
pocket_col = "pocket"
pocket_site_col = "pocket_site"
pocket_volume_col = "pocket_volume"
pocket_score_col = "pocket_score"
pocket_status_col = "pocket_status"
pocket_type_col = "pocket_type"
pocket_cont_col = "pocket_contacts"

bio_lig_col = "biological_ligand"
ion_lig_col = "ion_ligand"
pharm_lig_col = "pharmacological_ligand"
chem_lig_col = "chemical_ligand"
mod_lig_col = "modification_ligand"
mem_lig_col = "membrane_ligand"
pocket_lig_col = "pocket_ligand"

gene_class_col = "gene_class"
nuc_class_col = "nucleotide_class"
mut_class_col = "mutation_class"
pharm_class_col = "pharmacological_class"
match_class_col = "match_class"
prot_class_col = "protein_class"
interf_class_col = "interface_class"
pocket_class_col = "pocket_class"

bound_lig_col = "bound_ligand"
bound_prot_col = "bound_protein"
bound_prot_swiss_id_col = "bound_protein_swiss_id"

bound_lig_cont_col = 'bound_ligand_contacts'
bound_prot_cont_col = 'bound_protein_contacts'

pharm_lig_site_col = "pharmacological_ligand_site"
pharm_lig_match_col = "pharmarmacological_ligand_match"
pharm_lig_smiles_col = "pharmacological_ligand_smiles"

bound_prot_pfam_col = "bound_protein_pfam"
bound_prot_site_col = "bound_protein_site"

bound_prot_chainid_col = "bound_protein_chainid"
bound_interf_chainid_col = "bound_interf_chainid"

edia_col = "edia_score"
b_factor_col = "b_factor"

phi_col = "phi"
psi_col = "psi"
omega_col = "omega"
chi1_col = "chi1"
altchi1_col = "altchi1"
chi2_col = "chi2"
altchi2_col = "altchi2"
chi3_col = "chi3"
chi4_col = "chi4"
chi5_col = "chi5"

complete_col = "complete"
rama_col = "rama"
rotamer_col = "rotamer"
bb_resid_col = "bb_resid"
sc_resid_col = "sc_resid"
bb_seq_col = "bb_seq"
sc_seq_col = "sc_seq"
bb_len_col = "bb_length"
sc_len_col = "sc_length"

cluster_col = "cluster"
nn_dist_col = "nn_dist"
constr_dist_col = "constraint_distance"

nn_cutoff_col = "nn_cutoff"
constr_cutoff_col = "constraint_cutoff"

total_col = "total"
total_chain_col = "total_chain"
total_entry_col = "total_entry"
total_cf_col = "total_cf"

loop_col = "loop"
silh_col = "silhouette_score"
simi_col = "similarity_score"
total_complete_col = "total_complete"
cluster_count_col = "cluster_count"
total_cluster_pre_col = "total_cluster_pre"
total_noise_pre_col = "total_noise_pre"
total_pruned_nn_col = "total_pruned_nn"
total_pruned_constr_col = "total_pruned_constr"
total_classified_nn_col = "total_classified_nn"
total_classified_constr_col = "total_classified_constr"
total_cluster_post_col = "total_cluster_post"
total_noise_post_col = "total_noise_post"
select_col = "select"

rep_col = "representative"
common_col = "common"
entropy_col = "entropy"
occup_col = "occupancy"

mean_col = "mean"
max_col = "max"

vect_1_col = "vector_1"
vect_2_col = "vector_2"
atom_dist_col = "atom_dist"
hb_status_col = "hydrogen_bond"
hb_angle_1_col = "hb_angle_1"
hb_angle_2_col = "hb_angle_2"
wmhb_angle_col = "wmhb_angle"
outlier_col = "outlier"

bond_col = "bond"
angle_col = "angle"

dih_col = "dihedral"
rmsd_col = "rmsd"

index_col = "index"
p_col = "p_val"
correct_p_col = "corrected_p_val"

corr_col = "correlation"

a_col = "a_val"
b_col = "b_val"
c_col = "c_val"
d_col = "d_val"
risk_ratio_col = "risk_ratio"
low_ci_col = "lower_ci"
up_ci_col = "upper_ci"
sig_col = "significance"

bb_col_lst = [phi_col, psi_col, omega_col]

sc_col_lst = [
    chi1_col,
    chi2_col,
    altchi1_col,
    altchi2_col,
    chi3_col,
    chi4_col,
    chi5_col,
]

dih_col_lst = bb_col_lst + sc_col_lst

reformat_col_lst = dih_col_lst

dist_col_lst = [
    nn_dist_col,
    constr_dist_col,
    vect_1_col,
    vect_2_col,
    atom_dist_col,
    hb_angle_1_col,
    hb_angle_2_col,
    wmhb_angle_col,
    outlier_col,
]

sum_col_lst = [common_col, entropy_col, occup_col]

data_col_lst = reformat_col_lst + dist_col_lst + sum_col_lst

class_col_lst = [
    gene_class_col,
    nuc_class_col,
    mut_class_col,
    pharm_class_col,
    match_class_col,
    prot_class_col,
    interf_class_col,
    pocket_class_col,
]

path_col_lst = [
    core_path_col,
    rcsb_path_col,
    renum_path_col,
    rcsb_assembly_path_col,
    renum_assembly_path_col,
    sifts_path_col,
    edia_path_col,
    interf_path_col,
    pocket_path_col,
]

count_col_dict = {
    pdb_id_col: total_chain_col,
    pdb_code_col: total_entry_col,
    cf_col: total_cf_col,
}

order_col_lst = [
    id_col,
    pdb_id_col,
    pdb_code_col,
    modelid_col,
    chainid_col,
    sw1_col,
    sw2_col,
    y32_col,
    y71_col,
    cf_col,
    cluster_col,
    total_col,
    total_chain_col,
    total_entry_col,
    total_cf_col,
    rep_col,
    rama_col,
    rotamer_col,
    nn_dist_col,
    constr_dist_col,
    interf_id_col,
    interf_col,
    pocket_id_col,
    pocket_col,
    prot_col,
    swiss_id_col,
    uniprot_id_col,
    pfam_col,
    mut_status_col,
    mut_pos_col,
    gene_class_col,
    nuc_class_col,
    mut_class_col,
    pharm_class_col,
    match_class_col,
    prot_class_col,
    interf_class_col,
    pocket_class_col,
    loop_col,
    silh_col,
    simi_col,
    total_complete_col,
    cluster_count_col,
    total_cluster_pre_col,
    total_noise_pre_col,
    total_pruned_nn_col,
    total_pruned_constr_col,
    total_classified_nn_col,
    total_classified_constr_col,
    total_cluster_post_col,
    total_noise_post_col,
    select_col,
    nn_cutoff_col,
    constr_cutoff_col,
    space_col,
    len_col,
    range_col,
    date_col,
    pmid_col,
    index_col,
    p_col,
    correct_p_col,
    corr_col,
    risk_ratio_col,
    low_ci_col,
    up_ci_col,
    a_col,
    b_col,
    c_col,
    d_col,
    resid_col,
    resname_col,
    atomid_col,
    interf_area_col,
    cb_cont_col,
    atomid_cont_col,
    total_cb_cont_col,
    total_atomid_cont_col,
    interf_cont_col,
    cb_dist_col,
    iso_col,
    pocket_site_col,
    pocket_volume_col,
    pocket_score_col,
    pocket_status_col,
    pocket_type_col,
    pocket_cont_col,
    phi_col,
    psi_col,
    omega_col,
    chi1_col,
    chi2_col,
    altchi1_col,
    altchi2_col,
    chi3_col,
    chi4_col,
    chi5_col,
    edia_col,
    b_factor_col,
    vect_1_col,
    vect_2_col,
    atom_dist_col,
    hb_status_col,
    hb_angle_1_col,
    hb_angle_2_col,
    wmhb_angle_col,
    outlier_col,
    bond_col,
    angle_col,
    rmsd_col,
    common_col,
    entropy_col,
    occup_col,
    bound_lig_col,
    bound_prot_col,
    bound_lig_cont_col,
    bound_prot_cont_col,
    bound_prot_swiss_id_col,
    pharm_lig_site_col,
    pharm_lig_match_col,
    pharm_lig_smiles_col,
    bound_prot_pfam_col,
    bound_prot_site_col,
    bound_prot_chainid_col,
    bound_interf_chainid_col,
    bio_lig_col,
    ion_lig_col,
    pharm_lig_col,
    chem_lig_col,
    mod_lig_col,
    mem_lig_col,
    pocket_lig_col,
    method_col,
    resolution_col,
    r_factor_col,
    seq_col,
    bb_resid_col,
    bb_seq_col,
    bb_len_col,
    sc_resid_col,
    sc_seq_col,
    sc_len_col,
    complete_col,
    len_a_col,
    len_b_col,
    len_c_col,
    ang_a_col,
    ang_b_col,
    ang_g_col,
    core_path_col,
    rcsb_path_col,
    renum_path_col,
    rcsb_assembly_path_col,
    renum_assembly_path_col,
    sifts_path_col,
    edia_path_col,
    interf_path_col,
    pocket_path_col,
]

rename_col_dict = {
    id_col: "ID",
    pdb_id_col: "PDB ID",
    pdb_code_col: "PDB Code",
    modelid_col: "Model",
    chainid_col: "Chain",
    sw1_col: "SW1 Conformation",
    sw2_col: "SW2 Conformation",
    y32_col: "Y32 Position",
    y71_col: "Y71 Position",
    hb_status_col:"H-Bond Substate",
    gene_class_col: "RAS Isoform",
    mut_status_col: "Mutation Status",
    nuc_class_col: "Nucleotide State",
    prot_class_col: "Bound Protein",
    pocket_class_col: "Inhibitor Site",
    match_class_col: "Inhibitor Chemistry",
    interf_class_col: "Homodimer Status",
    bound_prot_col: "Bound Protein Name",
    bound_prot_swiss_id_col: "Bound Protein SwissProt ID",
    bound_prot_chainid_col: "Bound Protein Chain",
    bio_lig_col: "Nucleotide",
    ion_lig_col: "Ion",
    pharm_lig_col: "Inhibitor",
    chem_lig_col: "Chemical",
    mod_lig_col: "Modification",
    mem_lig_col: "Membrane",
    pocket_lig_col: "Pocket",
    space_col: "Space Group",
    method_col: "Experiment Method",
    resolution_col: "Resolution",
    r_factor_col: "R-Factor",
    date_col: "Deposit Date",
    pocket_score_col: "Druggability Scores",
    pocket_volume_col: "Pocket Volume"
}


def get_dist_col(x_resid, y_resid, x_atomid=None, y_atomid=None, ext=None):

    if ext is None:
        ext = atom_dist_col

    dist_col = str(x_resid)

    if str(x_atomid) != "None":
        dist_col += f"({x_atomid})"

    dist_col += ":"
    dist_col += str(y_resid)

    if str(y_atomid) != "None":
        dist_col += f"({y_atomid})"

    if ext is not None:
        dist_col += f"_{ext}"

    return dist_col
