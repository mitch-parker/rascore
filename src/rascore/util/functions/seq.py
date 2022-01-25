# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2022 Mitchell Isaac Parker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline

from .path import delete_path


def load_record_lst(path):

    return SeqIO.parse(path, "fasta")


def load_record_dict(path):

    return SeqIO.to_dict(load_record_lst(path))


def get_record_id(record):

    return str(record.id)


def get_record_seq(record):

    return str(record.seq)


def get_record_desc(record):

    return record.description


def build_record(seq, id, name=None, desc=None):

    return SeqRecord(Seq(seq), id=id, name=name, description=desc)


def clustal_omega(fasta_path, msa_path):

    delete_path(msa_path)

    cline = ClustalOmegaCommandline(infile=fasta_path, outfile=msa_path)
    cline()


def pair_seq_aln(
    ref_seq,
    mob_seq,
    **kwargs,
):

    matrix = kwargs.get("matrix", substitution_matrices.load("BLOSUM62"))
    gap_open = kwargs.get("gap_open", -10.0)
    gap_extend = kwargs.get("gap_extend", -0.5)

    return pairwise2.align.globalds(
        ref_seq, mob_seq, matrix, gap_open, gap_extend, penalize_end_gaps=(False, False)
    )


def calc_seq_id(ref_seq, mob_seq, gap=True, aln=True):

    if aln:
        aln = pair_seq_aln(ref_seq, mob_seq)

        best_aln = aln[0]

        ref_seq = best_aln[0]
        mob_seq = best_aln[1]

    seq_len = min(len(ref_seq), len(mob_seq))
    match_lst = [ref_seq[i] == mob_seq[i] for i in range(seq_len)]

    if gap:
        seq_id = (100 * sum(match_lst)) / seq_len
    else:
        gap_less_seq_len = sum(
            [1 for i in range(seq_len) if (ref_seq[i] != "-" and mob_seq[i] != "-")]
        )
        seq_id = (100 * sum(match_lst)) / gap_less_seq_len

    return seq_id
