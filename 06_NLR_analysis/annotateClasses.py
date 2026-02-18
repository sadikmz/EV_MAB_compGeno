#!/usr/bin/env python3
"""
Classify NLR proteins into subclasses (CNL, TNL, NL, etc.) based on MEME
motif patterns from NLR-Parser output.

Adapted from the NLR-Parser output summarizer by Phillip Bayer:
  https://gist.github.com/philippbayer/0052f5ad56121cd2252a1c5b90154ed1

Output format (tab-separated):
  SeqName         NLR_class
  scaffold10886_nlr_1  CNL
  scaffold2236_nlr_1   CNL
  scaffold2236_nlr_2   CNL

Usage:
  annotateClasses.py <nlrparser_output.txt> [<nlrparser_output2.txt> ...]

A summary of classified and unclassified sequences is printed to stderr.
"""

import argparse
import sys

# Motif-to-domain mapping from:
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-75
MOTIFS_D = {
    1:  'nb_arc_cnl_or_tnl',
    2:  'nb_arc_cnl',
    3:  'nb_arc_cnl_or_tnl',
    4:  'nb_arc_cnl_or_tnl',
    5:  'nb_arc_cnl_or_tnl',
    6:  'nb_arc_cnl',
    7:  'linker_cnl_or_tnl',
    8:  'linker_cnl_or_tnl',
    9:  'lrr_cnl_or_tnl',
    10: 'nb_arc_cnl_or_tnl',
    11: 'lrr_cnl_or_tnl',
    12: 'nb_arc_cnl_or_tnl',
    13: 'tir_tnl',
    14: 'monocot',
    15: 'tir_tnl',
    16: 'prenb_cnl',
    17: 'prenb_cnl',
    18: 'tir_tnl',
    19: 'lrr_cnl_or_tnl',
    20: 'monocot',
}

# Domain-combination to NLR-class assignment
CLASS_DICT = {
    frozenset(['lrr_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']): 'TNL',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['lrr_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'TCNL',
    frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CNL',
    frozenset(['nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'monocot', 'prenb_cnl', 'nb_arc_cnl']): 'CNL',
    frozenset(['monocot', 'nb_arc_cnl_or_tnl']): 'N',
    frozenset(['linker_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CN',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']): 'TCNL',
    frozenset(['monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'CN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'TNL',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'NL',
    frozenset(['nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'CN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'prenb_cnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']): 'TN',
    frozenset(['nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl']): 'N',
    frozenset(['tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']): 'TCN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'NL',
    frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'CNL',
    frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'nb_arc_cnl_or_tnl']): 'TCN',
    frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'TCN',
    frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl']): 'TN',
    frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CNL',
    frozenset(['monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl']): 'TN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'tir_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']): 'TCNL',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['nb_arc_cnl_or_tnl', 'prenb_cnl']): 'CN',
    frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'TCNL',
    frozenset(['lrr_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl']): 'TNL',
    frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl']): 'TN',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'TNL',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'TCN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl']): 'N',
    frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl']): 'NL',
    frozenset(['lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl']): 'NL',
    frozenset(['monocot', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'TCN',
    frozenset(['linker_cnl_or_tnl', 'tir_tnl', 'nb_arc_cnl_or_tnl', 'lrr_cnl_or_tnl']): 'TNL',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'lrr_cnl_or_tnl']): 'CNL',
    frozenset(['nb_arc_cnl_or_tnl']): 'N',
    frozenset(['linker_cnl_or_tnl', 'monocot', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'CN',
    frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'TCNL',
    frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl']): 'CNL',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CN',
    frozenset(['monocot', 'nb_arc_cnl_or_tnl', 'prenb_cnl', 'nb_arc_cnl']): 'CN',
    frozenset(['tir_tnl', 'nb_arc_cnl_or_tnl', 'linker_cnl_or_tnl', 'lrr_cnl_or_tnl', 'monocot', 'nb_arc_cnl']): 'TCNL',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl']): 'N',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'CN',
    frozenset(['linker_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'prenb_cnl']): 'CN',
    frozenset(['lrr_cnl_or_tnl', 'nb_arc_cnl_or_tnl', 'nb_arc_cnl']): 'CNL',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Classify NLR proteins into subclasses based on MEME motif "
            "patterns from NLR-Parser output."
        )
    )
    parser.add_argument(
        "input_files",
        nargs="+",
        metavar="NLR_PARSER_OUTPUT",
        help="One or more NLR-Parser tab-delimited output files.",
    )
    return parser.parse_args()


def classify_sequence(seqname, motif_field):
    """
    Resolve the NLR class for a single sequence.

    Returns a tuple of (seqname, nlr_class) where nlr_class is 'UNKNOWN'
    when the motif combination or a motif number is not in the lookup tables.
    Warnings for unrecognised entries are written to stderr.
    """
    raw_motifs = motif_field.split(",")
    domains = set()

    for token in raw_motifs:
        token = token.strip()
        try:
            motif_num = int(token)
        except ValueError:
            print(
                f"WARNING: non-integer motif token '{token}' for {seqname} -- skipping token",
                file=sys.stderr,
            )
            continue

        if motif_num not in MOTIFS_D:
            print(
                f"WARNING: motif number {motif_num} not in motif table for {seqname} -- skipping motif",
                file=sys.stderr,
            )
            continue

        domains.add(MOTIFS_D[motif_num])

    domain_key = frozenset(domains)

    if domain_key not in CLASS_DICT:
        print(
            f"WARNING: unknown domain combination {sorted(domain_key)} for {seqname}",
            file=sys.stderr,
        )
        return seqname, "UNKNOWN"

    return seqname, CLASS_DICT[domain_key]


def process_file(filepath):
    """
    Yield (seqname, nlr_class) for every data line in an NLR-Parser output file.
    """
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            seqname = fields[1]
            motif_field = fields[-1]
            yield classify_sequence(seqname, motif_field)


def main():
    args = parse_args()

    print("SeqName\tNLR_class")

    total = 0
    unknown = 0

    for filepath in args.input_files:
        for seqname, nlr_class in process_file(filepath):
            print(f"{seqname}\t{nlr_class}")
            total += 1
            if nlr_class == "UNKNOWN":
                unknown += 1

    classified = total - unknown
    print(
        f"Summary: {total} sequences processed, "
        f"{classified} classified, "
        f"{unknown} UNKNOWN (novel motif combinations).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
