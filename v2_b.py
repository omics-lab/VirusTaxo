# accession with complete species, genus and family are taken
# 12612 Accession, 8528 species, 1933 genus and 242 family
# ./Dataset/Accession_Species_Genus_Family_12612_meta.csv

import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
from tqdm import tqdm
from util import save_object


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', required=True,
                        help='File path of metadata')
    parser.add_argument('--seq', required=True,
                        help='Fasta sequence file path')
    parser.add_argument('--k', required=True, type=int,
                        help='Length of k-mer')
    parser.add_argument('--saving_dir', required=True,
                        help='Directory path to save model files')

    args = parser.parse_args()

    k = args.k
    df = pd.read_csv(args.meta)

    # Define the ranks to process
    ranks = ['Species', 'Genus', 'Family']

    for rank in ranks:
        print(f"Processing rank: {rank}")

        # Build accession-to-rank mapping for the current rank
        accession_to_rank = {row['Accession']: row[rank] for _, row in df.iterrows() if rank in row}

        db = defaultdict(set)

        # Parse sequences and build k-mer database
        for record in tqdm(SeqIO.parse(args.seq, format='fasta'), desc=f"Building database for {rank}"):
            if record.id not in accession_to_rank:
                continue
            seq = str(record.seq)
            rank_value = accession_to_rank[record.id]
            for idx in range(len(seq) - k + 1):
                db[seq[idx: idx + k]].add(rank_value)

        # Discard ambiguous k-mers with multiple rank annotations
        discard_candidate = set()
        for key, val in db.items():
            if len(val) > 1:
                discard_candidate.add(key)

        for key in discard_candidate:
            db.pop(key)

        # Save the database for the current rank
        saving_path = f"{args.saving_dir}/{rank}_database.pkl"
        save_object(db, saving_path)
        print(f"Database for {rank} saved at {saving_path}")

# import argparse
# import pandas as pd
# from Bio import SeqIO
# from collections import defaultdict
# from tqdm import tqdm
# from util import save_object


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--meta', required=True,
#                         help='File path of metadata')
#     parser.add_argument('--seq', required=True,
#                         help='Fasta sequence file path')
#     parser.add_argument('--k', required=True, type=int,
#                         help='Length of k-mer')
#     parser.add_argument('--saving_path', required=True,
#                         help='Absolute or relative path for saving model')

#     args = parser.parse_args()

#     k = args.k

#     df = pd.read_csv(args.meta)

#     accession_to_rank = {row['Accession']: row['rank'] for _, row in df.iterrows()}

#     db = defaultdict(set)

#     for record in tqdm(SeqIO.parse(args.seq, format='fasta')):
#         if record.id not in accession_to_genus:
#             continue
#         seq = str(record.seq)
#         rank = accession_to_rank[record.id]
#         for idx in range(len(seq) - k + 1):
#             db[seq[idx : idx + k]].add(rank)

#     # print('Length before discard:', len(db))
#     discard_candidate = set()

#     for key, val in db.items():
#         if len(val) > 1:
#             discard_candidate.add(key)

#     for key in discard_candidate:
#         db.pop(key)

#     save_object(db, args.saving_path)