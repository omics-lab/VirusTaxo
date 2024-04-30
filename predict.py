import pickle
from Bio import SeqIO
from collections import defaultdict
import numpy as np
from util import entropy
import argparse


def predict(fasta_file, model_path, entropy, enrichment):
    # print(f'Loading the model...')

    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    # print(f'Loading Done.')

    print('Id\tLength\tGenus\tEntropy')
    # k = params['k']
    # for key in model.keys(): ## so all keys are of same length
    #     k = len(key)
    #     break
    k = len(model.keys()[0])
    total_kmers = len(model.keys())

    # print('value of k', k)
    for record in SeqIO.parse(fasta_file, "fasta"):
        
        def get_genus(seq):
            count = defaultdict(int)
            # Iterate over k-mers in the sequence
            for idx in range(len(seq) - k + 1):
                kmer = seq[idx:idx + k]
                # Check if k-mer is in the model
                if kmer in model:
                    matched_genera = model[kmer]
                    # Found a match so increase the count for all the genus
                    for genus in matched_genera:
                        count[genus] += 1
            # each tuple is a genus and count
            genus_score_tuple = list(count.items())

            # No Genus was Classified
            if not genus_score_tuple:
                E = 1.0
                return -1, -1, E
            else:
                E = entropy(np.array([x for _, x in genus_score_tuple]))
                prediction = max(genus_score_tuple, key=lambda x: x[1])[0]
                cnt = max(genus_score_tuple, key=lambda x: x[1])[1]
                return prediction, cnt, E

        prediction_1, cnt_1, E_1 = get_genus(record.seq)
        prediction_2, cnt_2, E_2 = get_genus(record.seq.reverse_complement())
        enrichment_1 = cnt_1/total_kmers
        enrichment_2 = cnt_2/total_kmers
        length = len(record.seq)

        if prediction_1 == -1 and prediction_2 == -1:
            print(f'{record.id}\t{length}\tUnclassified')
            print(f'{record.id}\t{length}\tNo matches found in the database')
        else:
            if E_1 < entropy or E_2 < entropy:
                if cnt_1 > cnt_2:
                    print(f'{record.id}\t{length}\t{prediction_1}\t{E_1}')
                    if enrichment_1 > enrichment:
                        print(f'The sequence has an enrichment of {enrichment_1}, exceeding the enrichment cutoff {enrichment}')
                    else: 
                        print(f'The sequence has an enrichment of {enrichment_1}, below the enrichment cutoff {enrichment}')
                else:
                    print(f'{record.id}\t{length}\t{prediction_2}\t{E_2}')
                    if enrichment_2 > enrichment:
                        print(f'The sequence has an enrichment of {enrichment_1}, exceeding the enrichment cutoff {enrichment}')
                    else: 
                        print(f'The sequence has an enrichment of {enrichment_1}, below the enrichment cutoff {enrichment}')
                
            else:
                print(f'{record.id}\t{length}\tUnclassified\t')
                print('The entropy of the predictions not higher than the threshold, given sequence is Unclassified')           
                


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--model_path', required=True,
                        help='Absolute or relative path of pre-built model')

    parser.add_argument('--seq', required=True,
                        help='Absolute or relative path of fasta sequence file')
    
    parser.add_argument('--entropy', default=0.5,
                        help='Entropy value to be used as a threshold. If entropy > threshold, VirusTaxo will return unclassified')

    parser.add_argument('--enrichment', default = 0.5,
                        help='Enrichment score to be used as a threshold for unique kmer mapped')

    args = parser.parse_args()

    predict(args.seq, args.model_path, args.entropy, args.enrichment)
