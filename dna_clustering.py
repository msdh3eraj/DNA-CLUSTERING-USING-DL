#group: 14
import numpy as np
import csv
import random
from itertools import product
import matplotlib.pyplot as plt
def random_base():
    return random.choice(['A', 'C', 'G', 'T'])

def transition_mimic(seq, P_t):
    probabily=np.random.rand()
    if(probabily>P_t):
        return seq
    idx = random.randint(0, len(seq)-1)
    base = seq[idx]
    new_base = random.choice(['A', 'C', 'G', 'T'])
    while new_base == base or (base in ['A', 'G'] and new_base not in ['G', 'A']) or (base in ['C', 'T'] and new_base not in ['C', 'T']):
        new_base = random.choice(['A', 'C', 'G', 'T'])
    new_seq = list(seq)
    new_seq[idx] = new_base
    return ''.join(new_seq)

def transversion_mimic(seq, P_v):
    probabily=np.random.rand()
    # print(probabily)
    if(probabily>P_v):
        return seq
    idx = random.randint(0, len(seq)-1)
    base = seq[idx]
    if base in ['A', 'G']:
        new_base = random.choice(['C', 'T'])
    else:
        new_base = random.choice(['A', 'G'])
    new_seq = list(seq)
    new_seq[idx] = new_base
    return ''.join(new_seq)

def generate_mimics(seq, P_t, P_v):
    mimic1 = transition_mimic(seq, P_t)
    mimic2 = transversion_mimic(seq, P_v)
    m3=transition_mimic(seq,P_t)
    mimic3=transversion_mimic(m3,P_v)
    return mimic1, mimic2, mimic3

def kmer_count(seq, k):
    kmerDict = {}
    for k_mer in product('ACGT', repeat=k):
        kmer = ''.join(k_mer)
        kmerDict[kmer] = 0
    idx = 0
    while idx < len(seq) - k:
        try:
            kmerDict[seq[idx:idx + k]] += 1
        except KeyError:
            pass
        idx += 1
    return list(kmerDict.values())


def pos_gen(kmer):
    k = len(kmer)
    posx = 2 ** k
    posy = 2 ** k
    for i in range(1, k + 1):
        bp = kmer[-i]
        if bp == 'C':
            posx = posx - 2 ** (k - i)
            posy = posy - 2 ** (k - i)
        elif bp == 'A':
            posx = posx - 2 ** (k - i)
        elif bp == 'G':
            posy = posy - 2 ** (k - i)
    return int(posx - 1), int(posy - 1)


def cgr_gen(probs, k):
    kamers = product('ACGT', repeat=k)
    mat = np.zeros((2 ** k, 2 ** k))
    for i, kmer in enumerate(kamers):
        x, y = pos_gen(kmer)
        mat[y][x] = probs[i]
    return mat


k=6





P_t=0.01
P_v=0.001

with open('output.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        DNA=row[0]
        print(DNA)
        mimic1, mimic2, mimic3 = generate_mimics(DNA, P_t, P_v)
        kk=kmer_count(DNA,k)
        kk = normalize_array(kk)
        ll=cgr_gen(kk,k)
        plt.title("original DNA")
        plt.imshow(ll, cmap='hot')
        plt.colorbar()
        plt.show()
        kk=kmer_count(mimic1,k)
        kk = normalize_array(kk)
        ll=cgr_gen(kk,k)
        plt.title("mimic1 DNA")
        plt.imshow(ll, cmap='hot')
        plt.colorbar()
        plt.show()
        kk=kmer_count(mimic2,k)
        kk = normalize_array(kk)
        ll=cgr_gen(kk,k)
        plt.title("mimic2 DNA")
        plt.imshow(ll, cmap='hot')
        plt.colorbar()
        plt.show()
        kk=kmer_count(mimic3,k)
        kk = normalize_array(kk)
        ll=cgr_gen(kk,k)
        plt.title("mimic3 DNA")
        plt.imshow(ll, cmap='hot')
        plt.colorbar()
        plt.show()
