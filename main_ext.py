import pandas as pd
import numpy as np
import math
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

whole_data_df = pd.read_csv('./data/2SMILES_Ratio_Tg_TS_only2.csv')
copolymer_df = pd.DataFrame(whole_data_df, columns=['index', 'SMILES1', 'SMILES2', 'C1', 'C2', 'Tg[K]'])
copolymer_list = copolymer_df.values.tolist()

copolymer_parameter_list = []
for i in range(int(len(copolymer_list) / 2)):
    copolymer_parameter_list.append([[], 0, 0, 0]) #SMILES1, SMILES2, gradient, first(0.0), last(1.0)

    x1 = copolymer_list[i*2][3] / 100
    y1 = copolymer_list[i*2][5]
    x2 = copolymer_list[i*2+1][3] / 100
    y2 = copolymer_list[i*2+1][5]

    gradient = (y2 - y1) / (x2 - x1)
    bias = y1 - gradient * x1
    last = bias + gradient

    copolymer_parameter_list[i][0] = [copolymer_list[i*2][1], copolymer_list[i*2][2]]
    copolymer_parameter_list[i][1] = gradient
    copolymer_parameter_list[i][2] = bias
    copolymer_parameter_list[i][3] = last
NUMBER_OF_COPOLYMER = int(len(copolymer_parameter_list))

copolymer_min_max_list = []
for i in range(NUMBER_OF_COPOLYMER):
    copolymer_min_max_list.append([])
    copolymer_min_max_list[i].append(copolymer_parameter_list[i][2])
    copolymer_min_max_list[i].append(copolymer_parameter_list[i][3])
copolymer_min_max_array = np.array(copolymer_min_max_list)

print("feasible Tg[K] range: {:.2f} ~ {:.2f}".format(np.min(copolymer_min_max_array), np.max(copolymer_min_max_array)))
print("\nenter the target tg: ", end='')

target_Tg = float(input())
number_of_candidate = 0
list_of_candidate = []
for i in range(NUMBER_OF_COPOLYMER):
    copolymer_min = copolymer_min_max_array[i][0]
    copolymer_max = copolymer_min_max_array[i][1]
    if (copolymer_min <= target_Tg and target_Tg <= copolymer_max) or (copolymer_max <= target_Tg and target_Tg <= copolymer_min):
        number_of_candidate += 1
        list_of_candidate.append(copolymer_parameter_list[i])

print("there are {} possible candidate for target Tg[K] {}".format(number_of_candidate, target_Tg))

print("\nenter a motif: ", end='')
input_smi = input()
input_mol = Chem.MolFromSmiles(input_smi)
input_fps = FingerprintMols.GetRDKFingerprint(input_mol)

candidate_name = sum([candidate[0] for candidate in list_of_candidate], [])
candidate_name_to_mol = [Chem.MolFromSmiles(x) for x in candidate_name]
candidate_name_to_fps = [FingerprintMols.GetRDKFingerprint(x) for x in candidate_name_to_mol]

similarity_list = [DataStructs.FingerprintSimilarity(input_fps, x) for x in candidate_name_to_fps]
most_similar_molecule = candidate_name[similarity_list.index(max(similarity_list))]

print("most similar molecule in copolymer dataset: ", most_similar_molecule)
print("fingerprint similarity: ", max(similarity_list))

recommended_copolymer = list_of_candidate[ math.floor(similarity_list.index(max(similarity_list)) / 2) ]
recommended_SMILES1 = recommended_copolymer[0][0]
recommended_SMILES2 = recommended_copolymer[0][1]
recommended_ratio = (target_Tg - recommended_copolymer[2]) / recommended_copolymer[1]

print("\nrecommended polymers are as follows: ")
print("SMILES1: {}".format(recommended_SMILES1))
print("SMILES2: {}".format(recommended_SMILES2))
print("ratio: {:.2f}:{:.2f}".format(recommended_ratio * 100, (1-recommended_ratio) * 100))

# candidate of target -> 350.5
# [[['*OCCCOC(=O)c1ccc2c(c1)ccc(c2)C(=O)*', 'O(CCCOC(=O)c1sc(cc1)C(=O)*)*'], 46.60000000000002, 309.25, 355.85], [['*C(=O)c1cc(c(c(c1)OC)OC)c1c(c(cc(c1)C(=O)OCc1ccc(cc1)CO*)OC)OC', 'C(=O)(c1cc(c(c(c1)OC)OC)c1c(c(cc(c1)C(=O)OCCCCO*)OC)OC)*'], -106.25000000000001, 417.775, 311.525]]
