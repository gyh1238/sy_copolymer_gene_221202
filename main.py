import pandas as pd
import numpy as np
import random

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
print("enter the target tg: ", end='')

target_Tg = float(input())
number_of_candidate = 0
list_of_candidate = []
for i in range(NUMBER_OF_COPOLYMER):
    copolymer_min = copolymer_min_max_array[i][0]
    copolymer_max = copolymer_min_max_array[i][1]
    if (copolymer_min <= target_Tg and target_Tg <= copolymer_max) or (copolymer_max <= target_Tg and target_Tg <= copolymer_min):
        number_of_candidate += 1
        list_of_candidate.append(copolymer_parameter_list[i])

print("\nthere are {} possible candidate for target Tg[K] {}".format(number_of_candidate, target_Tg))

recommended_copolymer = random.choice(list_of_candidate)
recommended_SMILES1 = recommended_copolymer[0][0]
recommended_SMILES2 = recommended_copolymer[0][1]
recommended_ratio = (target_Tg - recommended_copolymer[2]) / recommended_copolymer[1]

print("\nrecommended polymers are as follows: ")
print("SMILES1: {}".format(recommended_SMILES1))
print("SMILES2: {}".format(recommended_SMILES2))
print("ratio: {:.2f}:{:.2f}".format(recommended_ratio * 100, (1-recommended_ratio) * 100))