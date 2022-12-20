from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd

# smiles example: 'C(CCOC(=O)C1CCC(CC1)C(=O)C)CO'

whole_data_df = pd.read_csv('./data/2SMILES_Ratio_Tg_TS_only2.csv')
copolymer_df = pd.DataFrame(whole_data_df, columns=['index', 'SMILES1', 'SMILES2', 'C1', 'C2', 'Tg[K]'])
copolymer_list = copolymer_df.values.tolist()

copolymer_name_1 = [x[1] for x in copolymer_list][0::2]
copolymer_name_2 = [x[2] for x in copolymer_list][0::2]
copolymer_name = sum([[x, y] for x, y in zip(copolymer_name_1, copolymer_name_2)], [])
copolymer_name_to_mol = [Chem.MolFromSmiles(x) for x in copolymer_name]
copolymer_name_to_fps = [FingerprintMols.GetRDKFingerprint(x) for x in copolymer_name_to_mol]

print("Enter a motif:", end='')
input_smi = input()
input_mol = Chem.MolFromSmiles(input_smi)
input_fps = FingerprintMols.GetRDKFingerprint(input_mol)

similarity_list = [DataStructs.FingerprintSimilarity(input_fps, x) for x in copolymer_name_to_fps]
most_similar_molecule = copolymer_name[similarity_list.index(max(similarity_list))]

print("Fingerprint similarity: ", max(similarity_list))
print("Most similar molecule in copolymer dataset: ", most_similar_molecule)