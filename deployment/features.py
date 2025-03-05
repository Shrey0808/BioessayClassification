from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import numpy as np
import pandas as pd

# --------------------------
# Define atom-type predicates
# --------------------------
def is_negative(atom):
    return atom.GetFormalCharge() < 0

def is_positive(atom):
    return atom.GetFormalCharge() > 0

def is_hbond_donor(atom):
    if atom.GetAtomicNum() in [7, 8]:
        return any(neigh.GetAtomicNum() == 1 for neigh in atom.GetNeighbors())
    return False

def is_hbond_acceptor(atom):
    if atom.GetAtomicNum() in [7, 8]:
        return atom.GetFormalCharge() <= 0
    return False

def is_aromatic(atom):
    return atom.GetIsAromatic()

def is_hydrophobic(atom):
    return (atom.GetAtomicNum() == 6) and (not atom.GetIsAromatic())

# --------------------------
# Helper function: Count pairs of atoms at a given bond distance
# --------------------------
def count_pairs(mol, pred1, pred2, bond_distance):
    count = 0
    num_atoms = mol.GetNumAtoms()
    for i in range(num_atoms):
        a1 = mol.GetAtomWithIdx(i)
        if not pred1(a1):
            continue
        for j in range(i+1, num_atoms):
            a2 = mol.GetAtomWithIdx(j)
            if not pred2(a2):
                continue
            path = Chem.GetShortestPath(mol, i, j)
            if len(path) - 1 == bond_distance:
                count += 1
    return count

# --------------------------
# Calculate features from SMILES
# --------------------------
def calculate_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None

    mol = Chem.AddHs(mol)
    features = {}

    # Standard descriptors
    features.update({
        'XLogP': Crippen.MolLogP(mol),
        'PSA': rdMolDescriptors.CalcTPSA(mol),
        'NumRot': Lipinski.NumRotatableBonds(mol),
        'NumHBA': Lipinski.NumHAcceptors(mol),
        'NumHBD': Lipinski.NumHDonors(mol),
        'MW': Descriptors.MolWt(mol)
    })

    # Pharmacophore-Based Descriptors
    for d in range(1, 8):
        features[f'NEG_{d:02d}_NEG'] = count_pairs(mol, is_negative, is_negative, d)
    for d in range(3, 8):
        features[f'NEG_{d:02d}_POS'] = count_pairs(mol, is_negative, is_positive, d)
    for d in range(1, 8):
        features[f'NEG_{d:02d}_HBD'] = count_pairs(mol, is_negative, is_hbond_donor, d)
    for d in range(3, 8):
        features[f'NEG_{d:02d}_HBA'] = count_pairs(mol, is_negative, is_hbond_acceptor, d)
    for d in range(2, 8):
        features[f'NEG_{d:02d}_ARC'] = count_pairs(mol, is_negative, is_aromatic, d)
    for d in range(2, 8):
        features[f'NEG_{d:02d}_HYP'] = count_pairs(mol, is_negative, is_hydrophobic, d)
    for d in range(3, 8):
        features[f'POS_{d:02d}_POS'] = count_pairs(mol, is_positive, is_positive, d)
    for d in range(2, 8):
        features[f'POS_{d:02d}_HBD'] = count_pairs(mol, is_positive, is_hbond_donor, d)
    for d in range(3, 8):
        features[f'POS_{d:02d}_HBA'] = count_pairs(mol, is_positive, is_hbond_acceptor, d)
    for d in range(2, 8):
        features[f'POS_{d:02d}_ARC'] = count_pairs(mol, is_positive, is_aromatic, d)
    for d in range(2, 8):
        features[f'POS_{d:02d}_HYP'] = count_pairs(mol, is_positive, is_hydrophobic, d)
    for d in range(3, 8):
        features[f'HBD_{d:02d}_HBD'] = count_pairs(mol, is_hbond_donor, is_hbond_donor, d)
    for d in range(3, 8):
        features[f'HBD_{d:02d}_HBA'] = count_pairs(mol, is_hbond_donor, is_hbond_acceptor, d)
    for d in range(2, 8):
        features[f'HBD_{d:02d}_ARC'] = count_pairs(mol, is_hbond_donor, is_aromatic, d)
    for d in range(2, 8):
        features[f'HBD_{d:02d}_HYP'] = count_pairs(mol, is_hbond_donor, is_hydrophobic, d)
    for d in range(3, 8):
        features[f'HBA_{d:02d}_HBA'] = count_pairs(mol, is_hbond_acceptor, is_hbond_acceptor, d)
    for d in range(3, 8):
        features[f'HBA_{d:02d}_ARC'] = count_pairs(mol, is_hbond_acceptor, is_aromatic, d)
    for d in range(2, 8):
        features[f'HBA_{d:02d}_HYP'] = count_pairs(mol, is_hbond_acceptor, is_hydrophobic, d)
    for d in range(1, 8):
        features[f'ARC_{d:02d}_ARC'] = count_pairs(mol, is_aromatic, is_aromatic, d)
    for d in range(2, 8):
        features[f'ARC_{d:02d}_HYP'] = count_pairs(mol, is_aromatic, is_hydrophobic, d)
    for d in range(1, 8):
        features[f'HYP_{d:02d}_HYP'] = count_pairs(mol, is_hydrophobic, is_hydrophobic, d)

    # Weighted Burden Number Features
    def get_bond_properties():
        bond_data = []
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            en1 = ELECTRONEGATIVITY.get(a1.GetAtomicNum(), 0)
            en2 = ELECTRONEGATIVITY.get(a2.GetAtomicNum(), 0)
            props = {
                'order': bond.GetBondTypeAsDouble(),
                'en_diff': abs(en1 - en2),
                'conjugated': int(bond.GetIsConjugated()),
                'aromatic': int(bond.GetIsAromatic())
            }
            bond_data.append(props)
        return bond_data

    ELECTRONEGATIVITY = {
        1: 2.20, 6: 2.55, 7: 3.04, 8: 3.44, 9: 3.98, 15: 2.19, 16: 2.58, 17: 3.16
    }
    
    bond_props = get_bond_properties()
    gc_values = [b['order'] for b in bond_props] if bond_props else [0]
    en_values = [b['en_diff'] for b in bond_props] if bond_props else [0]
    lp_values = [b['conjugated'] + b['aromatic'] for b in bond_props] if bond_props else [0]
    
    for q in [0.25, 0.50, 0.75, 1.00]:
        features[f'WBN_GC_L_{q:.2f}'] = np.quantile(gc_values, max(q-0.25, 0))
        features[f'WBN_GC_H_{q:.2f}'] = np.quantile(gc_values, q)
        features[f'WBN_EN_L_{q:.2f}'] = np.quantile(en_values, max(q-0.25, 0))
        features[f'WBN_EN_H_{q:.2f}'] = np.quantile(en_values, q)
        features[f'WBN_LP_L_{q:.2f}'] = np.quantile(lp_values, max(q-0.25, 0))
        features[f'WBN_LP_H_{q:.2f}'] = np.quantile(lp_values, q)

    # Special Features
    bad_groups = ['[Cl,Br,I]', '[S;D2](-O)-O', '[N+](=O)[O-]', '[B-](F)(F)F', '[PH4+]']
    def count_bad(smarts):
        pat = Chem.MolFromSmarts(smarts)
        return len(mol.GetSubstructMatches(pat)) if pat else 0
    features['BadGroup'] = sum(count_bad(p) for p in bad_groups)
    features['BBB'] = int((features['PSA'] < 90) and (features['MW'] < 500))

    columns = [
        'NEG_01_NEG', 'NEG_02_NEG', 'NEG_03_NEG', 'NEG_04_NEG', 'NEG_05_NEG', 'NEG_06_NEG', 'NEG_07_NEG',
        'NEG_03_POS', 'NEG_04_POS', 'NEG_05_POS', 'NEG_06_POS', 'NEG_07_POS',
        'NEG_01_HBD', 'NEG_02_HBD', 'NEG_03_HBD', 'NEG_04_HBD', 'NEG_05_HBD', 'NEG_06_HBD', 'NEG_07_HBD',
        'NEG_03_HBA', 'NEG_04_HBA', 'NEG_05_HBA', 'NEG_06_HBA', 'NEG_07_HBA',
        'NEG_02_ARC', 'NEG_03_ARC', 'NEG_04_ARC', 'NEG_05_ARC', 'NEG_06_ARC', 'NEG_07_ARC',
        'NEG_02_HYP', 'NEG_03_HYP', 'NEG_04_HYP', 'NEG_05_HYP', 'NEG_06_HYP', 'NEG_07_HYP',
        'POS_03_POS', 'POS_04_POS', 'POS_05_POS', 'POS_06_POS', 'POS_07_POS',
        'POS_02_HBD', 'POS_03_HBD', 'POS_04_HBD', 'POS_05_HBD', 'POS_06_HBD', 'POS_07_HBD',
        'POS_03_HBA', 'POS_04_HBA', 'POS_05_HBA', 'POS_06_HBA', 'POS_07_HBA',
        'POS_02_ARC', 'POS_03_ARC', 'POS_04_ARC', 'POS_05_ARC', 'POS_06_ARC', 'POS_07_ARC',
        'POS_02_HYP', 'POS_03_HYP', 'POS_04_HYP', 'POS_05_HYP', 'POS_06_HYP', 'POS_07_HYP',
        'HBD_03_HBD', 'HBD_04_HBD', 'HBD_05_HBD', 'HBD_06_HBD', 'HBD_07_HBD',
        'HBD_03_HBA', 'HBD_04_HBA', 'HBD_05_HBA', 'HBD_06_HBA', 'HBD_07_HBA',
        'HBD_02_ARC', 'HBD_03_ARC', 'HBD_04_ARC', 'HBD_05_ARC', 'HBD_06_ARC', 'HBD_07_ARC',
        'HBD_02_HYP', 'HBD_03_HYP', 'HBD_04_HYP', 'HBD_05_HYP', 'HBD_06_HYP', 'HBD_07_HYP',
        'HBA_03_HBA', 'HBA_04_HBA', 'HBA_05_HBA', 'HBA_06_HBA', 'HBA_07_HBA',
        'HBA_03_ARC', 'HBA_04_ARC', 'HBA_05_ARC', 'HBA_06_ARC', 'HBA_07_ARC',
        'HBA_02_HYP', 'HBA_03_HYP', 'HBA_04_HYP', 'HBA_05_HYP', 'HBA_06_HYP', 'HBA_07_HYP',
        'ARC_01_ARC', 'ARC_02_ARC', 'ARC_03_ARC', 'ARC_04_ARC', 'ARC_05_ARC', 'ARC_06_ARC', 'ARC_07_ARC',
        'ARC_02_HYP', 'ARC_03_HYP', 'ARC_04_HYP', 'ARC_05_HYP', 'ARC_06_HYP', 'ARC_07_HYP',
        'HYP_01_HYP', 'HYP_02_HYP', 'HYP_03_HYP', 'HYP_04_HYP', 'HYP_05_HYP', 'HYP_06_HYP', 'HYP_07_HYP',
        'WBN_GC_L_0.25', 'WBN_GC_H_0.25', 'WBN_GC_L_0.50', 'WBN_GC_H_0.50', 'WBN_GC_L_0.75', 'WBN_GC_H_0.75',
        'WBN_GC_L_1.00', 'WBN_GC_H_1.00',
        'WBN_EN_L_0.25', 'WBN_EN_H_0.25', 'WBN_EN_L_0.50', 'WBN_EN_H_0.50', 'WBN_EN_L_0.75', 'WBN_EN_H_0.75',
        'WBN_EN_L_1.00', 'WBN_EN_H_1.00',
        'WBN_LP_L_0.25', 'WBN_LP_H_0.25', 'WBN_LP_L_0.50', 'WBN_LP_H_0.50', 'WBN_LP_L_0.75', 'WBN_LP_H_0.75',
        'WBN_LP_L_1.00', 'WBN_LP_H_1.00',
        'XLogP', 'PSA', 'NumRot', 'NumHBA', 'NumHBD', 'MW', 'BBB', 'BadGroup'
    ]

    df = pd.DataFrame([features], columns=columns)
    df = df.fillna(0).replace([np.inf, -np.inf], 0)
    return df

if __name__ =="__main__":
    example_smiles = "CN1CCC2=C(C1)[C@H](C3=CC(=CC=C3F)OC2)OCO"
    features = calculate_features(example_smiles)
    print(features)