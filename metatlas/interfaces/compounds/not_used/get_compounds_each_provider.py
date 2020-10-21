import sys,os
import pandas as pd
#import get_compounds_from_wikidata as wd
sys.path.append('/global/project/projectdirs/openmsi/jupyterhub_libs/anaconda/lib/python2.7/site-packages')
from rdkit import Chem
import numpy as np
from rdkit.Chem import PandasTools

# MolToSmiles( (Mol)mol [, (bool)isomericSmiles=False 
# http://www.rdkit.org/Python_Docs/rdkit.Chem.rdmolfiles-module.html#MolToSmiles
#         - isomericSmiles: (optional) include information about stereochemistry in
#           the SMILES.  Defaults to false.

# https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/RZJQGNCSTQAWON-UHFFFAOYSA-N/synonyms/json


def desalt(mol):
    #input is an rdkit mol
    #returns an rdkit mol keeping the biggest component
    #returns original mol if only one component
    #returns a boolean indicated if cleaning was necessary
    d = Chem.rdmolops.GetMolFrags(mol) #these are atom indices
    if len(d) == 1: #If there are fragments or multiple molecules this will be greater than 1 
        return mol,False
    my_smiles=Chem.MolToSmiles(mol,True)
    parent_atom_count=0;
    disconnected=my_smiles.split('.')
    #With GetMolFrags, we've already established that there is more than one disconnected structure
    status = False
    for s in disconnected:
        little_mol=Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(s)))
        if little_mol is not None:
            count = little_mol.GetNumAtoms()
            if count > parent_atom_count:
                parent_atom_count = count
                parent_mol = little_mol
                status = True
    return parent_mol,status

def get_wikidata(terms_to_keep):
    prop_file = '/project/projectdirs/openmsi/projects/compound_data/wikidata/wikidata_compound_properties.xlsx'
    result = wd.get_wikidata(prop_file)
    df = pd.DataFrame(result)
    df.rename(columns={'canonicalSMILES': 'smiles'}, inplace=True)
    df.rename(columns={'InChI': 'inchi'}, inplace=True)
    df.rename(columns={'compoundLabel': 'common_name'}, inplace=True)
    df['source_database'] = 'wikidata'
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df

    
def get_img(terms_to_keep):
    df = pd.read_csv('/project/projectdirs/openmsi/projects/compound_data/img_abc/NPlist32763_30-mar-2016.xls',delimiter='\t')
    df.rename(columns={'SMILES': 'smiles'}, inplace=True)
    df.rename(columns={'InChl': 'inchi'}, inplace=True)
    df.rename(columns={'SM ID': 'img_abc_id'}, inplace=True)
    df.rename(columns={'Secondary Metabolite (SM) Name': 'common_name'}, inplace=True)
    df['source_database'] = 'img'
    df['ROMol'] = ''
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df
   
    
def get_enzo(terms_to_keep):
    df = pd.read_csv('/project/projectdirs/openmsi/projects/compound_data/enzo/BML-2865.txt',delimiter='\t')
    df.rename(columns={'SMILES': 'smiles'}, inplace=True)
    df.rename(columns={'Name': 'common_name'}, inplace=True)
    df['inchi'] = np.nan
    df['source_database'] = 'enzo'
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df


def get_msmls(terms_to_keep):
    df = pd.read_excel('/project/projectdirs/openmsi/projects/compound_data/msmls/MSMLS map - mz overlap edit sk.xlsx')
    df.rename(columns={'SMILES': 'smiles'}, inplace=True)
    df.rename(columns={'CNAME': 'common_name'}, inplace=True)
    df.rename(columns={'PC_CID': 'pubchem_compound_id'}, inplace=True)
    df['source_database'] = 'msmls'
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df

def dequote(s):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s

def get_metacyc(terms_to_keep):
    df = pd.read_excel('/project/projectdirs/openmsi/projects/compound_data/metacyc/MetAtlas_Export_MetaCyc_Compounds.xlsx', encoding='utf-8')
    df.rename(columns={'InChI': 'inchi'}, inplace=True)
    df.rename(columns={'KEGG': 'kegg_id'}, inplace=True)
    df.rename(columns={'PubChem': 'pubchem_compound_id'}, inplace=True)
    df.rename(columns={'Common-Name': 'common_name'}, inplace=True)
    df.rename(columns={'Names': 'synonyms'}, inplace=True) # reduced DCPIP // "2,6-dichloro-4-[(4-hydroxyphenyl)amino]phenol" // "reduced dichloroindophenol" // "reduced 2,6-dichlorophenolindophenol" // "reduced DCIP"
    df.rename(columns={'Object ID': 'metacyc_id'}, inplace=True)
    #df.synonyms.astype(str,inplace=True)
    #df['synonyms'].astype(basestring)
    df.loc[:,'synonyms'] = [[ s.strip() for s in mystr.split('//')] for mystr in df['synonyms'].astype(str).tolist() ]
    df.loc[:,'synonyms'] = [[ dequote(s) for s in mystr] for mystr in df['synonyms'].tolist() ]

    df['source_database'] = 'metacyc'
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df

# terms_to_keep = ['smiles','inchi','source_database','ROMol','common_name','synonyms','pubchem_compound_id','lipidmaps_id','metacyc_id','hmdb_id','img_abc_id','chebi_id','kegg_id']
def get_gnps(terms_to_keep):
    from pyteomics import mgf
    gnps = [s['params'] for s in mgf.read('/project/projectdirs/openmsi/projects/compound_data/gnps/ALL_GNPS (1).mgf')]
    df = pd.DataFrame(gnps)
    df['source_database'] = 'gnps'
#     df.rename(columns={'name': 'name'}, inplace=True)
#name has adduct "Hoiamide B M+H"
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df


def get_dr_dukes():
    df = pd.read_csv('/project/projectdirs/openmsi/projects/compound_data/dr_dukes_phytochemicals/CHEMICALS.csv',delimiter=',')
    df['source_database'] = 'dr_dukes'
    print df.keys()
    return df

def get_lipid_maps(terms_to_keep):
    df = PandasTools.LoadSDF('/project/projectdirs/openmsi/projects/compound_data/lipidmaps/LMSDFDownload28Jun15FinalAll.sdf')
    df['source_database'] = 'lipidmaps'
    df.rename(columns={'KEGG_ID': 'kegg_id'}, inplace=True)
    df.rename(columns={'PUBCHEM_CID': 'pubchem_compound_id'}, inplace=True)
    df.rename(columns={'COMMON_NAME': 'common_name'}, inplace=True)
    df.rename(columns={'SYNONYMS': 'synonyms'}, inplace=True)
#     Decanohydroxamic acid; caprinohydroxamic acid; n-Decanohydroxamic acid
    df.loc[:,'synonyms'] = [[ s.strip() for s in mystr.split(';')] for mystr in df['synonyms'].astype(str).tolist() ]
    df.rename(columns={'ID': 'lipidmaps_id'}, inplace=True) 
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df
    
    
def get_hmdb(terms_to_keep):
    df = PandasTools.LoadSDF('/project/projectdirs/openmsi/projects/compound_data/hmdb/structures.sdf')
    df['source_database'] = 'hmdb'
    df.rename(columns={'GENERIC_NAME': 'common_name'}, inplace=True)
    df.rename(columns={'SYNONYMS': 'synonyms'}, inplace=True)
    df.loc[:,'synonyms'] = [[ s.strip() for s in mystr.split(';')] for mystr in df['synonyms'].astype(str).tolist() ]

#     2-(8S,9S,13S,14S)-3-Hydroxy-2-methoxy-13-methyl-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]phenanthren-17-one; 2-Hydroxyestrone 2-methyl ether; 2-Methoxy-17-oxoestra-1,3,5(10)-trien-3-ol; 2-Methoxy-3-hydroxyestra-1,3,5(10)-trien-17-one; 3-Hydroxy-2-methoxy-Estra-1,3,5(10)-trien-17-one; 3-Hydroxy-2-methoxyestra-1,3,5(10)-trien-17-one; Methoxy-Estrone	
    df.rename(columns={'HMDB_ID': 'hmdb_id'}, inplace=True)
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df


def get_chembl(terms_to_keep):
    sdf_file = '/project/projectdirs/openmsi/projects/compound_data/chembl/chembl_21.sdf.gz'
    df = PandasTools.LoadSDF(sdf_file)
    df['source_database'] = 'chembl'
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df

def get_chebi(terms_to_keep):
    df = PandasTools.LoadSDF('/project/projectdirs/openmsi/projects/compound_data/chebi/ChEBI_complete.sdf.gz')
#     df = PandasTools.LoadSDF('/project/projectdirs/openmsi/projects/compound_data/chebi/ChEBI_complete_3star.sdf.gz')
    for index, row in df.iterrows():
        mol = row['ROMol']
        try:
            df.loc[index,'inchi'] = Chem.MolToInchi(mol)
        except:
            pass
    df['source_database'] = 'chebi'
    df.rename(columns={'KEGG COMPOUND Database Links': 'kegg_id'}, inplace=True)
    df.rename(columns={'ChEBI Name': 'common_name'}, inplace=True)
    df.rename(columns={'Synonyms': 'synonyms'}, inplace=True)
    df.loc[:,'synonyms'] = [[ s.strip() for s in mystr.split('\n')] for mystr in df['synonyms'].astype(str).tolist() ]

    # (-)-Epicatechin\n(-)-Epicatechol\n(2R,3R)-(-)-Epicatechin\n(2R,3R)-2-(3,4-dihydroxyphenyl)-3,4-dihydro-2H-1-benzopyran-3,5,7-triol\n3,3',4',5,7-Pentahydroxyflavane\nEpicatechol\nEpigallocatechin\nL(-)-Epicatechin\nL-Acacatechin\nL-Epicatechin\nL-Epicatechol\nalpha-Catechin
    df.rename(columns={'ChEBI ID': 'chebi_id'}, inplace=True)
    k = df.keys()
    for t in terms_to_keep:
        if not t in k:
            df[t] = ''
    return df


""" contribution from Hans de Winter """
def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y)) for x,y in patts]

_reactions=None
def NeutraliseCharges(mol, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
#     mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            # print Chem.MolToSmiles(mol,True)
            # print Chem.MolToSmiles(mol), Chem.MolToSmarts(reactant),Chem.MolToSmiles(product)
            rms = Chem.AllChem.ReplaceSubstructs(mol, reactant, product, replaceAll=True)
            # rms_smiles = Chem.MolToSmiles(rms[0],True)
            # mol = Chem.MolFromSmiles(rms_smiles)
            mol = rms[0]
            # print Chem.MolToSmiles(mol,True)
    if replaced:
        return (mol, True) #Chem.MolToSmiles(mol,True)
    else:
        return (mol, False)
    
# _reactions=None
# def NeutraliseCharges(smiles, reactions=None):
#     global _reactions
#     if reactions is None:
#         if _reactions is None:
#             _reactions=_InitialiseNeutralisationReactions()
#         reactions=_reactions
#     mol = Chem.MolFromSmiles(smiles)
#     replaced = False
#     for i,(reactant, product) in enumerate(reactions):
#         while mol.HasSubstructMatch(reactant):
#             replaced = True
#             rms = AllChem.ReplaceSubstructs(mol, reactant, product)
#             mol = rms[0]
#     if replaced:
#         return (Chem.MolToSmiles(mol,True), True)
#     else:
#         return (smiles, False)
    
