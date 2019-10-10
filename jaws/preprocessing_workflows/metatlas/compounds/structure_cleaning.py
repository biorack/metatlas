import sys
sys.path.append('/global/project/projectdirs/openmsi/jupyterhub_libs/anaconda/lib/python2.7/site-packages')
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

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

def desalt(mol):
    # This molecule escaped my patterns: InChI InChI=1S/2C6H11NO5.O.V/c2*1-3(5(8)9)7(12)4(2)6(10)11;;/h2*3-4,12H,1-2H3,(H,8,9)(H,10,11);;/q;;;+2/p-2/t2*3-,4-;;/m00../s1 gave an error Molecule must be fully connected by covalent bonds.:

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
        little_mol=Chem.MolFromInchi(Chem.MolToInchi(Chem.MolFromSmiles(s,sanitize=False)))
        #Sanitize=True will fail for choline sulfate.  Can't sanitize the radical.
        if little_mol is not None:
            count = little_mol.GetNumAtoms()
            if count > parent_atom_count:
                parent_atom_count = count
                parent_mol = little_mol
                status = True
    return parent_mol,status

def desalt_compounds_in_dataframe(x):
    '''
    df.ROMol = df.ROMol.apply(desalt_compounds_in_dataframe)
    '''
    if x:
        if x.GetNumAtoms()>1:
            c = desalt(x)
            if c[1]:
                return c[0]
            else:
                return x
    

def neutralize_compounds_in_dataframe(x):
    '''
    df.ROMol = df.ROMol.apply(neutralize_compounds_in_dataframe)    
    '''
    if x:
        if x.GetNumAtoms()> 0:
            neutral_mol = []
            try:
                c = NeutraliseCharges(x)
                neutral_mol = c[0]
            except:
                pass
            if neutral_mol:
                return neutral_mol
            
def calculate_num_radicals_in_dataframe(x):
    num_radicals = 0.0
    if x:
        num_radicals = Descriptors.NumRadicalElectrons(x)        
    return num_radicals

def calculate_formula_in_dataframe(x):
    formula = ''
    if x:
        formula = rdMolDescriptors.CalcMolFormula(x)
    return formula

def calculate_monoisotopic_mw_in_dataframe(x):
    mw = 0.0
    if x:
        mw = Descriptors.ExactMolWt(x)
    return mw

def calculate_inchi_in_dataframe(x):
    inchi = ''
    if x:
        try:
            inchi = Chem.MolToInchi(x)
        except:
            pass#This fails when can't kekulize mol
    return inchi

def calculate_flattened_inchi_in_dataframe(x):
    flattened_inchi = ''
    if x:
        sm = Chem.MolToSmiles(x).replace('@','')
        flattened_rdkit_mol = Chem.MolFromSmiles(sm)
        try:
            flattened_inchi = Chem.MolToInchi(flattened_rdkit_mol)
        except:
            pass#This fails when can't kekulize mol
    return flattened_inchi

def calculate_inchikey_in_dataframe(x):
    ik = ''
    if x:
        try:
            ik = Chem.InchiToInchiKey(x)
        except:
            pass#This fails when can't kekulize mol.  Carbo-cations are the culprit usually. 
    return ik

def calculate_charge_in_dataframe(x):
    if x:
        my_charge = Chem.GetFormalCharge(x)
        return my_charge