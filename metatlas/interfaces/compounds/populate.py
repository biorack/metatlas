"""Populate compound fields"""

import time
from typing import List, Optional
from urllib.parse import quote

import pandas as pd
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from metatlas.datastructures import metatlas_objects as metob
from metatlas.plots import dill2plots as dp
from metatlas.tools import cheminfo


def generate_template_atlas(
    raw_file_name: str, confidence_levels: List[str], polarity: str, name: str, mz_tolerance: float = 10
) -> metob.Atlas:
    data = pd.read_csv(raw_file_name, sep="\t")
    acceptable = data[data["confidence_category"].isin(confidence_levels)]
    by_polarity = acceptable[acceptable["polarity"] == polarity]
    atlas = dp.make_atlas_from_spreadsheet(
        by_polarity, name, filetype="dataframe", polarity=polarity, store=False, mz_tolerance=mz_tolerance
    )
    inchi_keys = [cid.compound[0].inchi_key for cid in atlas.compound_identifications]
    pubchem_results = query_pubchem(inchi_keys)
    for cid in atlas.compound_identifications:
        cid.compound[0] = fill_fields(cid.compound[0], pubchem_results)
    return atlas


def count_non_empty(compound: metob.Compound) -> int:
    return sum([v != "" for v in compound.trait_values().values()])


def flatten_inchi(mol: Chem.rdchem.Mol) -> str:
    smiles = Chem.MolToSmiles(mol).replace("@", "")
    flattened_rdkit_mol = Chem.MolFromSmiles(smiles)
    try:
        return Chem.MolToInchi(flattened_rdkit_mol)
    except Exception:  # This fails when can't kekulize mol
        return ""


def chunks(data: list, num: int):
    """Yield successive num-sized chunks from data."""
    for i in range(0, len(data), num):
        yield data[i : i + num]


def query_pubchem(inchi_key_list: List[str], items_per_query: int = 50):
    out = []
    counter = 1
    for inchi_key_sub_list in chunks(inchi_key_list, items_per_query):
        out.extend(pcp.get_compounds(inchi_key_sub_list, "inchikey"))
        if counter % 5 == 0:
            time.sleep(1)
        counter = counter + 1
    return out


def get_pubchem_compound(inchi_key: str, pub_chem_results: List[pcp.Compound]) -> Optional[pcp.Compound]:
    for compound in pub_chem_results:
        if compound.inchikey == inchi_key:
            return compound
    return None


def convert_id(input_id_type: str, output_id_type: str, query: str) -> str:
    base_url = "https://cts.fiehnlab.ucdavis.edu/rest/convert/"
    url = f"{base_url}{quote(input_id_type)}/{quote(output_id_type)}/{quote(query)}"
    result = requests.get(url)
    result.raise_for_status()
    return result.json()[0]["results"][0]


def set_id(rec: metob.Compound, metatlas_name: str, cts_name: str, base_url: str, inchi_key: str) -> None:
    id_attr_name = f"{metatlas_name}_id"
    url_attr_name = f"{metatlas_name}_url"
    try:
        rec.set_trait(
            id_attr_name, rec.trait_values()[id_attr_name] or convert_id("InChIKey", cts_name, inchi_key)
        )
        rec.set_trait(
            url_attr_name,
            rec.trait_values()[url_attr_name] or f"{base_url}{rec.trait_values()[id_attr_name]}",
        )
    except IndexError:
        pass


def set_all_ids(compound: metob.Compound, inchi_key: str) -> metob.Compound:
    ids = [
        {
            "metatlas_name": "hmdb",
            "cts_name": "Human Metabolome Database",
            "base_url": "https://www.hmdb.ca/metabolites/",
        },
        {
            "metatlas_name": "chebi",
            "cts_name": "ChEBI",
            "base_url": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId=",
        },
        {
            "metatlas_name": "lipidmaps",
            "cts_name": "LipidMAPS",
            "base_url": "https://www.lipidmaps.org/databases/lmsd/",
        },
        {
            "metatlas_name": "kegg",
            "cts_name": "KEGG",
            "base_url": "https://www.genome.jp/dbget-bin/www_bget?",
        },
    ]
    for id_type in ids:
        set_id(compound, id_type['metatlas_name'], id_type['cts_name'], id_type['base_url'], inchi_key)
    return compound


# pylint: disable=invalid-name
def fill_fields(c: metob.Compound, pubchem_results: List[pcp.Compound]) -> metob.Compound:
    """
    Populate bank fields that can be infered from other fields.
    Does not overwrite any existing values that are not None or ''"""
    mol = Chem.inchi.MolFromInchi(c.inchi)
    if mol is None:
        return c
    if c.neutralized_inchi:
        norm_mol = Chem.inchi.MolFromInchi(c.neutralized_inchi)
    else:
        norm_mol = cheminfo.normalize_molecule(mol)
    c.formula = c.formula or Chem.rdMolDescriptors.CalcMolFormula(mol)
    c.mono_isotopic_molecular_weight = c.mono_isotopic_molecular_weight or ExactMolWt(mol)
    c.permanent_charge = c.permanent_charge or Chem.GetFormalCharge(mol)
    c.number_components = c.number_components or 1  # type: ignore
    c.num_free_radicals = c.num_free_radicals or Chem.Descriptors.NumRadicalElectrons(mol)
    c.inchi_key = c.inchi_key or Chem.inchi.InchiToInchiKey(c.inchi)
    if not c.neutralized_inchi:
        norm_mol = cheminfo.normalize_molecule(mol)
        c.neutralized_inchi = Chem.inchi.MolToInchi(norm_mol)
    c.neutralized_inchi_key = c.neutralized_inchi_key or Chem.inchi.InchiToInchiKey(c.neutralized_inchi)
    if not c.neutralized_2d_inchi:
        if not norm_mol:
            norm_mol = Chem.inchi.MolFromInchi(c.neutralized_inchi)
        c.neutralized_2d_inchi = metob.MetUnicode(flatten_inchi(norm_mol))
    c.neutralized_2d_inchi_key = c.neutralized_2d_inchi_key or Chem.inchi.InchiToInchiKey(
        c.neutralized_2d_inchi
    )
    pubchem = get_pubchem_compound(c.inchi_key, pubchem_results)
    if pubchem is not None:
        c.pubchem_compound_id = c.pubchem_compound_id or pubchem.cid
        c.pubchem_url = c.pubchem_url or metob.MetUnicode(
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{c.pubchem_compound_id}"
        )
        c.synonyms = c.synonyms or metob.MetUnicode("///".join(pubchem.synonyms))
        c.iupac_name = c.iupac_name or pubchem.iupac_name
    if c.name in ["", "Untitled"]:
        c.name = c.synonyms.split("///")[0] or c.iupac_name
    return set_all_ids(c, c.inchi_key)


def create_c18_template_atlases() -> None:
    c18_data = "/global/u2/w/wjholtz/c18_atlas_creation.tab"
    for polarity in ["positive", "negative"]:
        name = f"C18_20220111_TPL_{polarity[:3].upper()}"
        new_atlas = generate_template_atlas(c18_data, ["Gold", "Platinum"], polarity, name)
        metob.store(new_atlas)
