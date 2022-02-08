"""Populate compound fields"""
# pylint: disable=missing-function-docstring

import logging
import time
from typing import List, Optional
from urllib.parse import quote

import pandas as pd
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

from metatlas.datastructures import metatlas_objects as metob
from metatlas.datastructures.object_helpers import MetUnicode
from metatlas.plots import dill2plots as dp
from metatlas.tools import cheminfo

logger = logging.getLogger(__name__)


def generate_template_atlas(
    raw_file_name: str, confidence_levels: List[str], polarity: str, name: str, mz_tolerance: float = 10
) -> metob.Atlas:
    data = pd.read_csv(raw_file_name, sep="\t")
    acceptable = data[data["confidence_category"].isin(confidence_levels)]
    by_polarity = acceptable[acceptable["polarity"] == polarity]
    by_polarity = by_polarity.assign(label=None)
    atlas = dp.make_atlas_from_spreadsheet(
        by_polarity, name, filetype="dataframe", polarity=polarity, store=False, mz_tolerance=mz_tolerance
    )
    inchi_keys = [cid.compound[0].inchi_key for cid in atlas.compound_identifications]
    pubchem_results = query_pubchem(inchi_keys)
    for cid in atlas.compound_identifications:
        fill_fields(cid.compound[0], pubchem_results)
        cid.name = cid.compound[0].name
    return atlas


def count_non_empty(compound: metob.Compound) -> int:
    return sum([v != "" for v in compound.trait_values().values()])


def flatten_inchi(mol: Chem.rdchem.Mol) -> str:
    smiles = Chem.MolToSmiles(mol).replace("@", "")
    flattened_rdkit_mol = Chem.MolFromSmiles(smiles)
    try:
        return Chem.MolToInchi(flattened_rdkit_mol)
    except Exception:  # This fails when can't kekulize mol # pylint: disable=broad-except
        logger.warning("failed to flatten a molecule")
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


def set_all_ids(comp: metob.Compound):
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
        set_id(comp, id_type["metatlas_name"], id_type["cts_name"], id_type["base_url"], comp.inchi_key)


def fill_neutralized_fields(comp: metob.Compound, mol: Chem.rdchem.Mol):
    try:
        norm_mol = cheminfo.normalize_molecule(mol)
    except Exception:  # pylint: disable=broad-except
        logger.warning("failed to normalized %s", comp.name)
        return
    assert norm_mol is not None
    if not comp.neutralized_inchi:
        comp.neutralized_inchi = Chem.inchi.MolToInchi(norm_mol)
    if not comp.neutralized_inchi_key:
        comp.neutralized_inchi_key = Chem.inchi.InchiToInchiKey(comp.neutralized_inchi)
    if not comp.neutralized_2d_inchi:
        comp.neutralized_2d_inchi = flatten_inchi(norm_mol)  # type: ignore
    if not comp.neutralized_2d_inchi_key:
        comp.neutralized_2d_inchi_key = Chem.inchi.InchiToInchiKey(comp.neutralized_2d_inchi)


def fill_calculated_fields(comp: metob.Compound, mol: Chem.rdchem.Mol):
    assert mol is not None
    comp.inchi_key = comp.inchi_key or Chem.inchi.InchiToInchiKey(comp.inchi)
    comp.formula = comp.formula or Chem.rdMolDescriptors.CalcMolFormula(mol)
    comp.mono_isotopic_molecular_weight = comp.mono_isotopic_molecular_weight or ExactMolWt(mol)
    comp.permanent_charge = comp.permanent_charge or Chem.GetFormalCharge(mol)
    comp.number_components = comp.number_components or 1  # type: ignore
    comp.num_free_radicals = comp.num_free_radicals or Chem.Descriptors.NumRadicalElectrons(mol)
    fill_neutralized_fields(comp, mol)


def first_all_ascii(list_of_strings: List[str]) -> str:
    for to_check in list_of_strings:
        if to_check.isascii():
            return to_check
    raise ValueError("No strings found with only ASCII characters")


def filter_out_strings_with_non_ascii(list_of_strings: List[str]) -> List[str]:
    return [s for s in list_of_strings if s.isascii()]


def fill_fields(comp: metob.Compound, pubchem_results: List[pcp.Compound]):
    """
    Populate blank fields that can be infered from other fields.
    Does not overwrite any existing values that are not None, '', or 'Untitled'"""
    mol = Chem.inchi.MolFromInchi(comp.inchi)
    if mol is None:
        return
    fill_calculated_fields(comp, mol)
    set_all_ids(comp)
    pubchem = get_pubchem_compound(comp.inchi_key, pubchem_results)
    if pubchem is not None:
        if not comp.pubchem_compound_id:
            comp.pubchem_compound_id = pubchem.cid
        if not comp.pubchem_url:
            comp.pubchem_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{comp.pubchem_compound_id}"
        if not comp.synonyms:
            comp.synonyms = "///".join(filter_out_strings_with_non_ascii(pubchem.synonyms))
        if not comp.iupac_name:
            comp.iupac_name = pubchem.iupac_name
    if comp.name in ["", "Untitled"] or "///" in comp.name:
        names = [first_all_ascii(comp.synonyms.split("///"))] + [comp.iupac_name]
        comp.name = names[0]


def create_c18_template_atlases():
    c18_data = "/global/u2/w/wjholtz/c18_atlas_creation.tab"
    for polarity in ["negative", "positive"]:
        name = f"C18_20220118_TPL_{polarity[:3].upper()}"
        new_atlas = generate_template_atlas(c18_data, ["Gold", "Platinum"], polarity, name)
        metob.store(new_atlas)


# pylint: disable=too-many-arguments
def generate_stds_atlas(
    raw_file_name: str,
    inchi_keys: List[str],
    polarity: str,
    name: str,
    mz_tolerance: float = 10,
    more_rows: Optional[pd.DataFrame] = None,
) -> metob.Atlas:
    data = pd.read_csv(raw_file_name, sep="\t")
    if more_rows is not None:
        data = data.append(more_rows)
    acceptable = data[data["inchi_key"].isin(inchi_keys)]
    by_polarity = acceptable[acceptable["polarity"] == polarity]
    by_polarity = by_polarity.assign(label=None)
    return make_atlas_from_df(by_polarity, name, polarity, mz_tolerance)


def fill_atlas_compound_fields(atlas):
    inchi_keys = [cid.compound[0].inchi_key for cid in atlas.compound_identifications]
    pubchem_results = query_pubchem(inchi_keys)
    for cid in atlas.compound_identifications:
        fill_fields(cid.compound[0], pubchem_results)
        cid.name = cid.compound[0].name
    return atlas


def make_atlas_from_df(df, name, polarity, mz_tolerance):
    atlas = dp.make_atlas_from_spreadsheet(
        df, name, filetype="dataframe", polarity=polarity, store=False, mz_tolerance=mz_tolerance
    )
    return fill_atlas_compound_fields(atlas)


def create_c18_stds_atlases(mz_tolerance: float = 10) -> None:
    c18_path = "/global/u2/w/wjholtz/c18_atlas_creation.tab"
    data = pd.read_csv(c18_path, sep="\t")
    std_inchi_keys = {
        "Phenylalanine": "COLNVLDHVKWLRT-QMMMGPOBSA-N",
        "L-Tryptophan": "QIVBCDIJIAJPQS-SECBINFHSA-N",
        "Salicylic acid": "YGSDEFSMJLZEOE-UHFFFAOYSA-N",
        "2-Amino-3-bromo-5-methylbenzoic acid": "LCMZECCEEOQWLQ-UHFFFAOYSA-N",  # this one will not be found in c18_data
    }
    abmba = "2-Amino-3-bromo-5-methylbenzoic acid"
    for polarity in ["positive", "negative"]:
        more_rows = pd.DataFrame(
            {
                "inchi_key": [std_inchi_keys[abmba]],
                "label": [abmba],
                "adduct": ["[M+H]+" if polarity == "positive" else "[M-H]-"],
                "polarity": [polarity],
                "rt_min": [4.5],
                "rt_peak": [4.7],
                "rt_max": [4.9],
                "mz": [228.97384 + (1.00727647 * (1 if polarity == "positive" else -1))],
                "confidence_category": "Platinum",
            }
        )
        if more_rows is not None:
            data = data.append(more_rows)
        acceptable = data[data["inchi_key"].isin(std_inchi_keys.values())]
        by_polarity = acceptable[acceptable["polarity"] == polarity]
        by_polarity = by_polarity.assign(label=None)
        by_polarity["rank"] = by_polarity["confidence_category"] == "Platinum"
        single = by_polarity.loc[by_polarity.groupby(["inchi_key"])["rank"].idxmax()]
        name = f"C18_20220208c_QC_{polarity[:3].upper()}"
        atlas = make_atlas_from_df(single, name, polarity, mz_tolerance)
        metob.store(atlas)
