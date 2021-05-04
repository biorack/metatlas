from metatlas.plots import dill2plots


def test_rt_range_overlaps_none():
    i = type('', (), {})()
    i.rt_min = 10
    i.rt_max = 20
    j = type('', (), {})()
    j.rt_min = 30
    j.rt_max = 40
    assert not dill2plots.rt_range_overlaps(i, j)
    assert not dill2plots.rt_range_overlaps(j, i)


def test_rt_range_overlaps_partial():
    i = type('', (), {})()
    i.rt_min = 10
    i.rt_max = 20
    j = type('', (), {})()
    j.rt_min = 15
    j.rt_max = 25
    assert dill2plots.rt_range_overlaps(i, j)
    assert dill2plots.rt_range_overlaps(j, i)


def test_rt_range_overlaps_inside():
    i = type('', (), {})()
    i.rt_min = 10
    i.rt_max = 20
    j = type('', (), {})()
    j.rt_min = 12
    j.rt_max = 18
    assert dill2plots.rt_range_overlaps(i, j)
    assert dill2plots.rt_range_overlaps(j, i)


def test_within_tolerance_yes():
    assert dill2plots.within_tolerance(99, 100, 0.02)


def test_within_tolerance_no():
    assert not dill2plots.within_tolerance(99, 100, 0.002)


def test_filter_metatlas_objects_by_list_remove_all():
    i = type('', (), {})()
    i.myattr = [1, 2]
    j = type('', (), {})()
    j.myattr = [2, 3]
    assert [] == dill2plots.filter_metatlas_objects_by_list([i, j], 'myattr', [4, 5])


def test_filter_metatlas_objects_by_list_remove_none():
    i = type('', (), {})()
    i.myattr = [1, 2]
    j = type('', (), {})()
    j.myattr = [2, 3]
    assert [i, j] == dill2plots.filter_metatlas_objects_by_list([i, j], 'myattr', [2])


def test_remove_metatlas_objects_by_list_remove_none():
    i = type('', (), {})()
    i.myattr = [1, 2]
    j = type('', (), {})()
    j.myattr = [2, 3]
    assert [i, j] == dill2plots.remove_metatlas_objects_by_list([i, j], 'myattr', [4, 5])


def test_remove_metatlas_objects_by_list_remove_all():
    i = type('', (), {})()
    i.myattr = [1, 2]
    j = type('', (), {})()
    j.myattr = [2, 3]
    assert [] == dill2plots.remove_metatlas_objects_by_list([i, j], 'myattr', [0, 2, 5])


def test_export_atlas_to_spreadsheet(atlas):
    expected = """{"chebi_id":{"0":"CHEBI:17256"},"chebi_url":{"0":"http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:17256"},"creation_time":{"0":1466212395.0},"description":{"0":"A purine 2'-deoxyribonucleoside having adenine as the nucleobase."},"formula":{"0":"C10H13N5O3"},"head_id":{"0":"60cd6743e56545c6a6cb066ec3553450"},"hmdb_id":{"0":"HMDB00101"},"hmdb_url":{"0":"http://www.hmdb.ca/metabolites/HMDB00101"},"img_abc_id":{"0":""},"inchi":{"0":"InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"},"inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"iupac_name":{"0":""},"kegg_id":{"0":"C00559"},"kegg_url":{"0":"http://www.genome.jp/dbget-bin/www_bget?C00559"},"last_modified":{"0":1612996604.0},"lipidmaps_id":{"0":""},"lipidmaps_url":{"0":""},"metacyc_id":{"0":"DEOXYADENOSINE"},"mono_isotopic_molecular_weight":{"0":251.101839276},"name":{"0":"2'-deoxyadenosine"},"neutralized_2d_inchi":{"0":"InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)"},"neutralized_2d_inchi_key":{"0":"OLXZPDWKRNYJJZ-UHFFFAOYSA-N"},"neutralized_inchi":{"0":"InChI=1S/C10H13N5O3/c11-9-8-10(13-3-12-9)15(4-14-8)7-1-5(17)6(2-16)18-7/h3-7,16-17H,1-2H2,(H2,11,12,13)/t5-,6+,7+/m0/s1"},"neutralized_inchi_key":{"0":"OLXZPDWKRNYJJZ-RRKCRQDMSA-N"},"num_free_radicals":{"0":0.0},"number_components":{"0":1.0},"permanent_charge":{"0":0.0},"prev_uid":{"0":"origin"},"pubchem_compound_id":{"0":"13730"},"pubchem_url":{"0":"http://pubchem.ncbi.nlm.nih.gov/compound/13730"},"source":{"0":"gnps///chebi///metacyc///hmdb"},"synonyms":{"0":"2'-deoxyadenosine"},"unique_id":{"0":"60cd6743e56545c6a6cb066ec3553450"},"username":{"0":"wjholtz"},"wikipedia_url":{"0":""},"label":{"0":"2'-deoxyadenosine"},"id_notes":{"0":"No description"},"ms1_notes":{"0":"keep"},"ms2_notes":{"0":"bad match to ref"},"identification_notes":{"0":"my id note"},"rt_min":{"0":1.6964640054},"rt_max":{"0":2.6964640054},"rt_peak":{"0":2.1964640054},"mz":{"0":252.1091393},"mz_tolerance":{"0":20.0},"adduct":{"0":"[M+H]+"},"polarity":{"0":"positive"}}"""  # noqa: E501
    assert expected == dill2plots.export_atlas_to_spreadsheet(atlas).to_json()
