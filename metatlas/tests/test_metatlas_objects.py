
from metatlas import metatlas_objects as mo
from metatlas.mzml_loader import get_test_data


def test_simple():
    test = mo.MetatlasObject()
    uid = test.unique_id
    test.store()
    assert test.unique_id == uid
    test.name = 'hello'
    test.store()
    assert test.unique_id != uid


def test_nested():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    assert len(test.items) == 2
    test.items[1].name = 'hello'
    orig_sub_version = test.items[1].unique_id
    assert len(test.items) == 2
    test.store()
    assert test.items[1].unique_id != orig_sub_version


def test_recover():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    top_version = test.unique_id
    sub_version = test.items[1].unique_id
    test.store()

    test.store()  # should have no effect
    assert len(test.items) == 2
    assert test.unique_id == top_version

    # make sure we can recover the previous version
    test.items = []
    assert test.unique_id == top_version
    test.retrieve()
    assert test.unique_id == top_version
    assert len(test.items) == 2, len(test.items)
    assert test.unique_id == top_version
    assert test.items[1].unique_id == sub_version


def test_unique_links():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    sub_version = test.items[1].unique_id
    test.items = [test.items[1]]
    test.store()

    test.items = []
    test.retrieve()
    assert len(test.items) == 1, len(test.items)
    assert test.items[0].unique_id == sub_version


def test_circular_reference():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    orig_id = test.unique_id
    test.items[0].items.append(test)
    test.store()
    test.items = []
    test.retrieve()
    assert len(test.items[0].items) == 2, len(test.items)
    assert test.items[0].items[1].unique_id == orig_id
    assert test.unique_id == orig_id


def test_simple_query():
    test1 = mo.LcmsRun(name='First')
    first_version = test1.unique_id
    test1.description = "Hey there"
    test1.store()
    assert test1.unique_id != first_version
    items = mo.queryDatabase('lcmsrun', name='First')
    assert items[-1].unique_id == test1.unique_id
    assert all([i.unique_id != first_version for i in items])


def test_glob_query():
    test1 = mo.LcmsRun(name='First')
    test2 = mo.LcmsRun(name='Second')
    test3 = mo.LcmsRun(name='Third')
    test1.store()
    test2.store()
    test3.store()
    items = mo.queryDatabase('lcmsrun', name='Fir%')
    assert items[-1].unique_id == test1.unique_id
    items = mo.queryDatabase('lcmsrun', name='%econd')
    assert items[-1].unique_id == test2.unique_id
    items = mo.queryDatabase('LcmsRuns', name='T%ir%')
    assert items[-1].unique_id == test3.unique_id


def test_escape_glob():
    test1 = mo.LcmsRun(description='Flow %')
    test1.store()
    items = mo.queryDatabase('lcmsrun', description='Flow %%')
    assert items[-1].unique_id == test1.unique_id


def test_select_reference_by_type():
    mz_refs = [mo.MzReference(name=str(i)) for i in range(3)]
    rt_refs = [mo.RtReference(name=str(i)) for i in range(4)]
    frag_refs = [mo.FragmentationReference(name=str(i)) for i in range(2)]

    compound_id = mo.CompoundId(name='test',
                                references=mz_refs + rt_refs + frag_refs)
    compound_id.store()

    assert mo.queryDatabase('mzreference', name='1')[0].name == '1'

    mz_select = compound_id.select_by_type('mz')
    assert mz_select[0] is mz_refs[0]
    assert len(mz_select) == 3

    assert len(compound_id.select_by_type('rt')) == 4
    assert len(compound_id.select_by_type('frag')) == 2


def test_load_lcms_files():
    paths = get_test_data().values()
    runs = mo.load_lcms_files(paths)
    for run in runs:
        assert run.mzml_file
        assert run.hdf5_file
        assert run.created
        assert run.created_by
        assert run.description
        assert run.name
        assert run.last_modified
        assert run.modified_by == run.created_by
        assert run.unique_id
        assert not run.prev_unique_id
        assert mo.queryDatabase('lcmsrun', unique_id=run.unique_id)


def test_id_grade_trait():
    e = mo.IdentificationGrade(name='E')
    e.store()
    cid = mo.CompoundId(identification_grade='e')
    assert cid.identification_grade.unique_id == e.unique_id
