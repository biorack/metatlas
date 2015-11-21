
from metatlas import metatlas_objects as mo
from metatlas.mzml_loader import get_test_data
import getpass
import dill


def test_simple():
    test = mo.Group()
    uid = test.unique_id
    mo.store(test)
    assert test.unique_id == uid
    assert test.prev_uid != ''
    test.name = 'hello'
    mo.store(test)
    assert test.unique_id == uid
    assert test.prev_uid != ''


def test_nested():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    assert len(test.items) == 2
    test.items[1].name = 'hello'
    orig_sub_version = test.items[1].unique_id
    assert len(test.items) == 2
    mo.store(test)
    assert test.items[1].unique_id == orig_sub_version


def test_recover():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    test.name = 'howdy'
    top_version = test.unique_id
    sub_version = test.items[1].unique_id

    mo.store(test)
    mo.store(test)  # should have no effect
    assert len(test.items) == 2
    assert test.unique_id == top_version

    # make sure we can recover the previous version
    test.items = []
    assert test.unique_id == top_version
    test = mo.retrieve('group', unique_id=top_version)[0]
    assert test.unique_id == top_version
    assert len(test.items) == 2, len(test.items)
    assert test.unique_id == top_version
    assert test.items[1].unique_id == sub_version


def test_unique_links():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    sub_version = test.items[1].unique_id
    test.items = [test.items[1]]
    mo.store(test)

    test.items = []
    test = mo.retrieve('group', unique_id=test.unique_id)[0]
    assert len(test.items) == 1, len(test.items)
    assert test.items[0].unique_id == sub_version


def test_circular_reference():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    orig_id = test.unique_id
    test.items[0].items.append(test)
    mo.store(test)
    test.items = []
    test = mo.retrieve('group', unique_id=test.unique_id)[0]
    sub0 = test.items[0]
    assert len(sub0.items) == 2, sub0.items
    assert sub0.items[1].unique_id == orig_id
    assert test.unique_id == orig_id


def test_simple_query():
    test1 = mo.LcmsRun(name='First')
    first_version = test1.unique_id
    test1.description = "Hey there"
    mo.store(test1)
    assert test1.unique_id == first_version
    items = mo.retrieve('lcmsrun', name='First')
    assert items[-1].unique_id == test1.unique_id
    assert all([i.unique_id != first_version for i in items[:-1]])


def test_glob_query():
    test1 = mo.LcmsRun(name='First')
    test2 = mo.LcmsRun(name='Second')
    test3 = mo.LcmsRun(name='Third')
    mo.store([test1, test2, test3])
    items = mo.retrieve('lcmsrun', name='Fir%')
    assert items[-1].unique_id == test1.unique_id
    items = mo.retrieve('lcmsrun', name='%econd')
    assert items[-1].unique_id == test2.unique_id
    items = mo.retrieve('LcmsRuns', name='T%ir%')
    assert items[-1].unique_id == test3.unique_id


def test_escape_glob():
    test1 = mo.LcmsRun(description='Flow %')
    mo.store(test1)
    items = mo.retrieve('lcmsrun', description='Flow %%')
    assert items[-1].unique_id == test1.unique_id


def test_load_lcms_files():
    paths = get_test_data().values()
    runs = mo.load_lcms_files(paths)
    for run in runs:
        assert run.mzml_file
        assert run.hdf5_file
        assert run.creation_time
        assert run.description
        assert run.name
        assert run.last_modified
        assert run.username
        assert run.unique_id
        assert mo.retrieve('lcmsrun', unique_id=run.unique_id)


def test_id_grade_trait():
    e = mo.IdentificationGrade(name='E')
    mo.store(e)
    cid = mo.CompoundIdentification(identification_grade='e')
    assert cid.identification_grade.unique_id == e.unique_id


def test_list_item_changed():
    g = mo.Group()
    g.items = []
    g._changed = False
    g.items.append(mo.Group())
    assert g._changed


def test_preserve_provenance():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    test2 = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    mo.store([test, test2])
    assert len(test.items) == 2
    test.items = []
    test2.items = []
    mo.store([test, test2])
    assert len(test.items) == 0
    print(test.unique_id)
    previous = mo.retrieve('group', unique_id=test.prev_uid)[0]
    assert len(previous.items) == 2, repr(previous)


def test_clone():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    test2 = test.clone()
    assert test2.unique_id != test.unique_id
    assert test2.prev_uid == test.unique_id
    assert test2.items[0].unique_id == test.items[0].unique_id

    test3 = test.clone(True)
    assert test3.unique_id != test.unique_id
    assert test3.prev_uid == test.unique_id
    assert test3.items[0].unique_id != test.items[0].unique_id
    assert test3.items[0].prev_uid == test.items[0].unique_id


def test_store_stubs():
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    mo.store(test)
    test = mo.retrieve('group', unique_id=test.unique_id)[0]
    assert isinstance(test.items[0], mo.Group)
    mo.store(test)


def test_get_latest():
    test = mo.Compound(name='hello')
    mo.store(test)
    test.name = 'goodbye'
    mo.store(test)
    test = mo.retrieve('compound', creation_time=test.creation_time)
    assert len(test) == 1, len(test)
    assert test[0].name == 'goodbye'


def test_user_preserve():
    run = mo.LcmsRun(username='foo')
    test = mo.Reference(name='hello', username='foo', lcms_run=run)
    orig_id = test.unique_id
    mo.store(test, _override=True)
    assert test.unique_id == orig_id
    mo.store(test)
    assert test.unique_id != orig_id
    items = mo.retrieve('reference', username='*', name='hello')
    username = getpass.getuser()
    assert items[-2].username == 'foo'
    assert items[-1].username == username
    assert items[-2].lcms_run.username == 'foo'
    assert items[-1].lcms_run.username == 'foo'
    run.name = 'hello'
    mo.store(test)
    items = mo.retrieve('reference', username='*',
                        creation_time=test.creation_time)
    return
    assert items[0].lcms_run.username == 'foo'
    assert items[1].lcms_run.username == username


def test_store_all():
    items = []
    for klass in mo.SUBCLASS_LUT.values():
        items.append(klass())
    mo.store(items)
    for klass in mo.SUBCLASS_LUT.values():
        name = klass.__name__
        assert len(mo.retrieve(name))


def test_stub_instance():
    run = mo.LcmsRun(username='foo')
    test = mo.Reference(name='hello', lcms_run=run)
    mo.store(test)
    item = mo.retrieve('reference', name='hello')[0]
    assert isinstance(item.lcms_run, mo.LcmsRun)


def test_floating_point():
    compound = mo.Compound(name='foo', MonoIsotopic_molecular_weight=1.0)
    mo.store(compound)
    compound.MonoIsotopic_molecular_weight = 1.000007
    mo.store(compound)
    test = mo.retrieve('compound', name='foo')[-1]
    assert test.MonoIsotopic_molecular_weight == 1.000007, test.MonoIsotopic_molecular_weight


def test_remove():
    compound = mo.Compound(name='foo', MonoIsotopic_molecular_weight=1.0,
                           reference_xrefs=[mo.ReferenceDatabase(name='baz')])
    sub_id = compound.reference_xrefs[0].unique_id
    mo.store(compound)
    db = mo.retrieve('referencedatabase', unique_id=sub_id)[0]
    assert db.unique_id == sub_id
    mo.remove('compound', name='foo', _override=True)
    test = mo.retrieve('compound', name='foo')
    assert not test
    test_sub = mo.retrieve('compounds_reference_xrefs', target_id=sub_id)
    assert not test_sub


def test_remove_objects():
    compound = mo.Compound(name='foo', MonoIsotopic_molecular_weight=1.0,
                           reference_xrefs=[mo.ReferenceDatabase(name='baz')])
    sub_id = compound.reference_xrefs[0].unique_id
    mo.store(compound)
    db = mo.retrieve('referencedatabase', unique_id=sub_id)[0]
    assert db.unique_id == sub_id
    mo.remove_objects(compound, _override=True)
    test = mo.retrieve('compound', name='foo')
    assert not test
    test_sub = mo.retrieve('compounds_reference_xrefs', target_id=sub_id)
    assert not test_sub


def test_dill():
    test = mo.Group(items=[mo.Group(description='hello')])
    blob = dill.dumps(test)
    new = dill.loads(blob)
    assert new.items[0].description == 'hello'


def test_retrieve_head():
    test = mo.LcmsRun(name='foo')
    mo.store(test)
    old = len(mo.retrieve('lcmsrun', name='foo'))
    test.description = 'bar'
    mo.store(test)
    new = len(mo.retrieve('lcmsrun', name='foo'))
    assert new == old
