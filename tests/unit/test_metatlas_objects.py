"""Test of metatlas objects"""
# pylint: disable=missing-function-docstring,protected-access,unused-argument,too-many-arguments

import getpass

import dill

from metatlas.datastructures import metatlas_objects as mo
from metatlas.datastructures import object_helpers as metoh

ADENOSINE_INCHI = 'InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)/t4-,6-,7-,10-/m1/s1'


def test_clone01():
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


def test_simple(sqlite):
    test = mo.Group()
    uid = test.unique_id
    mo.store(test)
    assert test.unique_id == uid
    assert test.prev_uid != ""
    test.name = "hello"
    mo.store(test)
    assert test.unique_id == uid
    assert test.prev_uid != ""


def test_nested(sqlite):
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    assert len(test.items) == 2
    test.items[1].name = "hello"
    orig_sub_version = test.items[1].unique_id
    assert len(test.items) == 2
    mo.store(test)
    assert test.items[1].unique_id == orig_sub_version


def test_recover(sqlite):
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    test.name = "howdy"
    top_version = test.unique_id
    sub_version = test.items[1].unique_id

    mo.store(test)
    mo.store(test)  # should have no effect
    assert len(test.items) == 2
    assert test.unique_id == top_version

    # make sure we can recover the previous version
    test.items = []
    assert test.unique_id == top_version
    test = mo.retrieve("group", unique_id=top_version)[0]
    assert test.unique_id == top_version
    assert len(test.items) == 2, len(test.items)
    assert test.unique_id == top_version
    assert test.items[1].unique_id == sub_version


def test_unique_links(sqlite):
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    sub_version = test.items[1].unique_id
    test.items = [test.items[1]]
    mo.store(test)

    test.items = []
    test = mo.retrieve("group", unique_id=test.unique_id)[0]
    assert len(test.items) == 1, len(test.items)
    assert test.items[0].unique_id == sub_version


def test_circular_reference(sqlite):
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    orig_id = test.unique_id
    test.items[0].items.append(test)
    mo.store(test)
    test.items = []
    test = mo.retrieve("group", unique_id=test.unique_id)[0]
    sub0 = test.items[0]
    assert len(sub0.items) == 2, sub0.items
    assert sub0.items[1].unique_id == orig_id
    assert test.unique_id == orig_id


def test_simple_query(sqlite):
    test1 = mo.LcmsRun(name="First")
    first_version = test1.unique_id
    test1.description = "Hey there"
    mo.store(test1)
    assert test1.unique_id == first_version
    items = mo.retrieve("lcmsrun", name="First")
    assert items[-1].unique_id == test1.unique_id
    assert all((i.unique_id != first_version for i in items[:-1]))


def test_glob_query(sqlite):
    test1 = mo.LcmsRun(name="First")
    test2 = mo.LcmsRun(name="Second")
    test3 = mo.LcmsRun(name="Third")
    mo.store([test1, test2, test3])
    items = mo.retrieve("lcmsrun", name="Fir%")
    assert items[-1].unique_id == test1.unique_id
    items = mo.retrieve("lcmsrun", name="%econd")
    assert items[-1].unique_id == test2.unique_id
    items = mo.retrieve("LcmsRuns", name="T%ir%")
    assert items[-1].unique_id == test3.unique_id


def test_escape_glob(sqlite):
    test1 = mo.LcmsRun(description="Flow %")
    mo.store(test1)
    items = mo.retrieve("lcmsrun", description="Flow %%")
    assert items[-1].unique_id == test1.unique_id


def test_id_grade_trait(sqlite):
    id_grade = mo.IdentificationGrade(name="E")
    mo.store(id_grade)
    cid = mo.CompoundIdentification(identification_grade="e")
    assert cid.identification_grade.unique_id == id_grade.unique_id


def test_list_item_changed():
    grp = mo.Group()
    grp.items = []
    grp._changed = False
    grp.items.append(mo.Group())
    assert grp._changed


def test_preserve_provenance(sqlite):
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    test2 = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    mo.store([test, test2])
    assert len(test.items) == 2
    test.items = []
    test2.items = []
    mo.store([test, test2])
    assert len(test.items) == 0
    previous = mo.retrieve("group", unique_id=test.prev_uid)[0]
    assert len(previous.items) == 2, repr(previous)


def test_store_stubs(sqlite):
    test = mo.Group(items=[mo.Group(items=[mo.LcmsRun()]), mo.LcmsRun()])
    mo.store(test)
    test = mo.retrieve("group", unique_id=test.unique_id)[0]
    assert isinstance(test.items[0], mo.Group)
    mo.store(test)


def test_user_preserve(sqlite):
    run = mo.LcmsRun(username="foo")
    test = mo.Reference(name="hello", username="foo", lcms_run=run)
    orig_id = test.unique_id
    mo.store(test, _override=True)
    assert test.unique_id == orig_id
    mo.store(test)
    assert test.unique_id != orig_id
    items = mo.retrieve("reference", username="*", name="hello")
    username = getpass.getuser()
    assert items[-2].username == "foo"
    assert items[-1].username == username
    assert items[-2].lcms_run.username == "foo"
    assert items[-1].lcms_run.username == "foo"


def test_store_all(sqlite):
    items = []
    for klass in metoh.Workspace.get_instance().subclass_lut.values():
        items.append(klass())
    mo.store(items)
    for klass in metoh.Workspace.get_instance().subclass_lut.values():
        name = klass.__name__
        assert len(mo.retrieve(name)) > 0


def test_stub_instance(sqlite):
    run = mo.LcmsRun(username="foo")
    test = mo.Reference(name="hello", lcms_run=run)
    mo.store(test, _override=True)
    item = mo.retrieve("reference", name="hello")[0]
    assert isinstance(item.lcms_run, mo.LcmsRun)


def test_floating_point(sqlite):
    compound = mo.Compound(name="foo", mono_isotopic_molecular_weight=1.0)
    mo.store(compound)
    compound.mono_isotopic_molecular_weight = 1.000007
    mo.store(compound)
    test = mo.retrieve("compound", name="foo")[-1]
    assert test.mono_isotopic_molecular_weight == 1.000007, test.mono_isotopic_molecular_weight


def test_remove(sqlite):
    group = mo.Group(name="foo", items=[mo.Group(name="baz", description="hello")])
    sub_id = group.items[0].unique_id
    mo.store(group)
    first = mo.retrieve("groups", unique_id=sub_id)[0]
    assert first.unique_id == sub_id
    mo.remove("groups", name="foo", _override=True)
    test = mo.retrieve("groups", name="foo")
    assert not test
    test_sub = mo.retrieve("groups_items", target_id=sub_id)
    assert not test_sub


def test_remove_objects(sqlite):
    group = mo.Group(name="foo", items=[mo.Group(name="baz", description="hello")])
    sub_id = group.items[0].unique_id
    mo.store(group)
    first = mo.retrieve("groups", unique_id=sub_id)[0]
    assert first.unique_id == sub_id
    mo.remove_objects(group, _override=True)
    test = mo.retrieve("groups", name="foo")
    assert not test
    test_sub = mo.retrieve("groups_items", target_id=sub_id)
    assert not test_sub


def test_dill():
    test = mo.Group(items=[mo.Group(description="hello")])
    blob = dill.dumps(test)
    new = dill.loads(blob)
    assert new.items[0].description == "hello"  # pylint: disable=no-member


def test_retrieve01(sqlite):
    compound = mo.Compound(name="foo", inchi=ADENOSINE_INCHI, inchi_key='foobar')
    mo.store(compound)
    assert mo.retrieve('Compounds', inchi_key=[], username="*") == []
    assert mo.retrieve('Compounds', inchi=[ADENOSINE_INCHI], username="*")[0].inchi == ADENOSINE_INCHI
