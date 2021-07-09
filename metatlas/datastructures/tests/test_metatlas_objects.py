from metatlas.datastructures import metatlas_objects as mo


def test_get_latest():
    test = mo.Compound(name='hello')
    mo.store(test)
    test.name = 'goodbye'
    mo.store(test)
    test = mo.retrieve('compound', creation_time=test.creation_time)
    assert len(test) == 1, len(test)
    assert test[0].name == 'goodbye'


def test_retrieve_head():
    test = mo.LcmsRun(name='foo')
    mo.store(test)
    old = len(mo.retrieve('lcmsrun', name='foo'))
    test.description = 'bar'
    mo.store(test)
    new = len(mo.retrieve('lcmsrun', name='foo'))
    assert new == old
