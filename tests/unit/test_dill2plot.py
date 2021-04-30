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
