import numpy as np
from .. import util


def test_sphere_dist_fast():
    assert(util.sphere_dist_fast(0, 0, 10, 0) == 10)
    assert(util.sphere_dist_fast(0, 0, 90, 0) == 90)
    assert(util.sphere_dist_fast(0, 0, 180, 0) == 180)
    # Accross meridian
    assert(util.sphere_dist_fast(355, 0, 5, 0) == 10)
    assert(util.sphere_dist_fast(315, 0, 45, 0) == 90)
    assert(util.sphere_dist_fast(270, 0, 90, 0) == 180)


def test_sphere_dist_fast_arrays():
    lon = np.array([10, 20, 30])
    assert((util.sphere_dist_fast(lon, 0, 0, 0) == lon).all())
    assert((util.sphere_dist_fast(0, 0, lon, 0) == lon).all())

    lon = np.array([355, 315])
    expected = np.array([10, 50])
    assert((util.sphere_dist_fast(lon, 0, 5, 0) == expected).all())
    assert((util.sphere_dist_fast(5, 0, lon, 0) == expected).all())

    lon1 = np.array([355, 315])
    lon2 = np.array([5, 45])
    expected = np.array([10, 90])
    assert((util.sphere_dist_fast(lon1, 0, lon2, 0) == expected).all())
    assert((util.sphere_dist_fast(lon2, 0, lon1, 0) == expected).all())
    lat = np.array([0, 0])
    assert((util.sphere_dist_fast(lon1, lat, lon2, lat) == expected).all())
    assert((util.sphere_dist_fast(lon2, lat, lon1, lat) == expected).all())
