from obarasaika.obara_saika import get_coulomb, get_overlap, get_kinetic, get_nuclear, get_moment


def test_get_coulomb():
    """A combined test for get_coulomb()."""

    za = 1.1
    zb = 1.2
    zc = 1.3
    zd = 1.4

    ra = [1.0, 0.0, 1.0]
    rb = [0.0, 1.0, 2.0]
    rc = [0.0, 0.0, 3.0]
    rd = [0.0, 0.0, 4.0]

    thresh = 1.0e-16

    ref = 0.08608517834596989
    integral = get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert abs(integral - ref) < thresh

    ref = -0.046888060557776973
    integral = get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    assert abs(integral - ref) < thresh

    ref = 1.71817807954e-05
    integral = get_coulomb(za, zb, zc, zd, ra, rb, rc, rd, [2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0])
    assert abs(integral - ref) < thresh


def test_get_overlap():
    """A combined test for get_overlap()."""

    za = 1.8
    zb = 2.8
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    thresh = 1.0e-16

    integral = get_overlap(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])
    assert abs(integral - 0.20373275913014607) < thresh

    integral = get_overlap(za, zb, ra, rb, [1, 0, 0, 0, 0, 0])
    assert abs(integral - 0.062005622343957505) < thresh

    integral = get_overlap(za, zb, ra, rb, [1, 1, 0, 1, 1, 0])
    assert abs(integral - -0.00043801221837779696) < thresh

    integral = get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0])
    assert abs(integral - -0.0002385994651113168) < thresh


def test_get_kinetic():
    """A combined test for get_kinetic()."""

    za = 1.8
    zb = 2.0
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]

    thresh = 1.0e-16

    integral = get_kinetic(za, zb, ra, rb, [0, 0, 0, 0, 0, 0])
    assert abs(integral - 0.3652714583525358) < thresh

    integral = get_kinetic(za, zb, ra, rb, [1, 0, 0, 0, 0, 0])
    assert abs(integral - 0.2514265587836556) < thresh

    integral = get_kinetic(za, zb, ra, rb, [2, 2, 2, 2, 2, 2])
    assert abs(integral - -7.40057384314e-05) < thresh


def test_get_nuclear():
    """A combined test for get_nuclear()."""

    za = 1.8
    zb = 2.0
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]
    rc = [0.5, 0.8, 0.2]

    # TODO why is this thresh lower now?
    thresh = 1.0e-15

    integral = get_nuclear(za, zb, ra, rb, rc, [0, 0, 0, 0, 0, 0])
    assert abs(integral - -0.49742209545104593) < thresh

    integral = get_nuclear(za, zb, ra, rb, rc, [1, 0, 0, 0, 0, 0])
    assert abs(integral - -0.15987439458254471) < thresh

    integral = get_nuclear(za, zb, ra, rb, rc, [2, 2, 2, 0, 0, 0])
    assert abs(integral - -0.003801373531942607) < thresh

    integral = get_nuclear(za, zb, ra, rb, rc, [1, 1, 1, 1, 1, 1])
    assert abs(integral - 8.8415484347060993e-5) < thresh


def test_get_moment():
    """A combined test for get_moment()."""

    za = 1.8
    zb = 2.0
    ra = [0.0, 0.0, 0.0]
    rb = [0.5, 0.8, -0.2]
    rc = [0.0, 0.0, 0.0]

    thresh = 1.0e-16

    integral = get_moment(za, zb, ra, rb, rc, [0, 0, 2, 0, 0, 0], [0, 0, 1])
    assert abs(integral - -0.01330515491323708) < thresh


if __name__ == '__main__':
    test_get_coulomb()
    test_get_overlap()
    test_get_kinetic()
    test_get_nuclear()
    test_get_moment()
