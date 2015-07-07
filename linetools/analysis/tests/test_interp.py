from ..interp import interp_Akima

def test_interp_Akima():
    x = np.sort(np.random.random(10) * 10)
    y = np.random.normal(0.0, 0.1, size=len(x))
    assert np.allclose(y, interp_Akima(x, x, y))
