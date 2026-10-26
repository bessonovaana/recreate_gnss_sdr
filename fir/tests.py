import unittest
from fir.fir import apply_fir

def fir_filter(x, taps, decim=1):
    y = []
    x2 = (len(taps) - 1) * [0, ] + x
    for i in range(0, len(x), decim):
        yi = 0
        for j in range(len(taps)):
            yi += taps[len(taps) - 1 - j] * x2[i + j]
        y.append(yi)
    return y


class TestMyFirAgainstReference(unittest.TestCase):

    def test_fff_decim1(self):
        decim = 1
        taps = 20 * [0.5, 0.5]
        src_data = 40 * [1, 2, 3, 4]
        expected = fir_filter(src_data, taps, decim)
        result = apply_fir(src_data, taps)

        self.assertEqual(len(expected), len(result))
        for e, r in zip(expected, result):
            self.assertAlmostEqual(e, r, places=5)


    def test_fsf_decim1(self):
        decim = 1
        taps = 20 * [0.5, 0.5]
        src_data = 40 * [1, 2, 3, 4]
        expected = [round(e) for e in fir_filter(src_data, taps, decim)]
        result = [round(e) for e in apply_fir(src_data, taps)]

        self.assertEqual(expected, result)


if __name__ == '__main__':
    unittest.main()