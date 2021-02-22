# code for testing our package
import math
import unittest
import psvpy.psvpy as psvpy

class TestPSVpy(unittest.TestCase):
    
    def test_steamPSV(self):
    # test case steam
    # From API 520, Part I, 2000 (section 3.7.2)
    # W = 69615 kg/h
    # P = 11032*1.1 + 101.3 = 12236 kPa
    # saturated
    # API gets 1100 mm2, and this is with a Napier factor of 51.5, we are using 51.45
        testFlow = 69615.0 # kg.h
        testP = 11032*1.1 + 101.3
        targetSize = 1100.0
        self.assertTrue(abs(psvpy.PSVsteamSize(testFlow, testP, "Sat") - targetSize) < 2.0)

    def test_vapourPSV(self):
    # test case vapour
    # From API 520, Part I, 2013 (section 5.6.3.2.3)
    # W = 24270.0 kg/h
    # P = 670.0 kPa
    # MW = 51
    # Z = 0.9
    # k = Cp/Cv = 1.11, or C = 328.
    # T = 348 - 273
    # API gets 3698 mm2
    # W, P, Tcelcius, MW, k, Z
        testFlow = 24270.0
        testP = 670.0
        testTC = 348.0 - 273.0
        testMW = 51.0
        testk = 1.11
        testZ = 0.9
        targetSize = 3698.0
        self.assertTrue(abs(psvpy.PSVvaporSize(testFlow, testP, testTC, testMW, testk, testZ) - targetSize) < 2.0)

    def test_liquid1(self):
    # test case liquid1
    # from API 520, Part I, 2013 (section 5.8.1.3)
    # flow = 6814 litres/min
    # density = 900 kg/m3
    # P = 1896 kpag (incudes 10%)
    # Pback = 345 kPag
    # viscosity 2000 SSU (API uses archaic units, the table gives 431.7 cSt)
    # API gets 3180 mm2 with a backpressure correction of 0.97, or 3085 without the backpressure correction
    # PSVliquidSize(W, P, Pback, d, mu)
    # PSVliquidSize(900.0*6814*60.0/1000.0, 1896.0, 345.0, 900.0, 431.7*0.900)
        testRho = 900.0
        testBackP = 345.0
        testVisc = 431.7*0.900
        testSetP = 1896.0
        targetFlow = testRho*6814*60.0/1000.0
        self.assertTrue(abs(psvpy.PSVliquidSize(targetFlow, testSetP, testBackP, testRho, testVisc) - 3085.0) < 10.0)

    def test_liquid2(self):
    # test case liquid2
    # from API 520, Part I, 2013 (section 5.8.1.3)
    # flow = 6814 litres/min
    # density = 900 kg/m3
    # P = 1896 kpag (incudes 10%)
    # Pback = 345 kPag
    # viscosity 2000 SSU (API uses archaic units, the table gives 431.7 cSt)
    # we will flip around the sizing case to get the rating case
    # API gets 3180 mm2 with a backpressure correction of 0.97, or 3085 without the backpressure correciton
    # PSVliquidSize(W, P, Pback, d, mu)
        testRho = 900.0
        testBackP = 345.0
        testVisc = 431.7*0.900
        testSetP = 1896.0
        targetFlow = testRho*6814*60.0/1000.0
        targetArea = psvpy.PSVliquidSize(targetFlow, testSetP, testBackP, testRho, testVisc)
        self.assertTrue(abs(psvpy.PSVliquidRate(targetArea, testSetP, testBackP, testRho, testVisc) - targetFlow) < 1.0)

if __name__ == '__main__':
    unittest.main()

