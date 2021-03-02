import unittest

from oasis.oasis_io import PrintFileLoader

class IoTest(unittest.TestCase):

    def setUp(self):
        self.excel_path = './test/test_artifact/test_excel.xls'
        self.excel_p = PrintFileLoader(self.excel_path)

    def test_load_header(self):
        header = self.excel_p.load_header()
        self.assertEqual(header.LayerThickness, 0.06)
    
    def test_load_velocity_profile(self):
        profiles = self.excel_p.load_velocity_profile()
        self.assertEqual(profiles[2].JumpDelay, 100)
    
    def test_load_seg_style(self):
        styles = self.excel_p.load_segment_style()
        self.assertEqual(styles['Hatch'].Traveler.Power, 340)
    
    def test_load_parts(self):
        parts_setting = self.excel_p.load_parts()
        self.assertTrue(len(parts_setting[0]) > 0)
    
