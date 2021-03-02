import unittest

from oasis.slicing import slice_object

class SliceTest(unittest.TestCase):

    def setUp(self):
        self.input_file = './test/test_artifact/3DBenchy.stl'

    def test_slicing(self):
        slice_object(self.input_file, 0.1, 0.1, 0.2, 2)