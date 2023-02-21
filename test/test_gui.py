import unittest

"""
TODO - implement meaningful tests here (how to do this?).
"""
class MyTestCase(unittest.TestCase):
    def test_something(self):
        self.assertEqual(True, True)

    @classmethod
    def setUpClass(cls):
        pass # todo
    @classmethod
    def tearDownClass(cls):
        pass

if __name__ == '__main__':
    unittest.main()
