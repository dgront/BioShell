import unittest
from pybioshell.src.cif import cif


class TestPyCifLoop(unittest.TestCase):

    def setUp(self):
     
        self.loop = PyCifLoop()

    def test_new_instance(self):
 
        self.assertIsInstance(self.loop, PyCifLoop)
        self.assertEqual(self.loop.count_rows(), 0)
        self.assertEqual(self.loop.count_columns(), 0)

    def test_add_column(self):
     
        self.loop.add_column("Column1")
        self.loop.add_column("Column2")
        self.assertEqual(self.loop.count_columns(), 2)
        self.assertEqual(self.loop.column_names(), ["Column1", "Column2"])

    def test_add_data_row(self):
  
        self.loop.add_column("Column1")
        self.loop.add_column("Column2")
        self.loop.add_data_row(["Data1", "Data2"])
        self.assertEqual(self.loop.count_rows(), 1)
        self.assertEqual(self.loop.row(0), ["Data1", "Data2"])

    def test_row_out_of_bounds(self):
        
        self.loop.add_column("Column1")
        self.loop.add_data_row(["Data1"])
        self.assertIsNone(self.loop.row(5))

    def test_column_index(self):
   
        self.loop.add_column("Column1")
        self.loop.add_column("Column2")
        self.assertEqual(self.loop.column_index("Column1"), 0)
        self.assertEqual(self.loop.column_index("Column2"), 1)
        self.assertIsNone(self.loop.column_index("Column3"))

    def test_empty_columns_and_rows(self):
      
        self.assertEqual(self.loop.column_names(), [])
        self.assertIsNone(self.loop.row(0))

    def test_add_incomplete_data_row(self):
      
        self.loop.add_column("Column1")
        self.loop.add_column("Column2")
        self.loop.add_data_row(["Data1"])  # Ligne incomplète
        
        self.assertEqual(self.loop.count_rows(), 1)
        self.assertEqual(self.loop.row(0), ["Data1"])

    def test_multiple_rows_and_columns(self):
      
        self.loop.add_column("Column1")
        self.loop.add_column("Column2")
        self.loop.add_data_row(["Data1", "Data2"])
        self.loop.add_data_row(["Data3", "Data4"])
        self.assertEqual(self.loop.count_rows(), 2)
        self.assertEqual(self.loop.count_columns(), 2)
        self.assertEqual(self.loop.row(0), ["Data1", "Data2"])
        self.assertEqual(self.loop.row(1), ["Data3", "Data4"])

if __name__ == "__main__":
    unittest.main()