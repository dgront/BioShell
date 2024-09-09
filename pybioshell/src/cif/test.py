# import unittest
# from pybioshell.src.cif import cif


# class TestPyCifLoop(unittest.TestCase):

#     def setUp(self):
     
#         self.loop = PyCifLoop()

#     def test_new_instance(self):
 
#         self.assertIsInstance(self.loop, PyCifLoop)
#         self.assertEqual(self.loop.count_rows(), 0)
#         self.assertEqual(self.loop.count_columns(), 0)

#     def test_add_column(self):
     
#         self.loop.add_column("Column1")
#         self.loop.add_column("Column2")
#         self.assertEqual(self.loop.count_columns(), 2)
#         self.assertEqual(self.loop.column_names(), ["Column1", "Column2"])

#     def test_add_data_row(self):
  
#         self.loop.add_column("Column1")
#         self.loop.add_column("Column2")
#         self.loop.add_data_row(["Data1", "Data2"])
#         self.assertEqual(self.loop.count_rows(), 1)
#         self.assertEqual(self.loop.row(0), ["Data1", "Data2"])

#     def test_row_out_of_bounds(self):
        
#         self.loop.add_column("Column1")
#         self.loop.add_data_row(["Data1"])
#         self.assertIsNone(self.loop.row(5))

#     def test_column_index(self):
   
#         self.loop.add_column("Column1")
#         self.loop.add_column("Column2")
#         self.assertEqual(self.loop.column_index("Column1"), 0)
#         self.assertEqual(self.loop.column_index("Column2"), 1)
#         self.assertIsNone(self.loop.column_index("Column3"))

#     def test_empty_columns_and_rows(self):
      
#         self.assertEqual(self.loop.column_names(), [])
#         self.assertIsNone(self.loop.row(0))

#     def test_add_incomplete_data_row(self):
      
#         self.loop.add_column("Column1")
#         self.loop.add_column("Column2")
#         self.loop.add_data_row(["Data1"])  # Ligne incomplète
        
#         self.assertEqual(self.loop.count_rows(), 1)
#         self.assertEqual(self.loop.row(0), ["Data1"])

#     def test_multiple_rows_and_columns(self):
      
#         self.loop.add_column("Column1")
#         self.loop.add_column("Column2")
#         self.loop.add_data_row(["Data1", "Data2"])
#         self.loop.add_data_row(["Data3", "Data4"])
#         self.assertEqual(self.loop.count_rows(), 2)
#         self.assertEqual(self.loop.count_columns(), 2)
#         self.assertEqual(self.loop.row(0), ["Data1", "Data2"])
#         self.assertEqual(self.loop.row(1), ["Data3", "Data4"])

# if __name__ == "__main__":
#     unittest.main()


import pybioshell

def test_cif_loop():
    # Create a new CifLoop instance with column names
    loop = pybioshell.CifLoop(['column1', 'column2'])
    
    # Add a new column
    loop.add_column('column3')
    
    # Add a new row of data
    loop.add_data_row(['data1', 'data2', 'data3'])
    
    # Get the list of rows
    rows = loop.rows()
    assert rows == [['data1', 'data2', 'data3']], f"Expected [['data1', 'data2', 'data3']], got {rows}"
    
    # Get the list of column names
    column_names = loop.column_names()
    assert column_names == ['column1', 'column2', 'column3'], f"Expected ['column1', 'column2', 'column3'], got {column_names}"
    
    # Get the number of rows
    num_rows = loop.count_rows()
    assert num_rows == 1, f"Expected 1, got {num_rows}"
    
    # Get the number of columns
    num_columns = loop.count_columns()
    assert num_columns == 3, f"Expected 3, got {num_columns}"
    
    # Find the index of a column by name
    column_index = loop.column_index('column2')
    assert column_index == 1, f"Expected 1, got {column_index}"
    
    # Check if a column name contains a specific substring
    contains = loop.column_name_contains('col')
    assert contains, "Expected column names to contain 'col'"
    
    # Get an entry mutably
    entry = loop.entry_mut(0, 'column3')
    assert entry == 'data3', f"Expected 'data3', got {entry}"

def test_cif_data():
    # Create a new CifData instance with a name
    data = pybioshell.CifData('test_data')
    
    # Add an item
    data.add_item('key1', 'value1')
    
    # Get an item
    value = data.get_item('key1')
    assert value == 'value1', f"Expected 'value1', got {value}"
    
    # Get all data items
    items = data.data_items()
    assert items == {'key1': 'value1'}, f"Expected {{'key1': 'value1'}}, got {items}"
    
    # Get mutable data items
    items_mut = data.data_items_mut()
    assert items_mut == {'key1': 'value1'}, f"Expected {{'key1': 'value1'}}, got {items_mut}"

def main():
    test_cif_loop()
    test_cif_data()
    print("All tests passed!")

if __name__ == "__main__":
    main()
