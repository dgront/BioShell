1) Basic information about a deposit
   - Print basic information from a file header:
   ```
   pdb_tool -i 1akd.cif --info
    ```
   - As above, but creates a tab-separated table for a whole folder of CIF files:
    ```
   for i in *.cif; do pdb_tool -i $i --info-table; done > table.dat
    ```
2) Convert between data formats

   - Attempt to convert the whole CIF file to PDB:
   ```
   pdb_tool -i file.cif -o file.pdb
   ```
   - Save only a selected chain to PDB:
   ```
   pdb_tool -i file.cif --select-chain AA -o file.pdb
   ```

