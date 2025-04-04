1) Basic information about a deposit
   - Print basic information from a file header:
   ```ignore
   pdb_tool -i 1akd.cif --info
    ```
   - As above, but creates a tab-separated table for a whole folder of CIF files:
    ```ignore
    for i in *.cif; do pdb_tool -i $i --info-table; done > table.dat
    ```
   - Print only the following columns in the table: id, resolution, methods
    ```ignore
    pdb_tool -i 1ehe.cif --info-table id resolution methods
    ```
   
2) Entities
   - List entities in a CIF file:
    ```ignore
    pdb_tool -i 1ehe.cif --entities
    ```
   - print the sequence of a selected entity as defined by the coordinates:
   ```ignore
   pdb_tool -i 1ehe.cif --select-entity 2 --out-fasta
   ```
   - print the full sequence of each entity as defined by the CIF file:
   ```ignore
   pdb_tool -i 1ehe.cif --entities --entity-sequence
   ```
   
3) Convert between data formats

   - Attempt to convert the whole CIF file to PDB:
   ```ignore
   pdb_tool -i file.cif -o file.pdb
   ```
   - Save only a selected chain to PDB:
   ```ignore
   pdb_tool -i file.cif --select-chain AA -o file.pdb
   ```

