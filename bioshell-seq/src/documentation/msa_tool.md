
1) **Basic information about an MSA**
   - Print basic information about the provided multiple sequence alignment:
   ```ignore
   msa_tool  -s test_files/4Fe-4S-example.sto --info
    ```
   - As above, but creates a tab-separated table for a whole folder of CIF files:
    ```ignore
    for i in *.cif; do pdb_tool -i $i --info-table; done > table.dat
    ```
   - Print only the following columns in the table: id, resolution, methods
    ```ignore
    pdb_tool -i 1ehe.cif --info-table id resolution methods
    ```
   
2) **Convert between data formats**

   - Convert the MSA in Stockholm format to fasta
   ```ignore
   msa_tool -s test_files/4Fe-4S-example.sto -o out.fasta
   ```
   - Convert the MSA in fasta format to Stockholm, print on the screen
   ```ignore
   msa_tool -f alignment.fasta --out-sto 
   ```

