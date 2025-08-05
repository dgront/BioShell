**Notice**: you must provide the path to the ``taxdump.tar.gz`` file downloaded from the NCBI website
with `-p` option

1. Print taxonomy information for a species given its name:
    ```ignore
    taxonomy -n "Homo sapiens" -p ./data/
    ```
 
2. Print taxonomy information for a given taxid:
    ```ignore
    taxonomy -p ./data/ -t 9606
    ```

3. Taxonomy can be also exported to JSON:
    ```ignore
    taxonomy -p ./data/ -t 9606 --json
    ```

4. Print the full lineage for a given species:
    ```ignore
    taxonomy -p ./data/ -n "Escherichia coli" --lineage
    ```

5. Export the full lineage as JSON:
    ```ignore
    taxonomy -p ./data/ -n "Escherichia coli" --lineage --json
    ```
