from bioshell.taxonomy import Taxonomy

taxonomy = Taxonomy.load_from_tar_gz("test_files/test_taxdump.tar.gz")
human_taxid = taxonomy.taxid("Homo sapiens")
assert human_taxid == 9606

human_order = taxonomy.rank(human_taxid, "Order")
assert human_order.name == "Primates"
assert str(human_order.rank) == "Order"

human_family = taxonomy.rank(human_taxid, "Family")
assert human_family.name == "Hominidae"
assert str(human_family.rank) == "Family"

lineage = taxonomy.lineage(human_taxid)
assert str(lineage[0].rank) == "Other"
assert str(lineage[-1].rank) == "Species"
assert len(lineage) == 32
