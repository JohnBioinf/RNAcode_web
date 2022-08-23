from ete3 import NCBITaxa
import os


print("Update ete3 taxonomy database")

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()

if os.path.isfile("./taxdump.tar.gz"):
    os.remove("./taxdump.tar.gz")
