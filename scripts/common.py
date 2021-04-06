# Set of common subroutines
#
# Most of these are to perform common tasks in parallel to maximise usage
# of the computers cores to reduce the duration of the tasks
#

if __name__ == "__main__":
    print("Do not directly call this file! Instead import it where you need it")
    exit(-1)

import Parallel
import Filesystem
import Download

def getID(record):
    idtype = None

    if("protein_id" in record.qualifiers):  idtype = "protein_id"
    elif("locus_tag" in record.qualifiers): idtype = "locus_tag"
    else: raise Exception(f"Failed to find suitable ID for:\n{record}")

    return record.qualifiers[idtype][0]
