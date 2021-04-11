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

with open("data/qc/minIds.txt") as file:
    __minimalSet = [id for id in file]

def getID(record):
    for idType in __minimalSet:
        if(idType in record.qualifiers):
            return record.qualifiers[idType][0]
    else: raise Exception(f"Failed to find suitable ID for:\n{record}")
