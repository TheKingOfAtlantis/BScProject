import os, io, itertools         # Various useful and necessary imports
import tqdm                      # Used to track progress
from zipfile import ZipFile      # Used to access the zipfile
from multiprocessing import Pool # Used to perform access in parallel fashion

if __name__ == "__main__":
    print("Do not directly call this file! Instead import it where you need it")
    exit(-1)

def __loadZip(x):
    zipPath, pos, op = x
    with ZipFile(zipPath) as zipFile:
        path = zipFile.namelist()[pos]
        file = io.TextIOWrapper(zipFile.open(path))
        return op(file)

def loadZip(path, operation, progress = True):
    with ZipFile(path) as zipFile:
        count = len(zipFile.namelist())
    with Pool(os.cpu_count()) as pool:
        if progress: return list(tqdm.tqdm(
            pool.imap(__loadZip, zip(
                itertools.repeat(path),
                range(1, count),
                itertools.repeat(operation)
            )),
            total=count - 1
        ))
        else: return pool.map(__loadZip, zip(
            itertools.repeat(path),
            range(1, count),
            itertools.repeat(operation)
        ))
