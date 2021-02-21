import os, io, itertools         # Various useful and necessary imports
import tqdm                      # Used to track progress
from zipfile import ZipFile      # Used to access the zipfile
from multiprocessing import Pool # Used to perform access in parallel fashion

if __name__ == "__main__":
    print("Do not directly call this file! Instead import it where you need it")
    exit(-1)

def __loadZip(x):
    # Each process calls this function
    # In here we get the file and pass it to the operation

    # Expand the parameters given to this function:
    # 1) Path to the zip file
    # 2) Index position of the file within the zip file
    # 3) Operation to be performed on the file
    zipPath, pos, op = x

    with ZipFile(zipPath) as zipFile:  # Open zip file
        path = zipFile.namelist()[pos] # Get path of file in zip file given index

        info = zipFile.getinfo(path)
        if info.is_dir(): # Check if path is to a directory
            return        # If so then return without doing anything

        file = io.TextIOWrapper(zipFile.open(path)) # Get the file and since given as bytes => Convert to Text
        return op(file)                             # Perform operation on file

def loadZip(path, operation, progress = True):
    with ZipFile(path) as zipFile:
        # Load the zip file just so we can get the number of files present
        count = len(zipFile.namelist())

    with Pool(min(
        os.cpu_count(),
        count
    )) as pool: # Lets not use up more than we need
        # Create a pool of process which are used to run operations on each file in zip file

        param = zip(        # Since we can only pass a single argument: tuple them together
            itertools.repeat(path),      # Use a repeat iterator so every process has path to the zip file
            range(0, count),             # Then pass each file an index to the file itself
            itertools.repeat(operation)  # Then pass the operation to be run to each process
        )

        # Check if we want to display the progress
        if progress: return list(tqdm.tqdm(
            pool.imap(__loadZip, param),
            total=count # Pass total count to tqdm but remove one as we ignore the root directory
        ))
        else: return pool.map(__loadZip, param)
