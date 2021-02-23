import itertools
from zipfile import ZipFile # Used to access the zipfile

if __name__ == "__main__":
    print("Do not directly call this file! Instead import it where you need it")
    exit(-1)

def __loadParallel(callable, param, count, progress = True):
    import os
    import tqdm # Used to track progress
    from multiprocessing import Pool # Used to perform access in parallel fashion

    # Create a pool of process which are used to run operations on each file
    with Pool(min( # Lets not use up more than we need
        os.cpu_count(),
        count
    )) as pool:

        # Check if we want to display the progress (using tqdm)
        if progress: result = list(tqdm.tqdm(
            pool.imap(callable, param),
            total=count # Pass number of files to process
        ))
        else: result = pool.map(callable, param)
        return [x for x in result if x is not None] # Remove any None values


def __loadZip(x):
    import io
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
        count = len(zipFile.namelist())

    # Load the zip file just so we can get the number of files present
    return __loadParallel(
        __loadZip,
        zip( # Since we can only pass a single argument: zip them together
            itertools.repeat(path),      # Use a repeat iterator so every process has path to the zip file
            range(0, count),             # Then pass each file an index to the file itself
            itertools.repeat(operation)  # Then pass the operation to be run to each process
        ), count
    )

def __loadGlob(x):
    # Each process calls this function
    # In here we get the file and pass it to the operation
    path, op = x

    with open(path) as file: # Open zip file
        return op(file)      # Perform operation on file
def loadGlob(pattern, operation, progress = True):
    import glob


    paths = glob.glob(pattern)
    return __loadParallel(
        __loadGlob,
        zip( # Since we can only pass a single argument: zip them together
            paths, # List of paths to process
            itertools.repeat(operation)  # Operation to be run in each process/file
        ), len(paths)
    )
