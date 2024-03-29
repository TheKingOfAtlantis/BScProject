
from Parallel import *

def mkdir(paths):
    import pathlib
    if(isinstance(paths, str)):
        pathlib.Path(paths).mkdir(parents=True, exist_ok=True)
    else:
        for path in paths:
            pathlib.Path(path).mkdir(parents=True, exist_ok=True)

def __loadGlob(x):
    # Each process calls this function
    # In here we get the file and pass it to the operation
    path, op, mode, args, kwargs = x

    with open(path, mode) as file: # Open zip file
        return op(file, *args, **kwargs)

def loadGlob(pattern, callable, mode = "r", desc = None, *args, **kwargs):
    """
        Loads various files in parallel and runs a given subroutine with the
        file as the first parameter

        Parameters:
            pattern:     -- Glob pattern to match files for processing
            callable     -- Callable subroutine to run on each file
            desc         -- Description to be used by tqdm
            args, kwargs -- Extra fixed arguments passed to every call
    """
    import glob

    paths = glob.glob(pattern)
    return loadParallel(
        __loadGlob,
        zip( # Since we can only pass a single argument: zip them together
            paths,                      # List of paths to process
            itertools.repeat(callable), # Callable to be run in each process/file
            itertools.repeat(mode),     # Mode by which the file shall be accessed
            itertools.repeat(args),     # Extra parameters to pass to operation
            itertools.repeat(kwargs)    # Extra parameters to pass to operation
        ), len(paths),
        desc = desc
    )
