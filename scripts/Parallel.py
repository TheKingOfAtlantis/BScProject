
import itertools, more_itertools


import pandas as pd

def getPool(jobCount = None, limit = None, **kwargs):
    """
        Brief: Creates a multiprocess pool to most efficiently tackle the job

        Creates a multiprocess pool to most efficiently tackle the job.
        If the number of jobs is less than the number of cores avaliable to us
        then we allocate the same number of processes as jobs to avoid any penalties
        associated with creating processes which won't do not get used
        If no count is given then we ensure that regardless of the job size enough processes
        will be avaliable to quickly process the job (i.e. the number of processes
        will be the same as the number of cores avaliable)

        Additionally, the user may passes a upper limit on the number of processes. If
        the number of jobs is lower or the number of cores is lower then these will
        define the upper limit, however if these would exceed it then the ceiling is set
        to the given limit

        Param:
            jobCount -- The number of jobs to be completed
            limit    -- User defined limit on number of processes
            **kargs  -- Optional arguments to pass to Pool
    """
    from multiprocessing import Pool
    return Pool(processes = __processesCountHeuristic(jobCount, limit), **kwargs)

def __processesCountHeuristic(jobCount = None, limit = None):
    import os
    count = os.cpu_count() if(jobCount is None) else min(os.cpu_count(), jobCount)
    if(limit is not None): count = min(count, limit)
    return count

def __chunkSizeHeuristic(jobCount, limit = None):
    # Iterables are split into chunks and feed to each process to then process
    # The problem can go two ways:
    #  1) Too small and each process incurs system call costs for retrieving the next batch
    #     BUT jobs more efficinetly distributed (avoiding idle processes)
    #  2) Too large and problems loading the data will occur
    #     BUT avoids issues of receiving many small data packets

    import math
    processes = __processesCountHeuristic(jobCount, limit)

    # If we have <1024 jobs or so per process => use default behaviour
    # 10,000 => With our 40-cores = 40,000 jobs before we use our chunksize
    if(jobCount/processes < 10000): return None
    # We start with batches of 10,000 / 4
    # Gently increase with job size such that we increase it by 10,000 / 4
    # with each order of magitude increase
    else: return math.floor(10000 * math.floor(math.log10(jobCount/processes) - 3) / 4)


def loadParallel(callable, param, count = None, **tqdmParam):
    from tqdm.auto import tqdm
    # Create a pool of process which are used to run operations on each file
    with getPool(count) as pool:
        # Check if we want to display the progress (using tqdm)
        result = list(tqdm(
            pool.imap(callable, param, __chunkSizeHeuristic(count)),
            total=count, # Pass number of files to process
            **tqdmParam
        ))
        return [x for x in result if x is not None] # Remove any None values

def __concatPreprocessor(x):
    data, preprocessor, args, kwargs = x
    return pd.concat(preprocessor(data), *args, **kwargs)

def __concat(x):
    data, args, kwargs = x
    return pd.concat(data, *args, **kwargs)

def concat(data, preprocessor = None, *args, **kwargs):
    '''
        data - Iterable with the data to be concat into a DataFrame
        preprocessor - Initial concat function: Allows for concat with non-DataFrame data
                  before working with list of dataframes
    '''
    with getPool() as pool:
        if(preprocessor != None):
            data = pool.map(
                __concatPreprocessor, zip(
                    more_itertools.chunked(data, 50),
                    itertools.repeat(preprocessor),
                    itertools.repeat(args),
                    itertools.repeat(kwargs)
                )
            )

        while(len(data) > 1):
            data = pool.map(
                __concat, zip(
                    more_itertools.chunked(data, 50),
                    itertools.repeat(args),
                    itertools.repeat(kwargs)
                )
            )

        return data[0]
