
import itertools, more_itertools

from multiprocessing import Pool
from tqdm.auto import tqdm

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
    import os
    processes = os.cpu_count() if(jobCount is None) else min(
        os.cpu_count(),
        jobCount
    )
    if(limit is not None):
        processes = min(
            processes,
            limit
        )

    return Pool(processes = processes, **kwargs)

def loadParallel(callable, param, count = None, **tqdmParam):
    # Create a pool of process which are used to run operations on each file
    with getPool(count) as pool:
        # Check if we want to display the progress (using tqdm)
        result = list(tqdm(
            pool.imap(callable, param),
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
