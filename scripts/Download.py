
import requests
from Parallel import *
from tqdm.auto import tqdm

import time
errPath = f"download-err-{time.time()}.log"

def __request(method, url, **tqdmParams):
    response   = method(url, stream=True)
    block_size = 1024
    file_size  = int(response.headers.get('Content-Length', 0))

    progress = tqdm(**tqdmParams)

    result = b''
    for data in response.iter_content(block_size):
        progress.update(len(data))
        result += data

    return result

def __requestFile(method, url, file, **tqdmParams):
    try:
        response   = method(url, stream=True)
        block_size = 1024
        file_size  = int(response.headers.get('Content-Length', 0))

        progress = tqdm(
            unit       = "iB",
            total      = file_size,
            unit_scale = True,
            **tqdmParams
        )

        with open(file, "wb") as file:
            for data in response.iter_content(block_size):
                progress.update(len(data))
                file.write(data)
    except Exception as e:
        with open(errPath, "a") as err:
            err.write(f"{url} failed with {e}")
def __processRequest(x):
    (method, url, pos, tqdmParams) = x
    return __request(method, url, pos, **tqdmParams)
def __processFileRequest(x):
    (method, url, file, pos, tqdmParams) = x

    if("desc" not in tqdmParams):
        tqdmParams["desc"] = f"Downloading file #{pos + 1}"

    return __requestFile(method, url, file, position=pos + 1, **tqdmParams)

def get(urls, limit = None, **tqdmParams):
    import requests
    """
        Brief: Performs a GET request with a progress bar for list of urls

        Parameters:
            url        -- Url for the GET request
            tqdmParams -- Parameters to customise tqdm
    """
    # If we were given a string then we only have 1 url to get
    # So we don't need to parallelise it
    if(isinstance(urls, str)):
        return __request(requests.get, urls, **tqdmParams)


    with getPool(len(urls), limit) as pool:
        return list(tqdm(
            pool.imap(
                __processRequest, zip(
                    itertools.repeat(requests.get),
                    urls,
                    range(len(urls)),
                    itertools.repeat(tqdmParams)
                )
            ),
            desc     = description,
            total    = len(urls),
            position = 0
        ))

def getFile(urls, paths, limit = None, description = None, **tqdmParams):
    import requests
    """
        Brief: Performs a GET request with a progress bar for list of urls

        Parameters:
            urls       -- Urls for the GET request
            tqdmParams -- Parameters to customise tqdm
    """
    # If we were given a string then we only have 1 url to get
    # So we don't need to parallelise it
    if(isinstance(urls, str)):
        __requestFile(requests.get, urls, paths, **tqdmParams)
    else:
        with getPool(len(urls), limit) as pool: list(tqdm(
            pool.imap(
                __processFileRequest, zip(
                    itertools.repeat(requests.get),
                    urls, paths,
                    range(len(urls)),
                    itertools.repeat(tqdmParams)
                )
            ),
            desc     = description,
            total    = len(urls),
            position = 0
        ))

def post(urls, limit = None, description = None, **tqdmParams):

    """
        Brief: Performs a GET request with a progress bar

        Parameters:
            url        -- Url for the GET request
            tqdmParams -- Parameters to customise tqdm
    """
    return __request(requests.post, url, tqdmParams)
