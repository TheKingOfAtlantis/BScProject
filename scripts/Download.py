
import requests
from Parallel import *
from tqdm.auto import tqdm

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
    response   = method(url, stream=True)
    block_size = 1024
    file_size  = int(response.headers.get('Content-Length', 0))

    progress = tqdm(**tqdmParams)

    with open(file, "wb") as file:
        for data in response.iter_content(block_size):
            progress.update(len(data))
            file.write(data)

def __processRequest(x):
    (method, url, tqdmParams) = x
    return __request(method, url, **tqdmParams)
def __processFileRequest(x):
    (method, url, file, tqdmParams) = x
    return __requestFile(method, url, **tqdmParams)

def get(urls, **tqdmParams):
    import requests
    """
        Brief: Performs a GET request with a progress bar for list of urls

        Parameters:
            url        -- Url for the GET request
            tqdmParams -- Parameters to customise tqdm
    """
    if(isinstance(urls, str)):
        return __request(requests.get, urls, **tqdmParams)

    with getPool(len(urls)) as pool:
        return list(tqdm(
            pool.imap(
                __processRequest, zip(
                    itertools.repeat(requests.get),
                    urls,
                    itertools.repeat(tqdmParams)
                )
            ),
            total = len(urls),
            position = 0
        ))

def getFile(urls, paths, **tqdmParams):
    import requests
    """
        Brief: Performs a GET request with a progress bar for list of urls

        Parameters:
            url        -- Url for the GET request
            tqdmParams -- Parameters to customise tqdm
    """
    if(isinstance(urls, str)):
        return __request(requests.get, urls, **tqdmParams)

    with getPool(len(urls)) as pool:
        return list(tqdm(
            pool.imap(
                __processFileRequest, zip(
                    itertools.repeat(requests.get),
                    urls,
                    paths,
                    itertools.repeat(tqdmParams)
                )
            ),
            total = len(urls),
            position = 0
        ))

def post(urls, **tqdmParams):

    """
        Brief: Performs a GET request with a progress bar

        Parameters:
            url        -- Url for the GET request
            tqdmParams -- Parameters to customise tqdm
    """
    return __request(requests.post, url, tqdmParams)
