from enum import Enum
from tqdm.auto import tqdm, trange

import math

class EntrezXML:
    class Action(Enum):
        search  = "esearch"
        summary = "esummary"
        info    = "einfo"
        fetch   = "efetch"
        post    = "epost"

    def __init__(self, email, api = None):
        self._email = email
        self._api   = api

    lastRequest = 0    # When we last requested data
    delay       = 0.37 # Can send 3 requests per sec w/o API (so delay for little)

    @staticmethod
    def __runAction(action, progress, **params):
        """
            Sends the request to Entrez

        Parameters:
            action -- Action by which to process the request
            params -- Parameters to pass to Entrez
        """

        import requests
        import time
        # Can only make 3 requests per second
        # So if we are making requests too quickly we need to slow down

        if(EntrezXML.lastRequest - time.time() > EntrezXML.delay):
            wait = (EntrezXML.lastRequest + EntrezXML.delay) - time.time()
            time.sleep(wait)
        EntrezXML.lastRequest = time.time()

        base      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        actionUrl = base + f"{action.value}.fcgi"

        if(action == EntrezXML.Action.post):
            return requests.post(actionUrl, data = params).text

        if(progress):
            response   = requests.get(actionUrl, params=params, stream=True)
            block_size = 1024
            file_size  = int(response.headers.get('Content-Length', 0))

            progress = tqdm(
                desc       = f"Entrez {action.name}",
                total      = file_size if (file_size > 0) else None,
                unit       = 'iB',
                unit_scale = True
            )

            result = b''
            for data in response.iter_content(block_size):
                progress.update(len(data))
                result += data

            return result.decode("utf-8")
        else: return requests.get(actionUrl, params=params).text

    @staticmethod
    def __parseNode(node):
        ''' Used to recursively parse the XML DOM '''
        if(len(list(node)) == 0):
            try:    return {node.tag: int(node.text)}
            except: return {node.tag: node.text}
        else:
            children = [EntrezXML.__parseNode(child) for child in node]
            attrib   = [] if node.attrib == {} else [node.attrib]

            return { node.tag: attrib + children }

    @staticmethod
    def __toXML(response):
        import xml.etree.ElementTree as xml
        return xml.fromstring(response)

    @staticmethod
    def __parse(response):
        ''' Converts the result into a JSON-like python object'''
        response_xml = EntrezXML.__toXML(response)
        as_dict = EntrezXML.__parseNode(response_xml)
        return as_dict

    @staticmethod
    def __process(action, result):
        ''' Performs additional processing on parsed results'''
        import collections
        if(action == EntrezXML.Action.search):
            [(rootName, values)] = result.items()
            return dict(collections.ChainMap(*values))
        else: return result

    def get(self, action, progress = True, output = "python", **params):
        """
        Sends a request to performa specific action using the Entrez Utility APIs

        Paramaters:
            action   -- The action to be used for this request, these actions include:
                            - Entrez.Action.search:  Sends a eSearch GET request
                            - Entrez.Action.summary: Sends a eSummary GET request
                            - Entrez.Action.info:    Sends a eInfo GET request
                            - Entrez.Action.fetch:   Sends a eFetch GET request
                            - Entrez.Action.post:    Sends a ePost POST request
            progress -- Whether or not to display the tqdm progress bar
                        while downloading the results
            output   -- Determines how the result is outputted:
                            - raw    => Raw result received are returned
                            - native => Skips the post-processing of the results into a JSON-like python object
                            - python => Creates a nested JSON-like python object
        """

        # Add additional parameters required by the Entrez API
        if("tool" not in params):
            params["tool"] = "python"
        if("email" not in params and self._email is not None):
            params["email"] = self._email
        if("api_key" not in params and self._api is not None):
            params["api_key"] = self._api

        if("id" in params):
            try: params["id"] = ",".join(iter(params["id"]))
            except: pass

        response = self.__runAction(action, progress, **params)

        if(output == "raw"): return response
        elif(output == "native"): return self.__toXML(response)
        elif(output == "python"):
            xml = self.__parse(response)
            return self.__process(action, xml)

    def getAll(self, action, count, progress = True, output = "python", **params):
        blockSize = {
            EntrezXML.Action.search:  100000,
            EntrezXML.Action.summary: 10000,
            EntrezXML.Action.fetch:   10000
        }

        return list(
            self.get(
                action   = action,
                progress = progress,
                output   = output,
                retstart = i * blockSize[action],
                retmax   = blockSize[action],
                **params
            ) for i in trange(math.ceil(count/blockSize[action]))
        )

Entrez = EntrezXML("ss2980@bath.ac.uk")
