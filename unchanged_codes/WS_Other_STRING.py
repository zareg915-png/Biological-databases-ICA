#!/usr/bin/env python3

# ########################################################################
#   This is a custom script to identify interactions from String        #
#   Notice how the table structure is formatted to allow the later join #
#    Simon Tomlinson  Nov 2026                                          #
# ########################################################################

import sys
import csv
import time
from typing import Dict, List, Optional, Any, Tuple
import json

import requests

BASE_STRING_QUERY="https://string-db.org/api/json/network?identifiers="
SPECIES = "9606"

#Generated from ELM
def http_get(
        url: str,
        params: Optional[Dict[str, Any]] = None,
        timeout: int = 60,
) -> Optional[requests.Response]:

    try:
        resp = requests.get(url, params=params, timeout=timeout)
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] GET {url} failed: {e}", file=sys.stderr)
        return None

    if resp.status_code != 200:
        print(
            f"[ERROR] GET {url} returned status {resp.status_code}: "
            f"{resp.text[:200]}",
            file=sys.stderr,
        )
        return None

    return resp

def build_string_query(EXAMPLE_HUMAN_INPUTS,species =9606):
    myquery = BASE_STRING_QUERY

    for idx in range(0,len(EXAMPLE_HUMAN_INPUTS)-1):
        myquery =myquery + EXAMPLE_HUMAN_INPUTS[idx]
        myquery = myquery + "%0D"
    #add the last element
    myquery = myquery + EXAMPLE_HUMAN_INPUTS[len(EXAMPLE_HUMAN_INPUTS)-1]

    # add the species
    myquery = myquery + "&species=" + SPECIES
    return myquery

def connect_string(EXAMPLE_HUMAN_INPUTS,string_output):
    qry = build_string_query(EXAMPLE_HUMAN_INPUTS)

    res=http_get(qry)
    #if resp is None:
      #      continue

    lines = res.text.strip()
    data =  json.loads(lines)
    print("completed connection")

    #first create a new dictionary for each gene in the input
    final_set = {}
    for ele in EXAMPLE_HUMAN_INPUTS:
        final_set[ele] = ""

    #iterate the data from the query and store the results in the simplified dictionary
    for aline in data:
        geneA =aline['preferredName_A']
        geneB =aline['preferredName_B']

        #stash the results
        if geneA in final_set.keys():
            final_set[geneA] = final_set[geneA]  + geneB + "| "
        if geneB in final_set.keys():
            final_set[geneB] = final_set[geneB]  + geneA + "| "

    #output set as a file...
    table_output_file = string_output + ".tsv"

    print("[INFO] Writing interaction table to " + table_output_file)
    with open(table_output_file, 'w') as f:
        f.write("Gene\tInteractions\n")#the header

        for k in final_set.keys():
            f.write(k + "\t" + final_set[k] + "\n")

    return final_set

#example function to drive the code...
#connect_string(gene_list,"tesoutput")
