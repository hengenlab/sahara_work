from sahara_work import Criticality_new as cr
import numpy as np
import json

def run_crit(json_file):
    with open(json_file) as j:
        params = json.load(j)

