import pickle
import os
import datetime
import hyperopt
from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
import time
import pandas as pd
import numpy as np


def getTrials(filename,restart = False ):
    if os.path.isfile(filename) and not restart:
        trials = pickle.load(open('deng_trials.pkl', 'rb'))
        print(f'Reload trials size :{len(trials)}')
    else:
        print('Creating new trials...')
        trials = Trials()
        
    return trials

def getResultsAsDf(trials, space):
    df= pd.DataFrame(columns=space.keys())
    results = []
    for t in trials.trials:
        # the values in vals are put inside arrays, convert this structure
        # to a json of values
        vals = {k:v[0] for k,v in list(t['misc']['vals'].items())}
        df.loc[df.shape[0]] = hyperopt.space_eval(space, vals)
        results.append(t['result']['loss'] if 'loss' in t['result'].keys() else np.nan )
    df['result'] = results
    return df
