import os
import sys
import auxsys
import importlib

import numpy as np

importlib.reload(auxsys)

def get_args(args):

    nproc = '4'

    for i, arg in enumerate(args):

        if arg == '--np':

            nproc = args[i + 1]

    return nproc

nproc = get_args(sys.argv[1:])

decay_rates = [26.5, 30.9, 41.3, 47]

B_sat = np.arange(200., 550., 50.)

for rate in decay_rates:

    filename = 'C22_' + str(rate)

    if os.path.isfile('./inp/' + filename):

        os.system('python spot.py --np ' + nproc + ' --inp ' + filename)

    else:

        auxsys.abort('Spot input file not found.')

    for B in B_sat:

        os.system('python facula.py --np ' + nproc + ' --D ' + str(rate) + ' --Bsat ' + str(B))
