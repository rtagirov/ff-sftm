import os
import sys
import auxsys
import importlib
import itertools
import socket

import numpy as np

importlib.reload(auxsys)

def get_args(args):

    nproc = '4'

    ar = 'spots'

    for i, arg in enumerate(args):

        if arg == '--np':

            nproc = args[i + 1]

        if arg == '--ar':

            ar = args[i + 1]

    return nproc, ar

nproc, ar = get_args(sys.argv[1:])

D = [26.5, 30.9, 41.3, 47.0]

if ar == 's':

    for d in D:

        fname = 'C22_' + str(d)

        if not os.path.isfile('./inp/' + fname):

            auxsys.abort('Spot input file not found.')

        os.system('python spot.py --np ' + nproc + ' --inp ' + fname)

elif ar == 'f':

    B_sat = np.arange(200., 550., 50.)

    pairs = list(itertools.product(D, B_sat))

    server = socket.gethostname()

    if server == 'pulpo':

        pairs = pairs[0 : 7]

    elif server == 'mojo':

        pairs = pairs[7 : 14]

    elif server == 'helios1':

        pairs = pairs[14 : 21]

    elif server == 'helios2':

        pairs = pairs[21 : 28]

    else:

        auxsys.abort('Server not recognized.')

    for p in pairs:

        os.system('python facula.py --np ' + nproc + ' --D ' + str(p[0]) + ' --Bsat ' + str(p[1]))

else:

    auxsys.abort('Active region type not recognized.')

