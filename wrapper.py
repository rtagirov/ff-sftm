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

#D = [26.5, 30.9, 41.3, 47.0]
#D = [26.5, 47.0]
D = [60, 80, 100]

if ar == 's':

    for d in D:

        fname = 'C22_' + str(d)

        if not os.path.isfile('./inp/' + fname):

            auxsys.abort('Spot input file not found.')

        os.system('python spot.py --np ' + nproc + ' --inp ' + fname)

elif ar == 'f':

#    B_sat = np.arange(200., 550., 50.)
#    B_sat = np.array([200.0, 500.0,])
    B_sat = np.array([400.0, 450.0, 500.0])

#    B_spot = np.array([800.0, 1200.0,])
    B_spot = np.array([1000.0])

#    psets = list(itertools.product(D, B_sat, B_spot))
    psets = list(itertools.product(D, B_sat, B_spot))

    server = socket.gethostname()

    if server == 'pulpo':

#        psets = psets[0 : 7]
        psets = psets[0 : 2]

    elif server == 'mojo':

        psets = psets[2 : 4]

    elif server == 'helios1':

        psets = psets[4 : 6]

    elif server == 'helios2':

        psets = psets[6 : 8]

    elif server == 'ph-rtagirov':

        psets = psets[8]

    else:

        auxsys.abort('Server not recognized.')

    for s in psets:

        os.system('python facula.py --np ' + nproc + \
                                  ' --D ' + str(s[0]) + \
                                  ' --Bsat ' + str(s[1]) + \
                                  ' --Bspot ' + str(s[2]))

else:

    auxsys.abort('Active region type not recognized.')

