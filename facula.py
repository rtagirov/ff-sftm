import importlib
import itertools
import auxfunc
import sys
import glob
import re
import sharedmem

import numpy as np

from tqdm import tqdm
from itertools import repeat

importlib.reload(auxfunc)

conv = np.pi / 180.0

# set position of observer and Bsat
x_c = 0
y_c = 0
B_sat = 484.0
B_spot = 1000.0

mu_low = [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075, 0.0]
mu_up = [1.0, 0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.075]

start = 170049

nproc = 4

if len(sys.argv) == 2:

    nproc = int(sys.argv[1])

norm = 90 * 90 * 4 / np.pi**2 * np.pi

magnetograms = sorted(glob.glob('./mag/CalcMagnetogram.2000.*'))

spot_mask = np.load('spot_mask.npy').item()

f = open('ff_fac.out','w')

px = list(itertools.product(range(180), range(360)))

for _, mag in enumerate(tqdm(magnetograms, \
                             ncols = auxfunc.term_width(), \
                             desc = 'Masking faculae', \
                             position = 0)):

    name = re.findall('2000.(\d+)', mag)

    data = np.loadtxt(mag)

    date = int(name[0])

    visibility = sharedmem.empty(1, dtype = 'f8')

    r = sharedmem.empty(11, dtype = 'f8')

    spot_x = np.concatenate((spot_mask[date]['xp'], spot_mask[date]['xn']))
    spot_y = np.concatenate((spot_mask[date]['yp'], spot_mask[date]['yn']))

    with sharedmem.MapReduce(np = nproc) as pool:
    
        def px_contrib(px):

            i = px[0]

            j = px[1]

            x_min = j
            x_max = j + 1

            y_min = i
            y_max = i + 1

            spot = spot_x[np.where((spot_x >= x_min) & (spot_x < x_max) & (spot_y >= y_min) & (spot_y < y_max))]

            B = abs(data[i][j])

            ff = 0

            if np.shape(spot) != (0, ):

                helper = B - B_spot * len(spot) * 0.1 * 0.1

                if helper > 0:

                    ff = (1 - len(spot) * 0.1 * 0.1) * helper / B_sat

            if np.shape(spot) == (0, ) and B < B_sat:

                ff = B / B_sat

            if np.shape(spot) == (0, ) and B >= B_sat:

                ff = 1.0

            x_rot = (j + 13.28 * (date - start)) % 359

            x_pos = 180.0 - x_rot

            y_pos = 90.0 - i

            delta_lambda = abs(x_pos - x_c)

            distance = np.arccos(np.sin(y_c * conv) * np.sin(y_pos * conv) + np.cos(y_c * conv) * \
                                 np.cos(y_pos * conv) * np.cos(delta_lambda * conv)) / conv

            vis = np.cos(distance * conv)

            idx = np.where((vis > mu_low) & (vis <= mu_up))

            r[idx] += ff * vis * np.cos(y_pos * conv)

            if distance <= 90.0:

                visibility[0] += ff * np.cos(distance * conv) * np.cos(y_pos * conv)

        pool.map(px_contrib, px)

#        pool.close()
#        pool.join()

#    p = Pool(processes = nproc)

#    p.starmap(scan_mag, zip(px, repeat(spot_x), repeat(spot_y), repeat(data), repeat(date)))

#    p.close()
#    p.join()

    r /= norm

    visibility /= norm

    f.write('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n' \
            %(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10], \
              date, sum(r), visibility))

f.close()

