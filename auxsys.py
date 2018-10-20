import os

from sys import exit

import time

def abort(message = 'Abort.'):

    print('\n' + message)

    exit()

def clean_dir(dir_path, mode = 'verbose'):

    if len(dir_path) == 0: abort('sysaux.clean_dir: directory name has to be provided. Abort.')

    if dir_path == '/':    abort('sysaux.clean_dir: directory can not be root. Abort.')

    if dir_path[len(dir_path) - 1] != '/': dir_path = dir_path + '/'

    if not os.path.exists(dir_path): abort('sysaux.clean_dir: directory ' + dir_path + ' does not exist. Abort.')

    if not os.path.isdir(dir_path):  abort('sysaux.clean_dir: ' + dir_path + ' is not a directory. Abort.')

    if mode != 'verbose' and mode != 'noverbose': abort('sysaux.clean_dir: mode is not recognized. Abort.')

    dir_files = os.listdir(dir_path)

    if mode == 'verbose':

       print('')

       for fname in dir_files:

           os.system('rm ' + dir_path + fname)

           print('sysaux.clean_dir: removed ' + fname + ' in ' + dir_path)

       print('')

    if mode == 'noverbose':

       for fname in dir_files: os.system('rm ' + dir_path + fname)
