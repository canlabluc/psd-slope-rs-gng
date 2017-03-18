"""
This script removes unneeded markers from the EMSE-exported .evt files.
These are the following: '[seg]', '222', '252', '223', '255', as well as
the unneeded headings at the top of the files, such as <bUseQID>
"""

import os
import sys
import glob


def get_filelist(import_path, extension):
    """
    Returns list of file paths from import_path with specified extension.
    """
    filelist = []
    for root, dirs, files in os.walk(import_path):
        filelist += glob.glob(os.path.join(root, '*.' + extension))
        return filelist


importpath = sys.argv[1]
exportpath = sys.argv[2]

files = get_filelist(importpath, 'evt')
for file in files:
    lines = open(file, 'r').readlines()
    lines[1] = ''
    lines[2] = ''
    lines[3] = ''
    for i in range(len(lines)):
        if ('<Name>[seg]</Name>' in lines[i] or
           '<Name>222</Name>' in lines[i] or
           '<Name>223</Name>' in lines[i] or
           '<Name>252</Name>' in lines[i] or
           '<Name>255</Name>' in lines[i]):
           lines[i] = ''
    new = open(exportpath + file.split('/')[-1], 'w')
    new.writelines(lines)
