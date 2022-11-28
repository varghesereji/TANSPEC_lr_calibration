# This code is to sort the filenames and save as org file.

import re
import os


def extract_filenumber(fname):
    """ Function which returns the key to sort SpecTANSPEC filename """
    filenumber = re.search('(\w+?)\.Z\.fits', fname).group(1)
    return filenumber


path = '/home/varghese/Desktop/DP2_TANSPEC/Data/tanspec_slopes/HIP41815'
files_list = os.listdir(path)

print('|-------|-------|')
for files in files_list:
    splitted = files.split('-')
    filenumber = extract_filenumber(splitted[1])
    splitted[1] = filenumber
    print('|', str(splitted[0]), '|', str(splitted[1]), '|')
print('|-------|-------|')
