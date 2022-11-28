import os
import shutil

from extracting_header_values import file_header_key

path = '/home/varghese/Desktop/DP2_TANSPEC/Data/tanspec_slopes'
files_list = os.listdir(path)


key = ['OBJECT']
object_list = ['HD37725', 'HIP41815']

for files in files_list:
    if files.endswith('.fits'):
        # print(files)
        file_path = os.path.join(path, files)
        header_key = file_header_key(key, file_path)
        for objects in object_list:
            if header_key[key[0]] == objects:
                destination = os.path.join(path, objects, files)
                shutil.copy(file_path, destination)
                print(files, 'go to', destination)
            # else:
            #     if header_key[key[0]] not in object_list[0]:
            #         destination = os.path.join(path, objects, files)
            #         shutil.copy(file_path, destination)
            #         print(files, 'go to', destination)
