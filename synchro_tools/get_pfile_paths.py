'''
This is the new way to generate a list of pfiles paths (paths.csv) that will be scp to dev in order to generate the symlinks.

For example:

1. On /home/cyrf0006/data/dev_database/pfiles
$ python get_pfile_paths.py (this script)

2. $ scp paths.csv cyrf@nl-bpo-dev.ent.dfo-mpo.ca:.

3. connect to dev and do:
$ python symlink_pfiles_onDev.py

4. Exit Dev and go back to /home/cyrf0006/data/dev_database/pfiles:
$ scp cyrf@nl-bpo-dev.ent.dfo-mpo.ca:symlinks_pfiles/* .

5. in /home/cyrf0006/data/dev_database/, generate the list:
$ ls pfiles/*.p2019 > 2019.list

6. Generate the yearly netCDF with:
$ python ~/github/AZMP-NL/python-toolbox/synchro_tools/yearly_netcdf_gen.py 
(the bin size can be specified in this script)


This script combines parts of  orignal extract_dev_pfiles.py 

Frederic.Cyr@dfo-mpo.gc.ca -January 2020

'''

import os
import pandas as pd

host = 'cyrf@nl-bpo-dev.ent.dfo-mpo.ca' 
path = '/data/seabird/profile_data'
dest_folder = '/home/cyrf0006/data/dev_database'
list_files = 'pfiles_and_time.list'
my_password = 'anai123!'

# For 2019:
#get_list = "ssh -q cyrf@nl-bpo-dev.ent.dfo-mpo.ca 'find /data/seabird/* -type f -name *.p2019 -exec ls -l {} \;' > pfiles_and_time.list"
get_list_command = "ssh -q " + host + " 'find " + path + "/* -type f -name *.p2019 -exec ls -l {} \;' > " + list_files
os.system(get_list_command)

## ----  First Part ---- ##

d = []
with open(list_files) as fp:  
   for idx, line in enumerate(fp):
        line = ' '.join(line.split())
        print(idx, line)
        n = line.strip().split(' ')[5:8]

        if ':' in n[2]:
            timestamp = pd.to_datetime(n[1]+'-'+n[0]+'-2019 ' +n[2], format='%d-%b-%Y %H:%M') # <--------- this is weak!!!
        else:
            timestamp = pd.to_datetime(n[1]+'-'+n[0]+'-'+n[2])
        
        d.append({'file_name':line.strip().split('/')[-1],\
        'last_update': timestamp,\
        'full_path':line.strip().split(' ')[-1]})
                   
df = pd.DataFrame(d)
df.to_pickle('file_names_all.pkl')

df = df.set_index('last_update')
df  = df.sort_index()
df = df.drop_duplicates(subset=['file_name'], keep='last', inplace=False)
df.to_pickle('file_names_unique.pkl')

# Write paths.csv (or path2019.csv)
path_file = 'paths.csv'
df['full_path'].to_csv(path_file, index=False, header=False, sep='\n')

