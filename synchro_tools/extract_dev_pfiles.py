'''
Original script that generates a list of file from Dev and copy them on local computer. I decided not to use this script because the files are copied one by one which is super slow. I decided finally to create symbolink links and scp them to my local machine with 'symlink_pfiles_onDev.py'

Frederic.Cyr@dfo-mpo.gc.ca - April 2018

'''

import os
import pandas as pd

host = 'cyrf@nl-bpo-dev.ent.dfo-mpo.ca' 
path = '/data/seabird'
dest_folder = '/home/cyrf0006/data/dev_database'
list_files = 'pfiles_and_time.list'
my_password = 'XXXXXXX'

# DO NOT CALL NOW
#get_list = "ssh -q cyrf@nl-bpo-dev.ent.dfo-mpo.ca 'find /data/seabird/* -type f -name *.p[1-2][0-9][0-9][0-9] -exec ls -l {} \;' > pfiles_and_time.list"
#get_list_command = "ssh -q " + host + " 'find " + path + "/* -type f -name *.p[1-2][0-9][0-9][0-9] -exec ls -l {} \;' > " + list_files
#os.system(get_list_command)

d = []
with open(list_files) as fp:  
   for idx, line in enumerate(fp):
        line = ' '.join(line.split())
        print idx, line
        n = line.strip().split(' ')[5:8]

        if ':' in n[2]:
            timestamp = pd.to_datetime(n[1]+'-'+n[0]+'-2017 ' +n[2], format='%d-%b-%Y %H:%M') # <--------- this is weak!!!
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

## ---> HERE YOU CAN JUST READ THE PICKLED FILE <--- ##
#df = pd.read_pickle('file_names_unique.pkl')
for file in df['full_path'].values:
        command = 'sshpass -p ' + my_password + ' scp ' + host + ':' + file + ' ' + dest_folder 
        print command
        os.system(command)

        
