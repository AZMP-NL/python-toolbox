import os
import pandas as pd

host = 'cyrf@nl-bpo-dev.ent.dfo-mpo.ca' 
path = '/data/seabird'
dest_folder = '/home/cyrf0006/data/dev_database'
list_files = 'pfiles_and_time.list'
my_password = 'orwe123!'

## ---> HERE YOU CAN JUST READ THE PICKLED FILE <--- ##
df = pd.read_pickle('file_names_unique_subset.pkl')
for file in df['full_path'].values:
        command = 'sshpass -p ' + my_password + ' scp ' + host + ':' + file + ' ' + dest_folder 
        print command
        os.system(command)

        
