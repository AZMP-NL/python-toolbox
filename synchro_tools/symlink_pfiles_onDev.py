import os

dest_folder = './symlinks_pfiles/'
list_ascii = 'paths.csv'

# Created this way:
#df = pd.read_pickle('file_names_unique.pkl')
#df['full_path'].to_csv('paths2018.csv', index=False, header=False, sep='\n')
# To copy: scp paths2018.csv cyrf@nl-bpo-dev.ent.dfo-mpo.ca:.

f = open(list_ascii, 'r')
for path in f:
    file = path.strip().split('/')[-1]
    command = 'ln -s ' + path.rstrip('\n') + ' ' + dest_folder + file
    print command
    os.system(command)
