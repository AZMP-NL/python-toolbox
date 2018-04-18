import re
import pandas as pd

eoh = '-- DATA --' # end-of-header key
in_header = True

# Read header
header = []
data = []
with open('39170001.p2016', 'r') as td:
    for idx, line in enumerate(td, 1):
        
        if re.match(eoh, line): # end-of-header            
            in_header = False
            continue # ignore line
        elif in_header: # read header
            header.append(line)
        else: # read data
            line = re.sub('\n','', line)
            line = line.strip()
            line = re.sub(' +',' ', line)
            data.append(line.split(' '))

# Read columns & remove last line of header
columns = header[-1]
header = header[:-1]           

# clean columns
columns = re.sub('\n','', columns)
columns = columns.strip()
columns = re.sub(' +',' ', columns)

# to dataFrame       
df = pd.DataFrame(data, columns=columns.split(' '))




        
keyboard
# HERE!!!!!!!
