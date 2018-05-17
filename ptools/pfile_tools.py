"""Some tools to process NAFC custom pfiles

References
----------

Atlantic Zone Monitoring Program @NAFC:
https://azmp-nl.github.io/

"""

__author__ = 'Frederic.Cyr@dfo-mpo.gc.ca'
__version__ = '0.1'

## to do:
# - option for dataFrame with z as columns
# - export multiple dataframe into a netCDF or H5
# - From header, xtract lat, lon, station, line, trip, ship and other attributes (meteo?)

import re
import pandas as pd
import pfiles_basics
import numpy as np
import time as tt
import netCDF4 as nc
import os
from sys import version_info


def pfile_variables(filename):
    """Returns the variables (columns) conmtained in a given pfile

    """

    eoh = pfiles_basics.eoh()

    # read line by line until finding eoh
    tmp = []
    with open(filename, 'r') as td:
        for line in td:
        
            if re.match(eoh, line): # end-of-header            
                break
            else:
                tmp.append(line)

    # Read columns & remove last line of header
    columns = tmp[-1]

    # clean columns
    columns = re.sub('\n','', columns)
    columns = columns.strip()
    columns = re.sub(' +',' ', columns)
    columns = columns.split(' ')

    return columns
        

def pfile_header(filename):
    """Returns the header of a given pfile

    """

    eoh = pfiles_basics.eoh()

    # Read header
    header = []
    with open(filename, 'r') as td:
        for line in td:
        
            if re.match(eoh, line): # end-of-header            
                break
            else:
                header.append(line)

    return header

def pfile_to_dataframe(filename):
    """Reads a pfile given in 'filename' as returns a Pandas.DataFrame with
       columns being the variables

    """
    eoh = pfiles_basics.eoh()
    in_header = True

    # Read header
    header = []
    data = []
    with open(filename, 'r') as td:
        for line in td:
        
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
    df = pd.DataFrame(data, columns=columns.split(' '), dtype=float)

    return df

def bin_pressure_from_dataframe(df, Pbin, var):
    """Extract variable 'var' from dataframe generated by pfile_to_dataframe and digitize it to Pbin

  """
    P = np.array(df['pres'])
    Ibtm = np.argmax(P)
    digitized = np.digitize(P[0:Ibtm], Pbin) 
        
    if var in df.columns: # This if part should be put in a function.
        X = df[var]
        X = np.array([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        X = X.reshape(1, len(X)) # row shape for appending in netCDF
    else:
        X = Pbin*np.nan
        X = X.reshape(1, len(X))

    return X

def pfiles_to_netcdf(infiles, nc_outfile, zbin=1, zmax=1500, zshrink=False): # pfiles_to_pannel
    """Given a list of pfiles stored in 'infiles', create a netCDF files with attributes

    Input params:
        - infiles: list of files (e.g. '2017pfiles.list')
          * '$ ls *.p2017 > 2017pfiles.list' would generate the appropriate file list in Linux
        - outfile: output file (e.g. 'AZMP2017.nc')
        - zbin: vertical averaging in final file (zbin=1 is default)
        - zmax: maximum possible depth of the final product (maximum depth may be smaller if zshrink=True, see below)
        - zshrink: if 'True', vertical dimension 'Z' will  be shrink to the first non-empty value (Default is 'False') 
          *NOTE: To open multiple nc files with xarray, Z dim must be the same!

    """
    
    # Check if outfile already exist
    if os.path.exists(nc_outfile):

        print nc_outfile + ' already exists!' 

        py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2

        response_isnt_good = True
        while response_isnt_good:
            if py3:
                response = input(" -> Do you want to remove it? (yes/no?): ")
            else:
                response = raw_input(" -> Do you want to remove it? (yes/no?): ")

            if response == 'yes':
                os.remove(nc_outfile)
                print ' -> ' + nc_outfile + ' removed!' 
                response_isnt_good = False
            elif response == 'no':
                print ' -> Nothing to be done then (an error will be raised)'
                break
            else:
                print ' -> Please answer "yes" or "no"'

    # Generate the list
    filelist = np.genfromtxt(infiles, dtype=str)
    filelist = np.reshape(filelist, filelist.size) # <--------- check if this doesn't cause error
    
    # Original binned depth vector
    Pbin = np.arange(zbin/2.0, zmax, zbin) #will be shrink after if zshrink=True

    
    ##### ------- Loop and store in lists ------- #####
    # Initialize lists
    cast_info_list = []
    cast_index = []
    Tlist = []
    Slist = []
    Clist = []
    SIGlist = []
    Flist = []
    O2list = []
    PHlist = []
    PARlist = []

    for fname in filelist:

        # check if file's there
        if os.path.isfile(fname) is False:
            print ' ->' + fname + ' not found! [skip]'
            continue
        
        # Else, do normal    
        print fname

        # get header
        header = pfile_header(fname)
        #check header
        if 'NAFC_Y2K_HEADER' not in header[0]:
            error_msg = 'Problem with file: header [skip]'
            print error_msg
            expr = 'echo ' + fname + ' : ' + error_msg +' >> .netcdfgen_log.txt'
            os.system(expr)
            continue 
        
        #get cast info and store the info (inspired from J Holden's pfile_IO.py
        cast_info = header[1]
        cast_info = cast_info.replace(',',' ')
        cast_id = cast_info[0:8]

        # Coordinate check
        cast_lat = np.float(cast_info[10:12]) + np.float(cast_info[13:18])/60.0
        cast_lon = np.sign(np.float(cast_info[19:23])) * (np.abs(np.float(cast_info[19:23])) + np.float(cast_info[24:29])/60.0)
        if ((np.int(cast_lon)==0) & (np.int(cast_lat)==0)):
            error_msg = 'Problem with file: (lat,lon) = (0,0) looks wrong [skip]'
            print error_msg
            expr = 'echo ' + fname + ' : ' + error_msg +' >> .netcdfgen_log.txt'
            os.system(expr)
            continue
        
    
        ## if '24:00' in cast_info: #correct 24:00 time
        ##     cast_info = cast_info.replace('24:00', '00:00')
        ##     expr = 'echo ' + fname + ' : time problem >> .netcdfgen_log.txt'
        ##     os.system(expr)
        ## elif '25:' in cast_info: #correct 25:00 time
        ##     cast_info = cast_info.replace('25:', '23:')
        ##     expr = 'echo ' + fname + ' : time problem >> .netcdfgen_log.txt'
        ##     os.system(expr)

        # time check
        if np.int(cast_info[44:46])>59:
            tmp = list(cast_info)
            tmp[44:46]=['0','0']
            #cast_info = cast_info.replace(cast_info[44:46], '00')
            cast_info = "".join(tmp)
            expr = 'echo ' + fname + ' : time problem >> .netcdfgen_log.txt'
            os.system(expr)
        elif np.int(cast_info[41:43])>23:
            tmp = list(cast_info)
            tmp[41:43]=['0','0']
            cast_info = "".join(tmp)
            #cast_info = cast_info.replace(cast_info[41:43], '00')
            expr = 'echo ' + fname + ' : time problem >> .netcdfgen_log.txt'
            os.system(expr)
        elif ((np.int(cast_info[35:37])>12) | (np.int(cast_info[38:40])>31)):
            error_msg = 'Problem with file: wrong date [skip]'
            print error_msg
            expr = 'echo ' + fname + ' : ' + error_msg +' >> .netcdfgen_log.txt'
            os.system(expr)
            continue
        
        # if tests passed, store the rest    
        cast_time = pd.Timestamp(cast_info[29:40] + ' ' + cast_info[40:46])
        cast_sounder = np.int(cast_info[46:51])
        cast_instid = cast_info[51:57]
        cast_instid = cast_instid.replace(' ','')
        cast_set = cast_info[57:61] # unused
        cast_set = cast_set.replace(' ','')
        cast_insttype = cast_info[62]
        cast_insttype = cast_insttype.replace(' ','')
        cast_comment = cast_info[64:78]
        cast_comment = cast_comment.replace(' ','')

        if 'S' in cast_insttype:
            cast_insttype = "V"
        elif 'XBT' in cast_insttype:
            cast_insttype = "F"
        elif 'CTD' in cast_insttype:
            cast_insttype = "V" 
        
        # To DataFrame (and check if empty)
        df = pfile_to_dataframe(fname)
        if df.empty:
            error_msg = 'Problem with file: empty [skip]'
            print error_msg
            expr = 'echo ' + fname + ' : ' + error_msg +' >> .netcdfgen_log.txt'
            os.system(expr)
            continue    
        
       ## ----- Fill wanted variables one by one ----- ##
        # Pressure
        if 'pres' in df.columns:
            P = np.array(df['pres'])            
            Ibtm = np.argmax(P)    
            digitized = np.digitize(P[0:Ibtm], Pbin) #<- this is awesome!
        elif 'depth' in df.columns:
            P = np.array(df['depth'])
            Ibtm = np.argmax(P)    
            digitized = np.digitize(P[0:Ibtm], Pbin) #<- this is awesome!
        else:
            error_msg = 'Problem with file, no pressure channel found [skip]'
            print error_msg
            expr = 'echo ' + fname + ' : ' + error_msg +' >> .netcdfgen_log.txt'
            os.system(expr)
            continue

        # Meta Data
        cast_info_list.append([cast_id, cast_lat, cast_lon, cast_sounder, cast_insttype, cast_instid, cast_comment])
        cast_index.append(cast_time)

        
        # Temperature
        if 'temp' in df.columns:
            X = np.array(df['temp'])
            Tlist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            Tlist.append(list(Pbin*np.nan))

        # Salinity
        if 'sal' in df.columns:
            X = np.array(df['sal'])
            Slist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            Slist.append(list(Pbin*np.nan))

        # Conductivity
        if 'cond' in df.columns:
            X = np.array(df['cond'])
            Clist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            Clist.append(list(Pbin*np.nan))

        # Sigma-t
        if 'sigt' in df.columns:
            X = np.array(df['sigt'])
            SIGlist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            SIGlist.append(list(Pbin*np.nan))

        # Fluorescence
        if 'flor' in df.columns:
            X = np.array(df['flor'])
            Flist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            Flist.append(list(Pbin*np.nan))

        # Oxygen
        if 'oxy' in df.columns:
            X = np.array(df['oxy'])
            O2list.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            O2list.append(list(Pbin*np.nan))

        # PAR
        if 'par' in df.columns:
            X = np.array(df['par'])
            PARlist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            PARlist.append(list(Pbin*np.nan))

        # PH
        if 'ph' in df.columns:
            X = np.array(df['ph'])
            PHlist.append([X[0:Ibtm][digitized == i].mean() for i in range(0, len(Pbin))])
        else:
            PHlist.append(list(Pbin*np.nan))
    #### -------------------------------------------------------- ####


    # List to array (for easier manipulation in next step)
    Tarray = np.array(Tlist)
    Sarray = np.array(Slist)
    Carray = np.array(Clist)
    SIGarray = np.array(SIGlist)
    Farray = np.array(Flist)
    O2array = np.array(O2list)
    PHarray = np.array(PHlist)
    PARarray = np.array(PARlist)

    # Resize Lists to unused maximum depth encountered (a bit weak the way I do it)
    if zshrink:
        idx_nan = np.squeeze(np.where(np.isnan(np.nanmean(Tarray, axis=0))==True))
        Tarray = Tarray[:,0:idx_nan[0]-1]
        Sarray = Sarray[:,0:idx_nan[0]-1]
        Carray = Carray[:,0:idx_nan[0]-1]
        SIGarray = SIGarray[:,0:idx_nan[0]-1]
        Farray = Farray[:,0:idx_nan[0]-1]
        O2array = O2array[:,0:idx_nan[0]-1]
        PHarray = PHarray[:,0:idx_nan[0]-1]
        PARarray = PARarray[:,0:idx_nan[0]-1]
        Pbin = Pbin[0:idx_nan[0]-1]

    # Dataframe content
    columns = ['cast_id', 'lat', 'lon', 'sounder_depth', 'instrument_type', 'instrument_id', 'cast_comment']
    df_info = pd.DataFrame(cast_info_list, index=cast_index, columns=columns)
    df_temp = pd.DataFrame(Tarray, index=cast_index, columns=Pbin)
    df_sali = pd.DataFrame(Sarray, index=cast_index, columns=Pbin)
    df_cond = pd.DataFrame(Carray, index=cast_index, columns=Pbin)
    df_sigt = pd.DataFrame(SIGarray, index=cast_index, columns=Pbin)
    df_fluo = pd.DataFrame(Farray, index=cast_index, columns=Pbin)
    df_oxyg = pd.DataFrame(O2array, index=cast_index, columns=Pbin)
    df_parr = pd.DataFrame(PARarray, index=cast_index, columns=Pbin)
    df_ph = pd.DataFrame(PHarray, index=cast_index, columns=Pbin)


    #### ------ Building netCDF file (inspired from MEOPAR & NCAR examples) ------ #####

    # File name + global attributes
    nc_out = nc.Dataset(nc_outfile, 'w')

    nc_out.Conventions = 'CF-1.6'
    nc_out.title = 'AZMP test netCDF file'
    nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
    nc_out.source = 'https://github.com/AZMP-NL/AZMP_python_toolbox'
    nc_out.references = 'https://azmp-nl.github.io/'
    nc_out.description = 'A test by Frederic.Cyr@dfo-mpo.gc.ca'
    nc.history = 'Created ' + tt.ctime(tt.time())
    ## nc_out.history = """
    ##     Should extract info like this from header:
    ##     [2013-10-30 13:18] Created netCDF4 zlib=True dataset.
    ##     [2013-10-30 15:22] Set depths between 0 and 4m to 4m and those >428m to 428m.
    ##     [2013-10-31 17:10] Algorithmic smoothing.
    ## """
    nc_out.comment = 'Just a trial at the moment, no distribution!'

    # Create dimensions
    time = nc_out.createDimension('time', None)
    level = nc_out.createDimension('level', len(df_temp.columns))

    # Create coordinate variables
    times = nc_out.createVariable('time', np.float64, ('time',))
    levels = nc_out.createVariable('level', np.int32, ('level',))
    # **** NOTE: Maybe consider using ID instead of time for dimension **** #

    # Create 1D variables
    latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
    longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
    cast_IDs = nc_out.createVariable('trip_ID', str, ('time'), zlib=True)
    comments = nc_out.createVariable('comments', str, ('time'), zlib=True)
    instrument_types = nc_out.createVariable('instrument_type', str, ('time'), zlib=True)
    instrument_IDs = nc_out.createVariable('instrument_ID', str, ('time'), zlib=True)
    sounder_depths = nc_out.createVariable('sounder_depth', np.float32, ('time'), zlib=True)

    # Create 2D variables
    temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    sal = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    cond = nc_out.createVariable('conductivity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    sigt = nc_out.createVariable('sigma-t', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    o2 = nc_out.createVariable('oxygen', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    fluo = nc_out.createVariable('fluorescence', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    par = nc_out.createVariable('irradiance', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    ph = nc_out.createVariable('ph', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    times.units = 'hours since 1900-01-01 00:00:00'
    times.calendar = 'gregorian'
    levels.units = 'dbar'
    levels.standard_name = "pressure"
    #levels.valid_range = np.array((0.0, 5000.0))
    levels.valid_min = 0
    temp.units = 'Celsius'
    temp.long_name = "Water Temperature" # (may be use to label plots)
    temp.standard_name = "sea_water_temperature"
    sal.long_name = "Practical Salinity"
    sal.standard_name = "sea_water_salinity"
    sal.units = "1"
    sal.valid_min = 0
    sigt.standard_name = "sigma_t"
    sigt.long_name = "Sigma-t"
    sigt.units = "Kg m-3"
    o2.long_name = "Dissolved Oxygen Concentration" ;
    o2.standard_name = "oxygen_concentration" ;
    o2.units = "mg L-1" ;
    cond.long_name = "Water Conductivity" ;
    cond.standard_name = "sea_water_conductivity" ;
    cond.units = "S m-1" ;
    fluo.long_name = "Chl-a Fluorescence" ;
    fluo.standard_name = "concentration_of_chlorophyll_in_sea_water" ;
    fluo.units = "mg m-3" ;
    par.long_name = "Irradiance" ;
    par.standard_name = "irradiance" ;
    par.units = "umol photons m-2 s-1" ;
    ph.long_name = "water PH" ;
    ph.standard_name = "PH" ;
    ph.units = "unitless" ;

    # Fill structure
    #print 'temp shape before adding data = ', temp.shape
    latitudes[:] = np.array(df_info.lat)
    longitudes[:] = np.array(df_info.lon)
    cast_IDs[:] = np.array(df_info.cast_id)
    comments[:] = np.array(df_info.cast_comment)
    instrument_types[:] = np.array(df_info.instrument_type)
    instrument_IDs[:] = np.array(df_info.instrument_id)
    sounder_depths[:] = np.array(df_info.sounder_depth)

    temp[:,:] = df_temp.values
    sal[:,:] = df_sali.values
    cond[:,:] = df_cond.values
    sigt[:,:] = df_sigt.values
    o2[:,:] = df_oxyg.values
    fluo[:,:] = df_fluo.values
    par[:,:] = df_parr.values
    ph[:,:] = df_ph.values

    # Fill time
    times[:] = nc.date2num(cast_index, units = times.units, calendar = times.calendar)
    levels[:] = Pbin
    #### -------------------------------------------------------- ####
    print 'Done!'

    nc_out.close()
    return None


def pfiles_to_netcdf_unlimitedz(infiles, nc_outfile, zbin=1): # pfiles_to_pannel
    """Given a list of pfiles stored in 'infiles', create a netCDF files with attributes

    Input params:
        - infiles: list of files (e.g. '2017pfiles.lis'**)
        - outfile: output file (e.g. 'AZMP2017.nc')
        - zbin: vertical averaging in final file (zbin=1 is default)

        ** '$ ls *.p2017 > 2017pfiles.list' would generate the appropriate file list in Linux
  """
    # Check if outfile already exist
    if os.path.exists(nc_outfile):

        print nc_outfile + ' already exists!' 

        py3 = version_info[0] > 2 #creates boolean value for test that Python major version > 2

        response_isnt_good = True
        while response_isnt_good:
            if py3:
                response = input(" -> Do you want to remove it? (yes/no?): ")
            else:
                response = raw_input(" -> Do you want to remove it? (yes/no?): ")

            if response == 'yes':
                os.remove(nc_outfile)
                print ' -> ' + nc_outfile + ' removed!' 
                response_isnt_good = False
            elif response == 'no':
                print ' -> Nothing to be done then (an error will be raised)'
                break
            else:
                print ' -> Please answer "yes" or "no"'

    # Generate the list
    filelist = np.genfromtxt(infiles, dtype=str)

    ##### ------- Initialize NetCDF ------- #####
    # File name + global attributes
    nc_out = nc.Dataset(nc_outfile, 'w')

    nc_out.Conventions = 'CF-1.6'
    nc_out.title = 'AZMP netCDF file'
    nc_out.institution = 'Northwest Atlantic Fisheries Centre, Fisheries and Oceans Canada'
    nc_out.source = 'https://github.com/AZMP-NL/AZMP_python_toolbox'
    #nc_out.references = 'a-certain-repo.ca'
    nc_out.description = 'A test by Frederic.Cyr@dfo-mpo.gc.ca'
    #nc_out.tripID = ''
    creation_date = tt.ctime(tt.time())
    nc.history = 'Created ' + creation_date+ ' from NACF Pfiles'
    ## nc_out.history = """
    ##     Should extract info like this from header:
    ##     [2017-10-15 13:18] Some seabird cleaning
    ##     [2017-10-20 15:22] transfer to pfiles
    ##     [2013-10-22 15:18] creating netCDFNetCDF
    ## """
    nc_out.comment = 'Just a trial at the moment, no distribution!'

    # Create dimensions
    time = nc_out.createDimension('time', None)
    level = nc_out.createDimension('level', None)

    # Create coordinate variables
    times = nc_out.createVariable('time', np.float64, ('time',))
    levels = nc_out.createVariable('level', np.int32, ('level',))

    # Create 1D variables
    latitudes = nc_out.createVariable('latitude', np.float32, ('time'), zlib=True)
    longitudes = nc_out.createVariable('longitude', np.float32, ('time'), zlib=True)
    cast_IDs = nc_out.createVariable('trip_ID', str, ('time'), zlib=True)
    comments = nc_out.createVariable('comments', str, ('time'), zlib=True)
    instrument_types = nc_out.createVariable('instrument_type', str, ('time'), zlib=True)
    instrument_IDs = nc_out.createVariable('instrument_ID', str, ('time'), zlib=True)
    sounder_depths = nc_out.createVariable('sounder_depth', np.float32, ('time'), zlib=True)

    # Create 2D variables
    temp = nc_out.createVariable('temperature', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    sal = nc_out.createVariable('salinity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    cond = nc_out.createVariable('conductivity', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    sigt = nc_out.createVariable('sigma-t', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    o2 = nc_out.createVariable('oxygen', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    fluo = nc_out.createVariable('fluorescence', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)
    par = nc_out.createVariable('irradiance', np.float32, ('time', 'level'), zlib=True, fill_value=-9999)

    # Variable Attributes
    latitudes.units = 'degree_north'
    longitudes.units = 'degree_east'
    times.units = 'hours since 1900-01-01 00:00:00'
    times.calendar = 'gregorian'
    levels.units = 'dbar'
    levels.standard_name = "pressure"
    #levels.valid_range = np.array((0.0, 5000.0))
    levels.valid_min = 0
    temp.units = 'Celsius'
    temp.long_name = "Water Temperature" # (may be use to label plots)
    temp.standard_name = "sea_water_temperature"
    sal.long_name = "Practical Salinity"
    sal.standard_name = "sea_water_salinity"
    sal.units = "1"
    sal.valid_min = 0
    sigt.standard_name = "sigma_t"
    sigt.long_name = "Sigma-t"
    sigt.units = "Kg m-3"
    o2.long_name = "Dissolved Oxygen Concentration" ;
    o2.standard_name = "oxygen_concentration" ;
    o2.units = "mg L-1" ;
    cond.long_name = "Water Conductivity" ;
    cond.standard_name = "sea_water_conductivity" ;
    cond.units = "S m-1" ;
    fluo.long_name = "Chl-a Fluorescence" ;
    fluo.standard_name = "concentration_of_chlorophyll_in_sea_water" ;
    fluo.units = "mg m-3" ;
    par.long_name = "Irradiance" ;
    par.standard_name = "irradiance" ;
    par.units = "umol photons m-2 s-1" ;
    ##### ------------------------------------- #####

    for idx, fname in enumerate(filelist):

        # check if file's there
        if os.path.isfile(fname) is False:
            print ' ->' + fname + ' not found! [skip]'
            continue

        # File present, good to go
        print fname
        # get header & cast info
        header = pfile_header(fname)
        cast_info = header[1]
        cast_info = re.sub('\n','', cast_info)
        cast_info = re.sub(' +',' ', cast_info)
        cast_info = cast_info.split(' ')
        cast_id = np.int(cast_info[0])
        cast_lat = np.float(cast_info[1]) + np.float(cast_info[2])/60.0
        cast_lon = np.sign(np.float(cast_info[3])) * (np.abs(np.float(cast_info[3])) + np.float(cast_info[4])/60.0)
        cast_time = pd.Timestamp(cast_info[5] + ' ' + cast_info[6])
        cast_sounder = np.int(cast_info[7])
        cast_type = cast_info[8]
        # This is weak, better solution to be found!
        if len(cast_info) >= 12:
            cast_station = cast_info[11]
        else:
            cast_station = 'n/a'

        # Fill cast info
        latitudes[idx] = cast_lat
        longitudes[idx] = cast_lon
        cast_IDs[idx] = cast_id
        stations[idx] = cast_station
        instruments[idx] = cast_type
        sounder_depths[idx] =cast_sounder

        # get data
        df = pfile_to_dataframe(fname)

       ## ----- Fill wanted variables one by one ----- ##
        # Find pressure
        if 'pres' in df.columns:
            P = np.array(df['pres'])
            Pbin = np.arange(1, np.max(P), zbin) 
        else:
            print 'Problem with file, no pressure channel found [skip]'
            continue

        # Check if Pbin not empty (in this case likely single surface measurement)
        #  NOTE: Depths at "0m" will be reassigned at 1m (minimum depth of dataset)
        if len(Pbin) == 0:
            print ' -> File with single measurement at ' + np.str(P) + 'm'
            idx = (np.abs(levels[:]-P)).argmin()
            Pbin = np.array([levels[idx]])
            
        # Fill dimension
        times[idx] = nc.date2num(cast_time, units = times.units, calendar = times.calendar)
        if len(Pbin) > len(levels): # update depth if this cast is deeper than present max depth
            levels[:] = np.array(Pbin)

            
        # Fill nc file with variables
        temp[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'temp')
        sal[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'sal')
        cond[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'cond')
        sigt[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'sigt')
        fluo[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'flor')
        o2[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'oxy')
        par[idx,:] = bin_pressure_from_dataframe(df, Pbin, 'par')
    # -------------------------------------------------------- #    
    print 'Done!'
    
    nc_out.close()

    return None
