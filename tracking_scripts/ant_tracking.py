#
# Script for loading one ant tracking file
# Loads a single csv from main_experiment/original_data/tracking/
# and puts it into a "flat" pandas DataFrame.
# 
# It's expected to use df.groupby(...) to later analyze 
# a single ant (with df.groupby('id') ) or by timestamp 
# (with df.groupby('t') )
#
# Resulting dataframe has a large number of rows with
# five columns: ['t', 'id', 'x', 'y', 'theta']
# corresponding to integer timestamp, ant tag ID, 
# x and y coordinates, and ant orientation theta (IN DEGREES).
#
# Author: Manuchehr Aminian
#

import csv
import pandas
import numpy as np

#
FOLDER = "C:/Users/maminian/Desktop/ANTS/main_experiment/original_data/tracking/"
FILENAME = "colony020_pathogen_PreTreatment.csv"

# pre-analysis just to figure out how many ants at how many time points we have.
nrows = 0
with open(FILENAME, 'r') as f:
    csvr = csv.reader(f)
    for line in csvr:
        nants = int( line[3] )
        nrows += nants
#

shft = 4    # column number for start of ant records

ants = np.zeros( (nrows,5), dtype=int )
#ants = []

rownum = 0
with open('colony020_pathogen_PreTreatment.csv', 'r') as f:
    csvr = csv.reader(f)
    for j,line in enumerate(csvr):
#        unix = float( line[0] )
        tidx = int( line[1] )
        na = int( line[3] )
        
        for k in range(na):  # loop over all ants; get information.
            aidx = int( line[shft + 4*k + 0] )
            x = int( line[shft + 4*k + 1] )
            y = int( line[shft + 4*k + 2] )
            theta = int(float( line[shft + 4*k + 3] ))
            
#            ants.append( [tidx,aidx,x,y,theta] )
            ants[rownum] = [tidx,aidx,x,y,theta]

            rownum += 1            
        #
#

df = pandas.DataFrame(data=ants, columns=['t','id','x','y','theta'])
