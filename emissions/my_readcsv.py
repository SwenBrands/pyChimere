#!/bin/env python

# import numpy to get this to work!

def my_readcsv(ifile_path,deltype):
    ifile = open(ifile_path, "rb")
    reader = csv.reader(ifile, delimiter=deltype)
    #create a list containing the csv's rows 
    ifile_content=[]
    for row in reader:
        #print row
        ifile_content.append(row)
    ifile.close()
    del(reader,row)
    return(ifile_content)
