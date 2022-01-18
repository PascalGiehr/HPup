'''
Created on Aug 27, 2015

@author: karln
'''

import subprocess
import os
import glob
import re
from lib.functions.aux_generics import collect_files_from_tree, syscall_in_out

def exists(output_file):
    #print("output_exists:")
    #print(output_file)
    assert os.path.exists(output_file), 'Given entity {} does not exist'.format(output_file)
    return output_file

def dummyFunc(inputFile,outputFile,cmd):
    print(inputFile)
    print(outputFile)
    print(cmd)
    
    return outputFile


def createWorkFolder(tmpdir,name):
    cmd= 'mktmp --tmpdir {tmpdir} -d {name}.XXXXXXXX'
    formatter= {'tmpdir': tmpdir, 'name':name}
    cmd= cmd.format(**formatter)
    wd= subprocess.getoutput(cmd)
    
    assert os.path.exists(wd)
    return wd



def getReadsToWd(input_file,
                  output_file,
                  outFolder,
                  fastqPattern,
                  syscall):
    
#     print(input_file)
#     print(output_file)
#     print(outFolder)
#     print(fastqPattern)
#     print(syscall)
    
    inputFiles=collect_files_from_tree(input_file[0], fastqPattern)
    pairs={}
    pairRegExp=re.compile(r"(.*)/([^/]+_)R[12](_[0-9]{3}).fastq.gz")
    
    print("start")
    
    #match paired files
    for file in inputFiles:
        match=pairRegExp.match(file)
        key=' '.join(match.groups())
        if not key in pairs:
            pairs[key]=[]
        pairs[key].append(match.group())
    
    fileList=[]
    i=1
#     outFolder=os.path.dirname(output_file[0])
    cmd='cp -L {inputfile} {outputfile}'
    
    #copy and rename files
    #TODO: return a list with the transformation
    for filePair in pairs.values():
        filePair.sort()
        oFiles=(outFolder + '/seq_R1_%04i.fastq.gz'%i)
        syscall_in_out(filePair[0], oFiles, cmd, syscall)
        if len(filePair)>1:
            oFiles= (oFiles, outFolder + '/seq_R2_%04i.fastq.gz'%i)
            syscall_in_out(filePair[1], oFiles[1], cmd, syscall)
        
        fileList.append(oFiles)
        i+=1
    
    return fileList
    


def locateFiles(input_file, output_file):
    return collect_files_from_tree(input_file[0], output_file)
    #return output_file
    
def trimReads(input_file, output_file):
    print(input_file,output_file)
    return output_file
        
 
    
