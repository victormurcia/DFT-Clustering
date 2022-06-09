#Code Written by Victor M. Murcia @ Washington State University
#5/27/2022
#To use this function adjust the path in line 18 for the string 'infile' to match the directory containing the output files of the 
#transition potential calculations.
#Then in line 74, change the two inputs to correspond to the MO values corresponding to the LUMO (input 1) and the final MO calculated by Stobe (input 2)
#Then in line 76, change the array of values to correspond to the excitation centers used for the calculations.
#Then in line 77, change the value to equal the number of atoms in your molecule. 
#Finally, run the command from a terminal: python moFinder.py
#This will extract the Mulliken Population Analysis for the desired MOs from StoBe for the desired atoms.
# -*- coding: utf-8 -*-
from itertools import islice
import glob, os
#This function reads the output files from StoBe and finds the results of the Mulliken Population analysis
#for a given molecular orbital. It then makes text files for each 
def findMOinAtom(atomNumber,MO,nAtoms):
    desktop = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')
  #  print(desktop)
    infile = 'E:/CuPc/CuPc DFT/CuPc_DFT_Fixed_5-14-2020/CuPc_2005_Isolated/Large Basis Sets/UHF/TP/c'+str(atomNumber)+'tp.out'
    outfile = desktop+'/MO'+str(MO)+'_c'+str(atomNumber) + '.txt'
    search_string = ' Orbital '+str(MO)+' (A   ): Spin alpha'
  #  print(infile)
 #   print(outfile)
    with open(infile,'r') as f:
        for line in f:
            if search_string in line:
               #for line in range(nAtoms+1):
     #           print(line)
                pop = ''.join(islice(f, nAtoms+2))
    #           print(''.join(islice(f, nAtoms+2)))
       #         print(pop)
                f2 = open(outfile, "w")
                f2.write(line + pop)
    f.close()
    f2.close()

#Wrapper function to process multiple MOs and atoms    
def makeMOFiles(moList,atomList,atomsInMolecule):
    
    for mo in moList:
        makeMOfolder(mo)
        for atom in atomList:
            findMOinAtom(atom,mo,atomsInMolecule)
            
def listMOfiles():
    desktop = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')
    os.chdir(desktop)
    #for file in sorted(glob.glob("MO*.txt"), key = last_6chars):
    #    print(file)
    sorted(glob.glob("MO*.txt"), key = lambda x:x[6])

def last_6chars(x):
    return(x[5:])

def makeMOfolder(mo):
    desktop = os.path.join(os.path.join(os.environ['USERPROFILE']), 'Desktop')
    path = desktop+'/'+str(mo)   
    isExist = os.path.exists(path)
   
    if not isExist:
        # Create a new directory because it does not exist 
        os.makedirs(path)
  #      print("The new directory is created!")

def MOvalueGenerator(iniMO,finMO):
    currentMOList = []
    mo = iniMO
    while mo <= finMO:
        currentMOList.append(mo)
        mo+=1

   # print(currentMOList)
    return currentMOList
        
moList = MOvalueGenerator(117,934)#934) #[117,125,255]
print(moList)
atomList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
atomsInMolecule = 57
makeMOFiles(moList,atomList,atomsInMolecule)    
#listMOfiles()     
