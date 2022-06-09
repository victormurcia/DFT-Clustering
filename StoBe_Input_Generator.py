#Code written by Victor M. Murcia Ruiz at Washington State University
import os, shutil
from os import path
import glob
from natsort import natsort_keygen, ns
import molConfig

DEFAULT_SYM = "C1"
DEFAULT_multiplicity = '1'
DEFAULT_GEOUNITS = 'ANGSTROMS'
DEFAULT_RUNTYPE = 'startup nooptimize'
DEFAULT_SCFTYPE = 'direct'
DEFAULT_POTENTIAL = 'nonlocal rpbe pbe'
DEFAULT_GRID = 'fine'
DEFAULT_CHARGE = '0'
DEFAULT_MAXCYCLES = '300'
DEFAULT_ECONVERGENCE = '0.001'
DEFAULT_DCONVERGENCE = '0.001'
DEFAULT_DMIXING = 'mdens 0.05'
DEFAULT_DIIS = 'new 7'
DEFAULT_ORBI = '5d'
DEFAULT_MULLIKEN = 'on full'
DEFAULT_VIRT = 'all'
DEFAULT_FSYMGND = 'scfocc'
DEFAULT_MOFILE = 'molden'
DEFAULT_FSYMEXC = 'scfocc excited'
DEFAULT_ALPHAOCC = '0 1 1 0.0'
DEFAULT_BETAOCC = '0 0'
DEFAULT_ALPHAOCCTP = '0 1 1 0.5'
DEFAULT_SPIN = "FULL"
DEFAULT_XRAY = "XAS"

def makeRunFile(fname, mname, aname, nAtoms, alphaGND, betaGND,element1a,element1b,element2,element3,element4,nElem1,nElem2,nElem3,nElem4,nElem5,nElem6,difElements,
               elem1_Abasis_a,elem1_Abasis_b,elem2_Abasis,elem3_Abasis,elem4_Abasis, elem5_Abasis, elem6_Abasis,
               elem1_Obasis_a, elem1_Obasis_b, elem2_Obasis, elem3_Obasis, elem4_Obasis,elem5_Obasis, elem6_Obasis,
               elem1_MCPbasis,
               title,
               sym, geom,runtype, scftype, potential, grid, multiplicity, charge, maxcycles, econvergence, dconvergence,
               dmixing, diis, orbi, mulliken, virt, spin, mofile, fsymGND,fsymEXC, fsymTP, alfaOcc, alfaOccTP, betaOcc, ftype):

    if sym is None:
        sym = DEFAULT_SYM

    if geom is None:
        geom = DEFAULT_GEOUNITS

    if runtype is None:
        runtype = DEFAULT_RUNTYPE

    if scftype is None:
        scftype = DEFAULT_SCFTYPE

    if potential is None:
        potential = DEFAULT_POTENTIAL

    if grid is None:
        grid = DEFAULT_GRID

    if multiplicity is None:
        multiplicity = DEFAULT_multiplicity

    if charge is None:
        charge = DEFAULT_CHARGE

    if maxcycles is None:
        maxcycles = DEFAULT_MAXCYCLES

    if econvergence is None:
        econvergence =DEFAULT_ECONVERGENCE

    if dconvergence is None:
        dconvergence = DEFAULT_DCONVERGENCE

    if dmixing is None:
        dmixing = DEFAULT_DMIXING

    if diis is None:
        diis = DEFAULT_DIIS

    if orbi is None:
        orbi = DEFAULT_ORBI

    if mulliken is None:
        mulliken = DEFAULT_MULLIKEN

    if virt is None:
        virt = DEFAULT_VIRT

    if fsymGND is None:
        fsymGND = DEFAULT_FSYMGND

    if mofile is None:
        mofile = DEFAULT_MOFILE
    
    if alfaOcc is None:
        alfaOcc = DEFAULT_ALPHAOCC

    if betaOcc is None:
        betaOcc = DEFAULT_BETAOCC
        
    if fsymEXC is None:
        fsymEXC = DEFAULT_FSYMEXC
    
    if fsymTP is None:
        fsymTP = DEFAULT_FSYMEXC
    
    if spin is None:
        spin = DEFAULT_SPIN
        
    i = 1
    #This loop is in charge of creating each .run file
    for i in range(1, nAtoms + 1):
        
        f = open(str(aname) + str(i) + ftype + ".run", "w+",newline="\n")
            
        f.write("#!/bin/csh -f\n")
        f.write("ln -s ~/STOBE/Basis/baslib.new7 fort.3\n")
        f.write("ln -s ~/STOBE/Basis/symbasis.new fort.4\n")
        
        f.write("cat >" + aname + str(i) + ftype + ".inp<</.\n")
            
        f.write("TITLE\n")
        f.write(title + " " + ftype.upper() + "\n")
        
        f.write("SYMMETRY " + sym + "\n")
        f.write("CARTESIAN " + geom + "\n")
        #Enter atomic coordinates from xyz file, effective nuclear charge and integration grid
        n1 = 1
        n2 = 1
        n3 = 1
        n4 = 1
        n5 = 1
        n6 = 1
        #Open the .xyz file to extract the atom labels and atom coordinates
        with open(fname + ".xyz") as f:
            #Open the current .run file using UNIX (LF) encoding to append the contents of the .xyz file 
            with open(aname + str(i) + ftype + ".run", "a+",newline="\n") as f1:
                 #This loop will append two columns to each line form the .xyz file. The first column represents
                 #the effective nuclear charge of the element and the next column represents the number of radial grid points
                 #used to evaluate the exchange/correlation potentials of the element
                for line in f.readlines():
                    #This will build the geometry input according to the number of different elements present in the system
                    if difElements == 6:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32) + '\n'
                            else:
                                line = line.strip() + '     ' + str(element1a) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n1 += 1
                        elif n2 <= nElem2:
                            line = line.strip() + '     ' + str(element2) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n2 += 1
                        elif n3 <= nElem3:
                            line = line.strip() + '     ' + str(element3) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n3 += 1
                        elif n4 <= nElem4:
                            line = line.strip() + '     ' + str(element4) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n4 += 1    
                        elif n5 <= nElem5:
                            line = line.strip() + '     ' + str(element5) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n5 += 1
                        elif n6 <= nElem6:
                            if n6 < nElem6:
                                line = line.strip() + '     ' + str(element6) + '     ' + str(32) + '\n'
                            elif n6 == nElem6:
                                line = line.strip() + '     ' + str(element6) + '     ' + str(32)
                            f1.write(line)
                            n6 += 1                       
                                
                    elif difElements == 5:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32) + '\n'
                            else:
                                line = line.strip() + '     ' + str(element1a) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n1 += 1
                        elif n2 <= nElem2:
                            line = line.strip() + '     ' + str(element2) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n2 += 1
                        elif n3 <= nElem3:
                            line = line.strip() + '     ' + str(element3) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n3 += 1
                        elif n4 <= nElem4:
                            line = line.strip() + '     ' + str(element4) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n4 += 1    
                        elif n5 <= nElem5:
                            if n5 < nElem5:
                                line = line.strip() + '     ' + str(element5) + '     ' + str(32) + '\n'
                            elif n5 == nElem5:
                                line = line.strip() + '     ' + str(element5) + '     ' + str(32)
                            f1.write(line)
                            n5 += 1

                    elif difElements == 4:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32) + '\n'
                            else:
                                line = line.strip() + '     ' + str(element1a) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n1 += 1
                        elif n2 <= nElem2:
                            line = line.strip() + '     ' + str(element2) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n2 += 1
                        elif n3 <= nElem3:
                            line = line.strip() + '     ' + str(element3) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n3 += 1
                        elif n4 <= nElem4:
                            if n4 < nElem4:
                                line = line.strip() + '     ' + str(element4) + '     ' + str(32) + '\n'
                            elif n4 == nElem4:
                                line = line.strip() + '     ' + str(element4) + '     ' + str(32)
                            f1.write(line)
                            n4 += 1
                                
                    elif difElements == 3:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32) + '\n'
                            else:
                                line = line.strip() + '     ' + str(element1a) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n1 += 1
                        elif n2 <= nElem2:
                            line = line.strip() + '     ' + str(element2) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n2 += 1
                        elif n3 <= nElem3:
                            if n3 < nElem3:
                                line = line.strip() + '     ' + str(element3) + '     ' + str(32) + '\n'
                            elif n3 == nElem3:
                                line = line.strip() + '     ' + str(element3) + '     ' + str(32)
                            f1.write(line)
                            n3 += 1
                        
                    elif difElements == 2:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32) + '\n'
                            else:
                                line = line.strip() + '     ' + str(element1a) + '     ' + str(32) + '\n'
                            f1.write(line)
                            n1 += 1
                        elif n2 <= nElem2:
                            if n2 < nElem2:
                                line = line.strip() + '     ' + str(element2) + '     ' + str(32) + '\n'
                            elif n2 == nElem2: 
                                line = line.strip() + '     ' + str(element2) + '     ' + str(32)
                            f1.write(line)
                            n2 += 1
                     
                    elif difElements == 1:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32) + '\n'
                            else:
                                line = line.strip() + '     ' + str(element1a) + '     ' + str(32) + '\n'
                                
                            if n1 == nElem1:
                                line = line.strip() + '     ' + str(element1b) + '     ' + str(32)
                                
                            f1.write(line)
                            n1 += 1    
        
        #The .xyz file will not be read anymore past this point
        f = open(str(aname) + str(i) + ftype + ".run", "a",newline="\n")
        
        f.write("\n")
        f.write("END\n")
        
        #Enter calculation parameters
        f.write("RUNTYPE " + runtype + "\n")
        f.write("SCFTYPE " + scftype + "\n")
        f.write("POTENTIAL " + potential + "\n")
        f.write("GRID " + grid + "\n")
        f.write("MULTIPLICITY " + multiplicity + "\n")
        f.write("CHARGE " + charge + "\n")
        f.write("MAXCYCLES " + maxcycles + "\n")
        f.write("ECONVERGENCE " + econvergence + "\n")
        f.write("DCONVERGENCE " + dconvergence + "\n")
        f.write("DMIXING " + dmixing + "\n")
        f.write("DIIS " + diis + "\n")
        f.write("ORBI " + orbi + "\n")
        f.write("MULLIKEN " + mulliken + "\n")
        f.write("VIRT " + virt + "\n")
        f.write("SPIN " + spin + "\n")
        
        #Depending on what type of .run file we are making choose one of the following ways to describe the electronic state input
        if ftype == "gnd":
            f.write("FSYM " + fsymGND + "\n")
            f.write("ALFA " + str(alphaGND) + "\n")
            f.write("BETA " + str(betaGND) + "\n")
            f.write("FILE " + mofile + "\n")
            f.write("END\n")
            
        elif ftype == "exc":
            f.write("FSYM " + fsymEXC + "\n")
            f.write("ALFA " + str(alphaGND + 1) + "\n")
            f.write("BETA " + str(betaGND) + "\n")
            f.write("SYM 1\n")
            f.write("ALFA " + alfaOcc + "\n")
            f.write("BETA " + betaOcc + "\n")
            f.write("END\n")
            f.write("MULLIKEN " + mulliken + "\n")
            f.write("FILE " + mofile + "\n")
            f.write("END\n")  
            
        elif ftype == "tp":
            f.write("FSYM " + fsymTP + "\n")
            f.write("ALFA " + str(alphaGND) + "\n")
            f.write("BETA " + str(betaGND) + "\n")
            f.write("SYM 1\n")
            f.write("ALFA " + alfaOccTP + "\n")
            f.write("BETA " + betaOcc + "\n")
            f.write("END\n")
            f.write("MULLIKEN on full\n")
            f.write("FILE " + mofile + "\n")
            f.write("XRAY xas\n")
            f.write("END\n")
            f.write("END\n")        
        

        #Enter auxiliary basis set input
        n1 = 1
        n2 = 1
        n3 = 1
        n4 = 1
        n5 = 1
        n6 = 1
        
        with open(str(aname) + str(i) + ftype + ".run",'r') as f:
            with open(str(aname) + str(i) + ftype + ".run", "a",newline="\n") as f1:
                
                for line in f.readlines():
                    if difElements == 6:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Abasis_b + '\n'
                            else:
                                line = elem1_Abasis_a + '\n'
                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Abasis + '\n'
                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Abasis + '\n'
                            f1.write(line)
                            n3 += 1

                        elif n4 <= nElem4:
                            line = elem4_Abasis + '\n'
                            f1.write(line)
                            n4 += 1

                        elif n5 <= nElem5:
                            line = elem5_Abasis + '\n'
                            f1.write(line)
                            n5 += 1

                        elif n6 <= nElem6:
                            line = elem6_Abasis + '\n'
                            f1.write(line)
                            n6 += 1

                    if difElements == 5:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Abasis_b + '\n'
                            else:
                                line = elem1_Abasis_a + '\n'
                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Abasis + '\n'
                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Abasis + '\n'
                            f1.write(line)
                            n3 += 1

                        elif n4 <= nElem4:
                            line = elem4_Abasis + '\n'
                            f1.write(line)
                            n4 += 1

                        elif n5 <= nElem5:
                            line = elem5_Abasis + '\n'
                            f1.write(line)
                            n5 += 1
                            
                    if difElements == 4:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Abasis_b + '\n'
                            else:
                                line = elem1_Abasis_a + '\n'
                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Abasis + '\n'
                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Abasis + '\n'
                            f1.write(line)
                            n3 += 1

                        elif n4 <= nElem4:
                            line = elem4_Abasis + '\n'
                            f1.write(line)
                            n4 += 1

                    if difElements == 3:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Abasis_b + '\n'
                            else:
                                line = elem1_Abasis_a + '\n'
                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Abasis + '\n'
                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Abasis + '\n'
                            f1.write(line)
                            n3 += 1

                    if difElements == 2:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Abasis_b + '\n'
                            else:
                                line = elem1_Abasis_a + '\n'
                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Abasis + '\n'
                            f1.write(line)
                            n2 += 1

                    if difElements == 1:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Abasis_b + '\n'
                            else:
                                line = elem1_Abasis_a + '\n'
                            f1.write(line)
                            n1 += 1
                            
        # Enter orbital basis set input
        n1 = 1
        n2 = 1
        n3 = 1
        n4 = 1
        n5 = 1
        n6 = 1
        with open(str(aname) + str(i) + ftype + ".run",'r') as f:
            with open(str(aname) + str(i) + ftype + ".run", "a",newline="\n") as f1:
                for line in f.readlines():
                    if difElements == 6:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Obasis_b + '\n'
                            else:
                                line = elem1_Obasis_a + '\n'

                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Obasis + '\n'

                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Obasis + '\n'

                            f1.write(line)
                            n3 += 1

                        elif n4 <= nElem4:
                            line = elem4_Obasis + '\n'

                            f1.write(line)
                            n4 += 1

                        elif n5 <= nElem5:
                            line = elem5_Obasis + '\n'

                            f1.write(line)
                            n5 += 1

                        elif n6 <= nElem6:
                            line = elem6_Obasis + '\n'

                            f1.write(line)
                            n6 += 1
                            
                    if difElements == 5:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Obasis_b + '\n'
                            else:
                                line = elem1_Obasis_a + '\n'

                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Obasis + '\n'

                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Obasis + '\n'

                            f1.write(line)
                            n3 += 1

                        elif n4 <= nElem4:
                            line = elem4_Obasis + '\n'

                            f1.write(line)
                            n4 += 1

                        elif n5 <= nElem5:
                            line = elem5_Obasis + '\n'

                            f1.write(line)
                            n5 += 1
                            
                    if difElements == 4:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Obasis_b + '\n'
                            else:
                                line = elem1_Obasis_a + '\n'

                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Obasis + '\n'

                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Obasis + '\n'

                            f1.write(line)
                            n3 += 1

                        elif n4 <= nElem4:
                            line = elem4_Obasis + '\n'

                            f1.write(line)
                            n4 += 1                                               
                            
                    if difElements == 3:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Obasis_b + '\n'
                            else:
                                line = elem1_Obasis_a + '\n'

                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Obasis + '\n'

                            f1.write(line)
                            n2 += 1

                        elif n3 <= nElem3:
                            line = elem3_Obasis + '\n'

                            f1.write(line)
                            n3 += 1  
                            
                    if difElements == 2:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Obasis_b + '\n'
                            else:
                                line = elem1_Obasis_a + '\n'

                            f1.write(line)
                            n1 += 1

                        elif n2 <= nElem2:
                            line = elem2_Obasis + '\n'

                            f1.write(line)
                            n2 += 1
                            
                    if difElements == 1:
                        if n1 <= nElem1:
                            if n1 != i:
                                line = elem1_Obasis_b + '\n'
                            else:
                                line = elem1_Obasis_a + '\n'

                            f1.write(line)
                            n1 += 1  
                            
        # Enter model core potential basis set input
        n1 = 1
        with open(str(aname) + str(i) + ftype + ".run",'r') as f:
            with open(str(aname) + str(i) + ftype + ".run", "a",newline="\n") as f1:
                for line in f.readlines():
                    if n1 < nElem1:
                        line = elem1_MCPbasis + '\n'

                        f1.write(line)
                        n1 += 1
        # Enter augmentation basis set input
        n1 = 1
        with open(str(aname) + str(i) + ftype + ".run", 'r') as f:
            with open(str(aname) + str(i) + ftype + ".run", "a",newline="\n") as f1:
                for line in f.readlines():
                    if n1 == i:
                        line = "X-FIRST\n"
                        f1.write(line)
                        n1+=1
                        break
                    else:
                        line = "X-DUMMY\n"
                        f1.write(line)
                        n1+=1

        #Designate output files and run conditions
        f = open(str(aname) + str(i) + ftype + ".run", "a",newline="\n")
        f.write("END\n")
        f.write("/.\n")
        f.write("~/STOBE/Source/StoBe.x <" + str(aname) + str(i) + ftype + ".inp>& " + str(aname) + str(i) + ftype + ".out\n")
        f.write("mv Molden.molf " + str(aname) + str(i) + ftype + ".molden\n")
        if ftype == "tp":
            f.write("mv fort.11 " + str(aname) + str(i) + ".xas\n")
        f.write("rm fort.*\n")
        i += 1
    
def makeXASrun(mname,aname,nAtoms,title):

	Erange = ""
	Ebroad = ""
	if aname == "C" or "c":
		Erange = "280 320"
		Ebroad = "0.5 12 288 320"
	elif aname == "N" or "n":
		Erange = "400 450"
		Ebroad = "0.5 12 418 450 "
	elif aname == "O" or "o":
		Erange = "530 580"
		Ebroad = "0.5 12 535 560"
	else:
		print("The energy range and broadening has not been implemented for this element yet.")

	print(Erange)
	i=1
	for i in range(1,nAtoms+1):

		f = open(str(aname)+str(i)+"xas.run","w+",newline="\n")

		f.write("#!/bin/csh -f\n")
		f.write("ln -s ~/STOBE/" + mname + "/" + aname +str(i) + "/" + aname +str(i) + ".xas fort.1\n")
		f.write("cat >" + aname + str(i)+"xas.inp<</.\n")
		f.write("title\n")
		f.write(title + " XAS \n")
		f.write("PRINT\n")
		f.write("RANGE " + Erange +"\n")
		f.write("POINTS 2000\n")
		f.write("WIDTH " + Ebroad +"\n")
		f.write("XRAY xas\n")
		f.write("TOTAL 1\n")
		f.write("END\n")
		f.write("/.\n")
		f.write("~/STOBE/Source/xrayspec.x <" + aname + str(i) + "xas.inp>& " +aname + str(i) +"xas.out\n")
		f.write("cp XrayT001.out " + aname + str(i) + ".out\n")
		f.write("rm fort.*\n")
		i+=1

def makeSEQrun(nAtoms,aname):
    
    i=1
    for i in range(1,nAtoms+1):
        f = open(str(aname)+str(i)+"seq.run","w+",newline="\n")
        
        f.write("chmod +x ./" + str(aname) + str(i) + "gnd.run\n")
        f.write("./" + str(aname) + str(i) + "gnd.run\n")
        f.write("\n")
        
        f.write("chmod +x ./" + str(aname) + str(i) + "exc.run\n")      
        f.write("./" + str(aname) + str(i) + "exc.run\n")
        f.write("\n")
        
        f.write("chmod +x ./" + str(aname) + str(i) + "tp.run\n")
        f.write("./" + str(aname) + str(i) + "tp.run\n")
        f.write("\n")
        
        f.write("chmod +x ./" + str(aname) + str(i) + "xas.run\n")
        f.write("./" + str(aname) + str(i) + "xas.run\n")
        i+=1

def makeFolders(aName,nAtoms,mname):

    i=1
    for i in range(1,nAtoms+1):
        dirName = str(aName) + str(i)
        if not os.path.exists(dirName):
            os.mkdir(dirName)
            print("Directory " + dirName + " has been created.")
        else:
            print("Directory " + dirName + " already exists.")
    
    if not os.path.exists(mname):
        os.mkdir(mname)
        print("Main calculation directory has been made.")
    else:
        print("Main calculation directory already exists.")
    
def listFiles(aname,nAtoms): #This function organizes the .run files for all atoms into folders

    iniPath = os.getcwd()
    runfiles = []
    folders  = []

    filenames= os.listdir (".")
    #find . -type f -print0 | xargs -0 dos2unix
    #os.system("find dosdir -type" + "*.run" + " -exec dos2unix -u {} \; ")
    #This populates the list runfiles with any files that have a .run extension
    for file in glob.glob("*.run"):
        #dos2unix file
        runfiles.append(file)
    
    #This populates the list folders with any folders that are present in the current working directory
    for file in filenames:
        if os.path.isdir(os.path.join(os.path.abspath("."), file)):
            folders.append(file)
    
    #Perform a natural sort of the files and folders (i.e. C1,C2,C3...C10,C11,etc, instead of C1,C10,C11...C2,C3)
    natSortKey = natsort_keygen(key=lambda y: y.lower(), alg=ns.IGNORECASE)
    
    runfiles.sort(key = natSortKey)
    folders.sort(key = natSortKey)
    
    #This loop moves the files to their respective atom folders and deletes the version in the original directory
    i = 1
    nf = 0
    for i in range(0,nAtoms + 1):
        atomDir = aname + str(i)
        for item in runfiles:
            if item.startswith(atomDir):
                if nf >= i*5:
                    break
                else:
                    shutil.copy(path.join(iniPath,item),atomDir)
                    os.remove(item)
                    nf += 1        
       
    print("The files are:")
    print(runfiles,"\n")
    
#    print("The folders are:")
#    print(folders,"\n")       
    
def makeAllRunFiles():

    print("********.run FILE GENERATOR FOR StoBe ***********")
    print("*******Written by Victor M. Murcia Ruiz**********")
    print("**********Washington State University************")
    print("*************************************************")
    print("\n")
    print("Please define the variable values to be used to build the files in a file called molConfig.py")
    print("The molConfig.py file should be in the same directory as this python code.")
    print("Ensure that the .xyz file containing the molecular geometry is in the current working directory!")
    print("Important preliminary considerations for the .xyz file.")
    print("1.  The xyz file provided should be sorted by element.") 
    print("2.  The first group of elements in the file should be the element that we are calculating a NEXAFS for.")
    print("3.  There should be a total of 4 columns in the xyz file corresponding to the atom label, and x, y, z coordinates.")
    print("For example, if calculating the Carbon-edge NEXAFS for Benzene (C6H6), the first 6 atoms in the xyz file should be Carbon atoms labeled as C1,C2,C3,C4,C5,C6")
    print("\n")
    
    fname     = molConfig.fname        #This is the name of the .xyz file. Do not include the .xyz extension as part of the name
    ##The provided .xyz file should only have the atom label and xyz coordinates.
    ##The atoms should be sorted by element. Basis set and numbering should match the sorting of the xyz file.
    mname     = molConfig.mname        #This is the name of the molecule we are going to calculate
    aname     = molConfig.aname        #This is the name of the element for which we are doing calculations
    nFiles    = molConfig.nFiles       #This is the number of atoms for which we are generating input files
    difElems  = molConfig.difElems     #This is the number of different ELEMENTS that are present in the molecule
    alpha     = molConfig.alpha        #This is the electronic occupation of the alpha orbitals
    beta      = molConfig.beta         #This is the electronic occupation of the beta orbitals
    element1a = molConfig.element1a    #This is the effective nuclear charge for the atom of element1  that we are simulating a core transition for
    element1b = molConfig.element1b    #This is the effective nuclear charge for the atoms of element1 that we are NOT simulating a core transition for
    element2  = molConfig.element2     #This is the effective nuclear charge for the atoms of element2
    element3  = molConfig.element3     #This is the effective nuclear charge for the atoms of element3
    element4  = molConfig.element4     #This is the effective nuclear charge for the atoms of element4
    element5  = molConfig.element5     #This is the effective nuclear charge for the atoms of element5
    element6  = molConfig.element6     #This is the effective nuclear charge for the atoms of element6
    nElem1    = molConfig.nElem1       #This is the number of atoms of element1 present in the molecule
    nElem2    = molConfig.nElem2       #This is the number of atoms of element2 present in the molecule
    nElem3    = molConfig.nElem3       #This is the number of atoms of element3 present in the molecule
    nElem4    = molConfig.nElem4       #This is the number of atoms of element4 present in the molecule
    nElem5    = molConfig.nElem5       #This is the number of atoms of element5 present in the molecule
    nElem6    = molConfig.nElem6       #This is the number of atoms of element6 present in the molecule

    #This is the input of the different basis sets for the different elements.
    elem1_Abasis_a = molConfig.elem1_Abasis_a       #This is the auxiliary basis set for the atoms of element1 that we are simulating a core transition for
    elem1_Abasis_b = molConfig.elem1_Abasis_b       #This is the auxiliary basis set for the atoms of element1 that we are NOT simulating a core transition for
    elem2_Abasis   = molConfig.elem2_Abasis         #This is the auxiliary basis set for the atoms of element2
    elem3_Abasis   = molConfig.elem3_Abasis         #This is the auxiliary basis set for the atoms of element3
    elem4_Abasis   = molConfig.elem4_Abasis         #This is the auxiliary basis set for the atoms of element4
    elem5_Abasis   = molConfig.elem5_Abasis         #This is the auxiliary basis set for the atoms of element5
    elem6_Abasis   = molConfig.elem6_Abasis         #This is the auxiliary basis set for the atoms of element6
    elem1_Obasis_a = molConfig.elem1_Obasis_a       #This is the orbital basis set for the atoms of element1 that we are simulating a core transition for
    elem1_Obasis_b = molConfig.elem1_Obasis_b       #This is the orbital basis set for the atoms of element1 that we are NOT simulating a core transition for
    elem2_Obasis   = molConfig.elem2_Obasis         #This is the orbital basis set for the atoms of element2
    elem3_Obasis   = molConfig.elem3_Obasis         #This is the orbital basis set for the atoms of element3
    elem4_Obasis   = molConfig.elem4_Obasis         #This is the orbital basis set for the atoms of element4
    elem5_Obasis   = molConfig.elem5_Obasis         #This is the orbital basis set for the atoms of element5
    elem6_Obasis   = molConfig.elem6_Obasis         #This is the orbital basis set for the atoms of element6
    elem1_MCPbasis = molConfig.elem1_MCPbasis       #This is the model core potential basis set for the atoms of element1

    #These are the parameters that define the electron populations/counts. This is part of how self consistency
    multiplicity = molConfig.multiplicity           #This is the multiplicity of your molecule
    alfaOcc      = molConfig.alfaOcc                #This is the electron occupation of the alpha orbitals for the excited state
    betaOcc      = molConfig.betaOcc                #This is the electron occupation of the beta orbitals for the excited and transition potential state
    alfaOccTP    = molConfig.alfaOccTP              #This is the electron occupation of the alpha orbitals of the transition potential state

    #These are parameters that won't need to be changed up too often from the defaults. The default values are at the beginning of this document.
    sym          = molConfig.sym
    geom         = molConfig.geom
    runtype      = molConfig.runtype
    scftype      = molConfig.scftype
    potential    = molConfig.potential
    grid         = molConfig.grid
    charge       = molConfig.charge
    maxcycles    = molConfig.maxcycles
    econvergence = molConfig.econvergence
    dconvergence = molConfig.dconvergence
    dmixing      = molConfig.dmixing
    diis         = molConfig.diis
    orbi         = molConfig.orbi
    mulliken     = molConfig.mulliken
    virt         = molConfig.virt
    spin         = molConfig.spin
    fsymGND      = molConfig.fsymGND
    mofile       = molConfig.mofile
    fsymEXC      = molConfig.fsymEXC 
    title        = molConfig.title
    print("Making GROUND STATE .run FILES")
    makeRunFile(fname, mname, aname, nFiles, alpha, beta,element1a,element1b,element2,element3,element4,nElem1,nElem2,nElem3,nElem4,nElem5,nElem6,difElems,
               elem1_Abasis_a,elem1_Abasis_b,elem2_Abasis,elem3_Abasis,elem4_Abasis, elem5_Abasis, elem6_Abasis,
               elem1_Obasis_a, elem1_Obasis_b, elem2_Obasis, elem3_Obasis, elem4_Obasis,elem5_Obasis, elem6_Obasis,
               elem1_MCPbasis,
               title,
               sym, geom,runtype, scftype, potential, grid, multiplicity, charge, maxcycles, econvergence, dconvergence,
               dmixing, diis, orbi, mulliken, virt, spin, mofile, fsymGND,fsymEXC, fsymEXC, alfaOcc, alfaOccTP, betaOcc,"gnd")       
    
    print("Making EXCITED STATE .run FILES")
    makeRunFile(fname, mname, aname, nFiles, alpha, beta,element1a,element1b,element2,element3,element4,nElem1,nElem2,nElem3,nElem4,nElem5,nElem6,difElems,
               elem1_Abasis_a,elem1_Abasis_b,elem2_Abasis,elem3_Abasis,elem4_Abasis, elem5_Abasis, elem6_Abasis,
               elem1_Obasis_a, elem1_Obasis_b, elem2_Obasis, elem3_Obasis, elem4_Obasis,elem5_Obasis, elem6_Obasis,
               elem1_MCPbasis,
               title,
               sym, geom,runtype, scftype, potential, grid, multiplicity, charge, maxcycles, econvergence, dconvergence,
               dmixing, diis, orbi, mulliken, virt, spin, mofile, fsymGND,fsymEXC, fsymEXC, alfaOcc, alfaOccTP, betaOcc,"exc")
    
    print("Making TRANSITION POTENTIAL .run FILES")
    makeRunFile(fname, mname, aname, nFiles, alpha, beta,element1a,element1b,element2,element3,element4,nElem1,nElem2,nElem3,nElem4,nElem5,nElem6,difElems,
               elem1_Abasis_a,elem1_Abasis_b,elem2_Abasis,elem3_Abasis,elem4_Abasis, elem5_Abasis, elem6_Abasis,
               elem1_Obasis_a, elem1_Obasis_b, elem2_Obasis, elem3_Obasis, elem4_Obasis,elem5_Obasis, elem6_Obasis,
               elem1_MCPbasis,
               title,
               sym, geom,runtype, scftype, potential, grid, multiplicity, charge, maxcycles, econvergence, dconvergence,
               dmixing, diis, orbi, mulliken, virt, spin, mofile, fsymGND,fsymEXC, fsymEXC, alfaOcc, alfaOccTP, betaOcc,"tp")

    print("Making XAS .run FILES")
    makeXASrun(mname,aname,nFiles,title)
    
    print("Making BATCH .run FILES")
    makeSEQrun(nFiles,aname)
    
    print("MAKING DIRECTORIES FOR .run FILES")
    makeFolders(aname,nFiles,mname)
    
    print("ORGANIZING .run FILES")
    listFiles(aname,nFiles)
    
    print("ALL OPERATIONS COMPLETED")
makeAllRunFiles()              