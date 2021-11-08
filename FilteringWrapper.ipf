#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

constant fwhmconversion = 2.355 	//FWHM = 2*sqrt(2*ln(2))*STD = 2.355*STD

Function filteringRoutines(fnum,iniFolder,baseFolderName,atomName,erange,tval,broad1,broad2,Eini,Efin,gRes,ovpMax,corrWave,dClusters,LUMOwave,IPwave,mol)

	Variable fnum,erange,tval,broad1,broad2,Eini,Efin,gRes,ovpMax,dClusters
	String iniFolder,baseFolderName,atomName,mol
	Wave corrWave
	Wave LUMOwave
	Wave IPWave
	
	//This wrapper function contains routines that do the following:
	//1. Determine the maximum oscillator strength (OS) coming from all atoms
	//2. Remove transitions that are below the OS threshold value
	//3. Remove transitions that are above a given transition energy
	//4. Makes 2D parameter wave for these transitions
	//5. Makes the peaks corresponding to the OS filtered transitions 
	//6. Calculates the overlap matrix for ALL transition pairs [OPTIONAL]
	//7. Organizes waves into folders
	print "Starting operations on original transitions..."
	initialTransitionsWrapper(fnum,atomName,erange,baseFolderName,Eini,Efin,gRes,broad1,broad2,tval,corrWave,LUMOwave,IPwave)
	
	//This wrapper function contains routines that do the following:
	//1. Makes 2D parameter wave for these transitions
	//2. Makes the peaks corresponding to the OS filtered transitions 
	//3. Calculates the overlap matrix for ALL transition pairs [OPTIONAL]
	//4. Organizes waves into folders
	print "Starting operations on OS filtered transitions..."
	osFilteredTransitionsWrapper(fnum,atomName,erange,baseFolderName,Eini,Efin,gRes,broad1,broad2,tval,corrWave,LUMOwave,IPwave)
	
	Variable F1Transitions
	Wave allParamsF1Sorted,ovpWaveAll
	Duplicate/O  ovpWaveAll,ovpOriginal
	
	//Identify the DFT transitions with the final transition clusters
	consolidateTransitionClusters(allParamsF1Sorted)
	
	//Start transition clustering based on symmetry
	print "Starting clustering of transitions using overlap threshold of " + num2str(ovpMax) + " %"
	String ffFolder = GetdataFolder(1)
	Variable cMethod = 2//
	clusteringTransitions(allParamsF1Sorted,ovpWaveAll,tval,ovpMax,Eini,Efin,"ALL",cMethod,mol)
	SetDataFolder $ffFolder
	
End

Function initialTransitionsWrapper(fnum,atomName,erange,baseFolderName,Eini,Efin,gRes,broad1,broad2,tval,corrWave,LUMOwave,IPwave)
	
	Variable fnum
	String atomName
	Variable erange
	String baseFolderName
	Variable Eini
	Variable Efin
	Variable gRes
	Variable broad1
	Variable broad2
	Variable tval
	Wave corrWave
	Wave LUMOwave
	Wave IPWave
	
	//Determine the maximum oscillator strength for each Carbon atom within 280-erange eV
	determineMaxOS(fnum,atomName,erange)
	Wave maxOS
	WaveStats/Q maxOS
	Variable/G globalMaxOS = V_max
	//Set minimum OS threshold for ALL atoms as a given percentage of maximum OS.
	//If tval is set to zero then consider ALL transitions within 280-erange eV
	Variable/G minOSThreshold = (tval/100)*V_max 
	print "The Oscillator Strength threshold value at " + num2str(tval) +"% is " + num2str(minOSThreshold)
	
	print "Making energy filtered waves..."
	//Make filtered OS,Energy, and Theta waves that correspond to threshold constraint 
	energyFilter(fnum,atomName,erange,minOSThreshold,LUMOwave)
	print "Energy filtered waves made"
	
	print "Organizing original data"
	//Organize Data from into Folders
	organizeEnergyFilter(baseFolderName)
	
	print "Making peaks for original transitions"
	//Prepare to make original peaks prior to filtering
	origPeaksPrep(baseFolderName,atomName,fnum,Eini,Efin,gRes,"f0")
	Wave nF0Trans
	WaveStats/Q nF0Trans
	
	print "Making 2D parameter wave for original peaks"
	//Make 2D parameter Waves for Original Peaks
	makePWave(fnum,broad1,broad2,erange,"f0",baseFolderName,corrWave,LUMOwave,IPwave)
	Wave allParams
	Wave allParamsSorted
	//Look for any transitions that have an oscillator strength of 0,NaN, or Inf and remove them from pWave
	removeZeroAmplitude(allParamsSorted)
	Wave allParamsSorted
	
	//Make peaks from each TDM component
	print "Making peaks derived from TDM components"
	makeNEXAFS(allParamsSorted,Eini,Efin,gRes,"f0")
	print "Organizing original peaks"
	//Organize original peaks
	organizePeaks(atomName,baseFolderName,"f0")
	print "Making 2D wave containing original peaks"
	//Make 2D wave containing original peaks
	//makeGPeaks2D(allParamsSorted,nF0Trans,Eini,Efin,gRes,"f0")
End	

Function osFilteredTransitionsWrapper(fnum,atomName,erange,baseFolderName,Eini,Efin,gRes,broad1,broad2,tval,corrWave,LUMOwave,IPwave)

	Variable fnum
	String atomName
	Variable erange
	String baseFolderName
	Variable Eini
	Variable Efin
	Variable gRes
	Variable broad1
	Variable broad2
	Variable tval
	Wave corrWave
	Wave LUMOwave
	Wave IPwave
	
	print "Preparing to make peaks from first filtering..."
	//Prepare to make peaks from first filtering (i.e. filtered by oscillator strength)
	origPeaksPrep(baseFolderName,atomName,fnum,Eini,Efin,gRes,"f1")
	Wave nF1Trans
	WaveStats/Q nF1Trans
	
	print "Making parameter wave for first filter peaks..."
	//Make 2D parameter Waves for First Filter Peaks
	makePWave(fnum,broad1,broad2,erange,"f1",baseFolderName,corrWave,LUMOwave,IPwave)	
	Wave allParamsF1Sorted
	Variable F1Transitions = DimSize(allParamsF1Sorted,0)
		
	//Calculate percent overlap between all filtered peaks using the average OS
	print "Determine overlap for OS filtered peaks"
	peakOverlap(allParamsF1Sorted,"All")

	//Make peaks from each TDM component
	print "Making peaks derived from TDM components"
	makeNEXAFS(allParamsF1Sorted,Eini,Efin,gRes,"f1")
	print "Organizing data..."
	//Organize the peaks from first filtering into a folder
	organizePeaks(atomName,baseFolderName,"f1")	
	
	print "Making 2D  Wave for Oscialltor Strength filtered peaks..."
	//Make 2D wave containing filtered peaks
	//makeGPeaks2D(allParamsF1Sorted,nF1Trans,Eini,Efin,gRes,"f1")
End

Function newThreshold(pWave,osThreshold,tdmSym,Emax)
	Wave pWave
	Variable osThreshold
	String tdmSym
	Variable Emax
	
	osThreshold = osThreshold/100
	String pWaveName ="pWave_" + tdmSym
	Duplicate/O pWave,$pWaveName
	Wave w = $pWaveName
	
	
	Variable minOSThreshold,j,x=DimSize(pWave,0),osCol 
	//Remove transitions that are above energy threshold
	for(j=x-1;j>=0;j-=1)
		Variable en = pWave[j][0]
		if(en > Emax)
			DeletePoints j,1,w
		endif
	endfor
	
	WaveStats/Q/RMD=[][1,1] pWave
	minOSThreshold = osThreshold*V_max
	osCol = 1

	
	x=DimSize(w,0)
		
	for(j=x-1;j>=0;j-=1)
		Variable cOSVal = pWave[j][osCol]
		if(cOSVal == 0)
			DeletePoints j,1,w
		elseif(cOSVal >= minOSThreshold)
			w[j] = pWave[j]	
		elseif(cOSVal < minOSThreshold)
			DeletePoints j,1,w
		endif
	endfor	
End

Function determineMaxOS(fnum,atomName,erange)
	
	Variable fnum,erange
	String atomName
	
	Variable i
	//Determine the maximum oscillator strength for each Carbon atom within 280-erange eV
	Make/O/N=(fnum) maxOS
	print "Upper Energy[eV]	Max OS		Energy of MaxOS[eV]"
	for(i=0;i<=fnum-1;i+=1)
		String cAtom = atomName + num2str(i+1)
		String osWaveName = "OS_" + cAtom
		String cEnName = "eV_" + cAtom
		Wave currOS = $osWaveName 
		Wave eWave = $cEnName
		Variable finalE = BinarySearch($cEnName,erange)	
		WaveStats/Q/R=[0,finalE] $osWaveName
		maxOS[i] = V_max
		Variable OSmaxEnergy = eWave[V_maxloc]
		print "     ",finalE,"   ",V_max,"     ",OSmaxEnergy
	endfor
End

Function energyFilter(fnum,atomName,erange,minOSThreshold,LUMOwave)
	Variable fnum,erange,minOSThreshold
	String atomName
	Wave LUMOwave
	
	Variable i,j,k=0,l
	print "Initial Transitions	OS Filtered Transitions"
	for(i=1;i<=fnum;i+=1)
		String cAtom = atomName + num2str(i)
		String f1OSname = "f1OS_" + cAtom
		String f1enName = "f1En_" + cAtom
		String f1SymName = "f1sym_" + cAtom
		String f1MOsName = "f1MOs_" + cAtom
		
		String osWaveName = "OS_" + cAtom
		String cEnName = "eV_" + cAtom
		String symName = "symmetryParams" +cAtom 
		
		Wave os = $osWaveName
		Wave en = $cEnName
		Wave sym = $symName
		Variable nTrans = numpnts(os)
		
		Duplicate/D/O os,$f1OSname
		Duplicate/D/O en,$f1enName
		Duplicate/D/O sym,$f1symName
		Wave f1OS = $f1OSname
		Wave f1En = $f1enName
		Wave f1sym = $f1SymName
		Make/O/N=(nTrans) $f1MOsName 
		Wave mos = $f1MOsName
		
		Variable nErange = BinarySearch(en,erange)
		DeletePoints nErange,nTrans-nErange,$f1OSname,$f1enName,$f1SymName,$f1MOsName 
		Variable lumoVal = LUMOwave[i-1]
		Variable nTrans2 =numpnts($f1OSname)
		
		for(l=0;l<nTrans2;l+=1)
			mos[l] = lumoVal
			lumoVal+=1
		endfor
		
		//If the OS value is greater than or equal to threshold value, keep the value, else delete it
		for(j=nTrans2-1;j>=0;j-=1)
			if(minOSThreshold > 0)
				if(f1OS[j] >= minOSThreshold)
					f1OS[j] = f1OS[j]	
				elseif(f1OS[j] < minOSThreshold)
					DeletePoints j,1,$f1OSname,$f1enName,$f1symName,$f1MOsName
				endif
			//If the OS threshold is set to 0, remove transitions that have an OS of 0.
			elseif(minOSThreshold == 0)
				if(f1OS[j] > minOSThreshold)
					f1OS[j] = f1OS[j]	
				elseif(f1OS[j] <= minOSThreshold)
					DeletePoints j,1,$f1OSname,$f1enName,$f1symName,$f1MOsName
				endif
			else
				print "Error while filtering oscillator strengths"
			endif
			
			if(k>j)
				break
			endif
		endfor	
		
		Variable nTransf1 = numpnts(f1OS)
	
		print "       ",nTrans,"          ", nTransf1
	endfor
End

Function organizeEnergyFilter(baseFolderName)
	
	String baseFolderName
	
	NewDataFolder/O originalPeaks
	NewDataFolder/O firstFilter
	
	String originalOS = WaveList("!f1*",";","")
	String firstFilterOS = WaveList("f1*",";","")
	
	Variable oOS = ItemsInList(originalOS)
	Variable ffOS = ItemsInList(firstFilterOS)
	
	Variable i
	for(i=0;oOS-1;i+=1)
		String cOriOS = StringFromList(i,originalOS)
		String cOriOSNewName = baseFolderName + ":" + "originalPeaks:" +StringFromList(i,originalOS)
		Wave w1 = $cOriOS
			
		if(WaveExists($cOriOS))
			Duplicate/O w1,$cOriOSNewName
			KillWaves/Z w1
		else 
			break
		endif		
	endfor
	
	for(i=0;ffOS-1;i+=1)
		String cffOS = StringFromList(i,firstFilterOS)
		String cffOSNewName =  baseFolderName + ":" + "firstFilter:" + StringFromList(i,firstFilterOS)
		Wave w2 = $cffOS	
			
		if(WaveExists($cffOS))
			Duplicate/O w2,$cffOSNewName
			KillWaves/Z w2
		else 
			break
		endif
	endfor
End

Function origPeaksPrep(baseFolderName,atomName,fnum,Eini,Efin,gRes,fStep)
	
	String baseFolderName,atomName,fStep
	Variable fnum,Eini,Efin,gRes
	
	String filterStage
	if(StringMatch(fStep,"f0")==1)
		filterStage = ":originalPeaks:"
	elseif(StringMatch(fStep,"f1")==1)
		filterStage = ":firstFilter:"
	endif
	
	//Define folder names
	String cFolder         = baseFolderName + filterStage
	String iniEnergies     = cFolder + "Energies"
	String iniOS           = cFolder + "OS"
	String iniSymmetryWave = cFolder + "symmetryWaves"
	String tdmCompSpecs    = cFolder + "tdmCompSpecs"
	String moFolder        = cFolder + "MOs"
	
	SetDataFolder cFolder
	
	//Make folders	
	NewDataFolder/O $iniEnergies
	NewDataFolder/O $iniOS
	NewDataFolder/O $iniSymmetryWave	
	newDataFolder/O $tdmCompSpecs
	NewDataFolder/O $moFolder
	
	Variable i
	String cAtom, OSname, nTransName,totalSpec
	//Determine how many peaks to make per atom
	for(i=1;i<=fnum;i+=1)
		cAtom = atomName + num2str(i)
		
		if(StringMatch(fStep,"f0")==1)
			OSname = "OS_" + cAtom
			nTransName = "nF0Trans"
			totalSpec = "Total_Specf0"
		elseif(StringMatch(fStep,"f1")==1)
			OSname = "f1OS_" + cAtom
			nTransName = "nF1Trans"
			totalSpec = "Total_Specf1"
		endif 
		
		//This segment determines how many transitions are derived from every atom and puts them into wave "x"
		Wave w = $OSname
		Make/O/N=(fnum) $nTransName
		Wave x = $nTransName
		Variable n = numpnts($OSname)
		x[i-1] = n
		
		//Check if the total NEXAFS made by adding all the atomic NEXAFS exists. If it does, set it to zero, if not, make it.
		if(WaveExists($totalSpec))
			Wave sumSpec = $totalSpec
			sumSpec = 0
		else
			Make/O/N=(gRes) $totalSpec
			SetScale/I x,Eini,Efin,"",$totalSpec
		endif
	endfor
End

Function organizePeaks(atomName,baseFolderName,algoStage)

	String atomName,baseFolderName,algoStage
	
	//This will determine the data folder to use as an anchor point
	String filterStage
	if(StringMatch(algoStage,"f0")==1)
		filterStage = ":originalPeaks:" 
	elseif(StringMatch(algoStage,"f1")==1)
		filterStage = ":firstFilter:" 
	endif
	
	//Definition of data folders
	String cFolder = baseFolderName + filterStage
	String iniEnergies     = cFolder + "Energies:"
	String iniOS           = cFolder + "OS:"
	String iniSymmetryWave = cFolder + "symmetryWaves:"
	String tdmCompSpecs    = cFolder + "tdmCompSpecs:"
	
	SetDataFolder cFolder
		
	String energies,oscillatorStrengths,symmetryWaves,totalSpecs,parameterWaves,moList
	//Decide what waves to look for depending on what stage 
	if(StringMatch(algoStage,"f0")==1)
		energies = WaveList("eV_*",";","")
		oscillatorStrengths = WaveList("OS_*",";","")
		symmetryWaves = WaveList("symmetryParams*",";","")
		totalSpecs = WaveList("Totalf0_*",";","")
		parameterWaves = WaveList("pWave*",";","")
		Wave Total_Specf0
	elseif(StringMatch(algoStage,"f1")==1)
		energies = WaveList("f1En_*",";","")
		oscillatorStrengths = WaveList("f1OS_*",";","")
		symmetryWaves = WaveList("f1sym_*",";","")
		totalSpecs = WaveList("Totalf1_*",";","")
		Wave Total_Specf1
	endif
	
	//Name of total NEXAFS generated from each tensor element
	String sumAVG = tdmCompSpecs + "Total_Spec" + algoStage
	
	//Number of Miscellaneous waves
	Variable nEnergies   = ItemsInList(energies)
	Variable nOSs        = ItemsInList(oscillatorStrengths)
	Variable nSymWaves   = ItemsInList(symmetryWaves)
	Variable nTotalSpecs = ItemsInList(totalSpecs)
	
	Variable i,nAtoms = nOSs
	String cAtom
	//Determine what waves to duplicate from each tensor element		
	if(StringMatch(algoStage,"f0")==1)
		Duplicate/O Total_Specf0,$sumAVG
	elseif(StringMatch(algoStage,"f1")==1)
		Duplicate/O Total_Specf1,$sumAVG
	endif
	
	//This will transfer the miscellaneous waves into the correct folder
	for(i=0;i<=nAtoms-1;i+=1)
		cAtom = atomName + num2str(i+1)
		String cEnergiesNameOrig = StringFromList(i,energies)
		String cEnergiesName     = iniEnergies + StringFromList(i,energies)
		String cOSNameOrig       = StringFromList(i,oscillatorStrengths)
		String cOSName           = iniOS + StringFromList(i,oscillatorStrengths)			
		String cSymWaveNameOrig  = StringFromList(i,symmetryWaves)
		String cSymWaveName      = iniSymmetryWave + StringFromList(i,symmetryWaves)	
		
		Wave en = $cEnergiesNameOrig
		Wave OS = $cOSNameOrig
		Wave symWave = $cSymWaveNameOrig
		
		if(WaveExists(en))
			Duplicate/O en,$cEnergiesName
			KillWaves/Z en
		endif
		
		if(WaveExists(OS))
			Duplicate/O OS,$cOSName
			KillWaves/Z OS
		endif
		
		if(WaveExists(symWave))
			Duplicate/O symWave,$cSymWaveName
			KillWaves/Z symWave
		endif		
	endfor
	
	//Consolidate data folders	
	SetDataFolder cFolder
	
	NewDataFolder/O Spectra
	String pathToSpecs = cFolder + "Spectra:"
	DuplicateDataFolder/O=1/Z $tdmCompSpecs,$pathToSpecs
	
	NewDataFolder/O Misc
	String pathToMisc = cFolder + "Misc:"
	DuplicateDataFolder/O=1/Z $iniEnergies,$pathToMisc
	DuplicateDataFolder/O=1/Z $iniOS,$pathToMisc
	DuplicateDataFolder/O=1/Z $iniSymmetryWave,$pathToMisc
	
	//Kill redundant folders
	KillDataFolder $tdmCompSpecs
	KillDataFolder $iniEnergies
	KillDataFolder $iniOS
	KillDataFolder $iniSymmetryWave
End

Function makeGPeaks2D(pWave,Eini,Efin,gRes,fStep) //Makes peaks for the transitions sorted by ascending energy
	
	Wave pWave
	Variable Eini,Efin,gRes
	String fStep
	
	String name2D
	if(StringMatch(fStep,"f0")==1)	
		name2D = "OriallPeaks2D"
	elseif(StringMatch(fStep,"f1")==1)
		name2D = "f1allPeaks2D"
	endif
	
	Make/O/N=(DimSize(pwave,0),gRes) $name2D
	setscale/i y,Eini,Efin,$name2D
	Wave w = $name2D

	Multithread w = (pWave[p][1])*Gauss(y,pWave[p][0],pWave[p][2])	//Columns 0,1 and 2 from pWave have the peak position,amplitude and width respectively
End

Function makePWave(fnum,broad1,broad2,ewid2,fStep,baseFolderName,corrWave,LUMOwave,IPwave)//Makes a 2D wave containing the transition energies, peak widths, OS and TDM
									//theta's for all atoms	
	Variable fnum,broad1,broad2,ewid2
	String fStep,baseFolderName
	Wave corrWave
	Wave LUMOwave
	Wave IPWave
	
	print "Making parameter wave..."
	
	Variable i,j,n=0
	broad1 = broad1/fwhmConversion
	broad2 = broad2/fwhmConversion
	String filterStep
	
	if(StringMatch(fStep,"f0")==1)
		filterStep = "originalPeaks"
	elseif(StringMatch(fStep,"f1")==1)
		filterStep = "firstFilter"
	endif
	
	String transWaveName = baseFolderName + ":" + filterStep + ":" + "n" + fStep + "Trans"
	
	Wave nFTrans =$transWaveName
	
	for(i=0;i<=fnum-1;i+=1)
		n += nFTrans[i]
	endfor
	
	String pWaveName
	
	if(StringMatch(fStep,"f0")==1)
		pWaveName = "allParams"
	elseif(StringMatch(fStep,"f1")==1)
		pWaveName = "allParamsF1"
	endif
	
	Wave allParams = $pWaveName
	
	Make/O/N=(n,17) allParams
	//Col0 is dft transition energy
	//col1 is the total oscillator strength,
	//col2 is the transition width, 
	//col3, col4, and col5 are the x,y, and z components of the transition dipole moment respectively
	//col6 is the angle between the TDM and the molecular z-axis
	//Col7 is the transition symmetry
	//Col8,9,10,11,12,13 are the xx,yy,zz,xy,xz,and yz tensor elements respectively
	//Col14 is the atom from which that transition derives from
	//Col15 is the molecular orbital associated with that transition
	//Col16 is the cluster that transition belongs to.
	for(i=1;i<=fnum;i+=1)
		if(StringMatch(fStep,"f0")==1)
			String energyList = SortList(WaveList("eV_*",";",""),";",16)
			String osList     = SortList(WaveList("OS_*",";",""),";",16)
			String symList    = SortList(WaveList("symmetry*",";",""),";",16)
		elseif(StringMatch(fStep,"f1")==1)
			energyList    = SortList(WaveList(fStep+"En_*",";",""),";",16)
			osList        = SortList(WaveList(fStep+"OS_*",";",""),";",16)
			symList       = SortList(WaveList(fStep+"sym_*",";",""),";",16)
			String moList = SortList(WaveList(fStep + "MOs*",";",""),";",16)
		endif
		
		Variable k=0 
		
		String currEnWave = StringFromList(i-1,energyList)
		String currOSWave = StringFromList(i-1,osList)
		String currSymWave = StringFromList(i-1,symList)
		if((StringMatch(fStep,"f1")==1))
			String currLUMOWave = StringFromList(i-1,moList)
			Wave z = $currLUMOWave
		endif
		
		Wave w = $currEnWave
		Wave x = $currOSWave
		Wave y = $currSymWave
		Variable tot
		Variable LUMOpos = LUMOwave[i-1]
		if(i==1)
			j=0
			tot = NFtrans[i-1]
		elseif(i>1)
			j = tot
			tot += NFtrans[i-1]
		endif
		
		//////Use IP wave////
		for(j=j;j<=tot-1;j+=1)
			allParams[j][0] = w[k]	//All Energies
			allParams[j][1] = x[k]	//All OS
			
			////5/8/2021
			Variable ewid1 = IPwave[i-1]	///Use IP to define width energy 1
			////////////
			if(w[k] <= ewid1)		//All widths
				allParams[j][2] = broad1
			elseif((ewid1 < w[k]) && (w[k] < ewid2))
				allParams[j][2] = (broad1) +(((broad2)-(broad1))/(ewid2-ewid1))*(w[k]-ewid1)
			elseif(w[k] >= ewid2)
				allParams[j][2] = broad2
			endif
			
			allParams[j][3]  = y[k][0]	//TDMx
			allParams[j][4]  = y[k][1]	//TDMy	
			allParams[j][5]  = y[k][2]	//TDMz
			allParams[j][6]  = y[k][3]	//theta
			allParams[j][7]  = y[k][4]	//Transition Symmetry
			allParams[j][8]  = y[k][5]	//OSxx
			allParams[j][9]  = y[k][6]	//OSyy
			allParams[j][10] = y[k][7]	//OSzz
			
			allParams[j][11] = y[k][8]	//OSxy
			allParams[j][12] = y[k][9]	//OSxz
			allParams[j][13] = y[k][10]	//OSyz
			
			allParams[j][14] = i			//Atom ID
			if(StringMatch(fStep,"f0")==1)
				allParams[j][15] = LUMOpos	//Molecular Orbital	
				LUMOpos += 1
			else
				allParams[j][15] = z[k]
			endif
			allParams[j][16] = 0			//Transition Cluster
			k+=1
		endfor		
	endfor
	
	String pWaveAtomName
	for(i=0;i<=fnum-1;i+=1)
		if(StringMatch(fStep,"f0")==1)
			pWaveAtomName = "allParams"+num2str(i)
		elseif(StringMatch(fStep,"f1")==1)
			pWaveAtomName = "allParamsF1"+num2str(i)
		endif
		Wave w = $pWaveAtomName
		if(WaveExists(w))
			KillWaves/Z w
		endif
	endfor
	
	String sortedPWaveName 
	if(StringMatch(fStep,"f0")==1)
			sortedPWaveName = "allParamsSorted"
		elseif(StringMatch(fStep,"f1")==1)
			sortedPWaveName = "allParamsF1Sorted"
		endif
	Duplicate/O allParams,$sortedPWaveName
	KillWaves/Z allParams
	
	Wave sortedPWave = $sortedPWaveName
	
	MDSort(sortedPWave,0,rn=0)
	
	energyCorrection(sortedPWave,corrWave)	
	
	MDSort(sortedPWave,0,rn=0)
	
	print "Parameter Wave containing transition energies, oscillator strengths, TDM theta, and peak widths has been generated"
End

Function energyCorrection(pWave,corrWave)	

	Wave pWave,corrWave
	
	Variable i,j, x=DimSize(pWave,0),y=numpnts(corrWave)
	Variable corr , atomID, currentAtomID
	for(i=0;i<=y-1;i+=1)
		corr = corrWave[i]
		for(j=0;j<=x-1;j+=1)
			atomID = i+1
			currentAtomID = pWave[j][14]
			if(currentAtomID == atomID)
				pWave[j][0] = pWave[j][0] + corr
			endif
		endfor
	endfor
End

Function makeNEXAFS(pw,Eini,Efin,gRes,fStep)	
	
	Wave pw
	Variable Eini,Efin,gRes
	String fStep
	
	Variable nTrans=DimSize(pw,0),i,j
	
	String xxName,yyName,zzName,xyName,xzName,yzName,totalName,sumName
	xxName    = "xxTotal" + fStep
	yyName    = "yyTotal" + fStep
	zzName    = "zzTotal" + fStep
	xyName    = "xyTotal" + fStep
	xzName    = "xzTotal" + fStep
	yzName    = "yzTotal" + fStep
	sumName   = "sumSpec" + fStep
	totalName = "Total_Spec" + fStep 	
	Make/O/N=(gRes)  $xxName,$yyName,$zzName,$xyName,$xzName,$yzName,$sumName,$totalName
	Wave xxTotalf0 = $xxName,yyTotalf0 = $yyName,zzTotalf0 = $zzName,xyTotalf0 = $xyName,xzTotalf0 = $xzName,yzTotalf0 = $yzName,sumSpec = $sumName,Total_Spec = $totalName   
	SetScale/I x,Eini,Efin,"",xxTotalf0,yyTotalf0,zzTotalf0,xyTotalf0,xzTotalf0,yzTotalf0,sumSpec,Total_Spec 
	xxTotalf0 = 0;yyTotalf0 = 0;zzTotalf0 = 0;xyTotalf0 = 0;xzTotalf0 = 0;yzTotalf0 = 0;sumSpec=0;Total_Spec = 0  
	for(j=0;j<nTrans;j+=1)
		Variable pos,wid,amp,ampXX,ampYY,ampZZ,ampXY,ampXZ,ampYZ
				
		pos     = pw[j][0]
		wid     = pw[j][2]
		amp     = pw[j][1]
		ampXX   = pw[j][8]
		ampYY   = pw[j][9] 
		ampZZ   = pw[j][10]
		ampXY   = pw[j][11] 
		ampXZ   = pw[j][12]
		ampYZ   = pw[j][13]
		
		xxTotalf0  += ampXX *Gauss(x,pos,wid)
		yyTotalf0  += ampYY *Gauss(x,pos,wid)
		zzTotalf0  += ampZZ *Gauss(x,pos,wid)
		xyTotalf0  += ampXY *Gauss(x,pos,wid)
		xzTotalf0  += ampXZ *Gauss(x,pos,wid)
		yzTotalf0  += ampYZ *Gauss(x,pos,wid)
		Total_Spec += amp*Gauss(x,pos,wid)
	endfor
	
End

Function MDsort(w,keycol,[reversed,rn])
    
   Wave w	//Wave to be sorted
   variable keycol	//Column of wave w that is to serve as the sorting key
   Variable reversed	//If this parameter is set to 1 then the sorting is done from high to low. Else, sorting is done from low to high
   	Variable rn	//Rename wave?
   	
  	variable type   
   	type = Wavetype(w)
   	
 	if(rn==1)
	 	String nName = NameOfWave(w) + "Sorted"
		Duplicate/O w,$nName
		Wave wnew = $nName
	endif
	
	make/Y=(type)/free/n=(dimsize(w,0)) key
	make/free/n=(dimsize(w,0)) valindex 
   
   if(type == 0)      
       if(rn == 1)
       	Wave/t indirectSource = wnew
       else
       	Wave/t indirectSource = w
       endif
       
       Wave/t output = key
       output[] = indirectSource[p][keycol]
   else    
       if(rn == 1)
       	Wave indirectSource2 = wnew
       else
       	Wave indirectSource2 = w
       endif
       
       multithread key[] = indirectSource2[p][keycol]
   endif
   
   valindex=p
   if(reversed)
    sort/a/r key,key,valindex
   	else
  		Sort/A key,key,valindex
   	endif
   	
   if(type == 0)
       duplicate/free indirectSource, M_newtoInsert
       Wave/t output = M_newtoInsert
       output[][] = indirectSource[valindex[p]][q]
       indirectSource = output
   else
       duplicate/free indirectSource2, M_newtoInsert
       multithread M_newtoinsert[][] = indirectSource2[valindex[p]][q]
       multithread indirectSource2 = M_newtoinsert
   endif  
End

Function removeZeroAmplitude(pWave)

	Wave pWave
	
	Variable nTrans = DimSize(pWave,0),i
	
	for(i=nTrans-1;i>=0;i-=1)
		Variable amp = pWave[i][1]
		if(amp <= 0)
			DeletePoints i,1,pWave
		elseif(numtype(amp) !=0)
			DeletePoints i,1,pWave
		endif
	endfor
End

Function fetchClusterPWave(iniFolder,tval,ovpMax)
	
	String iniFolder
	Variable tval,ovpMax
	
	String pwFolder = iniFolder + ":firstFilter:cluster_" + replaceString(".",num2str(tval),"p") + "_OS_" + replaceString(".",num2str(ovpMax),"p") + "_OVP_ALL:"
	String ffFolder = iniFolder + ":firstFilter:"
	
	SetDataFolder pwFolder
	Wave combClusterPWAll
	String pwName = ffFolder + "combClusterPWAll"
	Duplicate/O combClusterPWAll,$pwName
End