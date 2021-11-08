#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function makeFittedOSWave(pWave,nPeaks)
	
	Wave pWave
	Variable nPeaks
	
	Variable i,j=1
	Make/O/N=(nPeaks) OS,En
	
	for(i=0;i<nPeaks;i+=1)
		En[i] = pWave[j]
		OS[i] = pWave[j+1]
		j+=16
	endfor
End

Function fitDFT_NXFS(pWave,yw,xw) : FitFunc
	Wave pWave //1d parameter wave containing the positions, widths, and amplitudes for al of the merged peaks
	Wave yw	//DFT NEXAFS to fit
	Wave xw	//Energy wave	
	
	yw = 0 
	Make/O/N=3 gpw
	Variable i
	for(i=1; i<numpnts(pWave); i+=16 )
		gpw[0] = pWave[i+1] //amplitude
		gpw[1] = pWave[i+0] //position
		gpw[2] = pWave[i+2] //width
		//yw+=sqrt(2*Pi)*gpw[0]*gpw[2]*gauss(xw,gpw[1],gpw[2])
		yw+=gpw[0]*gauss(xw,gpw[1],gpw[2])
	endfor
End

Function/WAVE makeMergeSpec(pWave,yw,xw)
	Wave pWave //1d parameter wave containing the positions, widths, and amplitudes for al of the merged peaks
	Wave yw	//DFT NEXAFS to fit
	Wave xw	//Energy wave	
	
	yw = 0 
	Make/O/N=3 gpw
	Variable i
	for(i=1; i<numpnts(pWave); i+=16 )
		gpw[0] = pWave[i+1] //amplitude
		gpw[1] = pWave[i+0] //position
		gpw[2] = pWave[i+2] //width
		yw+=gpw[0]*gauss(xw,gpw[1],gpw[2])
	endfor
	
	return yw
End

Function getAmpFitWaves(baseFolderName,tval,ovpMax)

	String baseFolderName
	Variable tval,ovpMax
	
	SetDataFolder $baseFolderName
	String newAmpFitFolder = "AmplitudeFitting"
	NewDataFolder/O $newAmpFitFolder
	
	//Location of Original/Unfiltered DFT NEXAFS
	String pathToOriSpecs = baseFolderName + ":originalPeaks:Spectra:tdmCompSpecs:"
	
	//Location of Merged/OS-Symmetry Filtered DFT NEXAFS
	String allClusterFolder = "cluster_" + replacestring(".",num2str(tval),"p") + "_OS_" + replacestring(".",num2str(ovpMax),"p") + "_OVP_All:"
	
	String mergeParamsLoc
	String pathToAllmergeParams
	mergeParamsLoc = baseFolderName + ":firstFilter:"
	pathToAllmergeParams = mergeParamsLoc + allClusterFolder
	
	Variable i
	//This will transfer the unfiltered DFT NEXAFS into the directory where the amplitude fitting will take place	
	SetDataFolder $pathToOriSpecs
	String ampFitLoc = ":AmplitudeFitting:"
	
	String oriSpecList = SortList(WaveList("*f0",";",""))
	Variable nOriSpec = ItemsInList(oriSpecList)
	for(i=0;i<=nOriSpec-1;i+=1)
		String currOriSpecName = StringFromList(i,oriSpecList)
		String oriSpecName = baseFolderName + ampFitLoc + StringFromList(i,oriSpecList)  
		Wave w = $currOriSpecName
		if(WaveExists(w))
			Duplicate/O w,$oriSpecName
		endif
	endfor
	
	//This defines the names of the 1D merged parameter waves
	String specPathList = pathToAllmergeParams
	String mergeAll = "combClusterPWAll"//"mergedParams_"  + replacestring(".",num2str(tval),"p") + "OS_" + replacestring(".",num2str(ovpMax),"p") + "_OVP_All"

	String specWaveList = mergeAll
	String specWaveList2 = ""
	String specPathList2 = ""
	String elemList = "ALL"
	String elemList2 = ""
	Variable j=0
	for(i=0;i<1;i+=1)
		String currFolder = StringFromList(i,specPathList)
		if(DataFolderExists(currFolder))
			specWaveList2 = SortList(AddListItem(StringFromList(i,specWaveList),specWaveList2,";"))
			specPathList2 = SortList(AddListItem(StringFromList(i,specPathList),specPathList2,";"))
			elemList2     = SortList(AddListItem(StringFromList(i,elemList),elemList2,";"))
			SetDataFolder $currFolder
			String currWave =  StringFromList(i,specWaveList)
			String transferName =   baseFolderName + ampFitLoc + replaceString(";",StringFromList(j,specWaveList2),"")
			j+=1
			if(WaveExists($currWave))
				Duplicate/O $currWave,$transferName
			endif
		endif
	endfor
	
	SetDataFolder $baseFolderName
	SetDataFolder $newAmpFitFolder
	//Organize Waves		
	j=0
	Variable nSymmetries = ItemsInList(specWaveList2)
	String fNameList =""
	for(i=0;i<nSymmetries;i+=1)
		String listItem = StringFromList(i,specWaveList2)
		fNameList = SortList(AddListItem(listItem,fNameList,";"))
	endfor
	//print fNameList
	for(i=0;i<nSymmetries;i+=1)
		String fName   = ReplaceString("combClusterPW",StringFromList(i,fNameList),"")
		String oriName = StringFromList(i,oriSpecList)
		String merName = StringFromList(j,specWaveList2)
		String fName2  = StringFromList(j,specPathList2)
		String cElem2  = StringFromList(j,elemList2)
		if(WaveExists($oriName))
			NewDataFolder/O $fName
			String oriSpecNameDuplicate =  baseFolderName + ampFitLoc + fName  + ":" + StringFromList(i,oriSpecList)
			Duplicate/O $oriName,$oriSpecNameDuplicate
		endif
		if(WaveExists($merName))
			String merSpecNameDuplicate =  baseFolderName + ampFitLoc + fName  + ":" + replaceString(";",StringFromList(j,specWaveList2),"")
			Duplicate/O $merName,$merSpecNameDuplicate
			SetDataFolder $fName2
			SetDataFolder $baseFolderName
			SetDataFolder $newAmpFitFolder
			j+=1
		endif
		KillWaves/Z $oriName,$merName
	endfor
	
	SetDataFolder $baseFolderName
	
End

Function removeUselessPeaks(pWave)

	Wave pWave
	
	Variable tol = 1e-6
	Variable n = numpnts(pWave)
	Variable i
	for(i=n-1;i>=0;i-=3)
		Variable amp = pWave[i]
		if(amp <= tol)
			DeletePoints i-2,3,pWave
		endif
	endfor

End

Function calcPercentDiff(oriSpec,mergeFitSpec)
	Wave oriSpec,mergeFitSpec
	
	Variable percentDiff,num,den,i
	
	for(i=0;i<=numpnts(oriSpec)-1;i+=1)
		num += abs(oriSpec[i]-mergeFitSpec[i])
		den += oriSpec[i]
	endfor
	
	percentDiff = (num/den)*100
	return percentDiff
End

Function ampFittingWrapper(baseFolderName,tval,ovpMax,Eini,Efin,gRes)
	
	String baseFolderName
	Variable tval,ovpMax,Eini,Efin,gRes
	
	//Get and organize the waves necessary for the amplitude fitting
	getAmpFitWaves(baseFolderName,tval,ovpMax)
	
	String ampDirectory = baseFolderName + ":AmplitudeFitting"

	SetDataFolder $ampDirectory
	Variable i
	
	String folderList = ReplaceString("\r",ReplaceString(",",ReplaceString(";FOLDERS:",SortList(dataFolderDir(1)),""),";"),"")///"ALL;XX;XY;XZ;YY;YZ;ZZ"
	Variable nFolders = ItemsInList(folderList)
	
	for(i=0;i<=nFolders;i+=1)
		String currentFolder = StringFromList(i,folderList)
		SetDataFolder $currentFolder
		String perDiffName = "percentDiff_" + currentFolder
		Make/O/N=1 $perDiffName
		Wave perDiffWave = $perDiffName
		String iniWaveName = WaveList("*f0",";","")
		String iniName = StringFromList(0,iniWaveName)
		Wave spec = $iniName
		String pWaveName = WaveList("combCluster*",";","")//WaveList("mergedParams*",";","")
		String pName = StringFromList(0,pWaveName)
		Wave pWave = $pName
		if(WaveExists(pWave))
			Variable perDiff = fitAmplitudes2(pWave,spec,Eini,Efin,currentFolder,tval,ovpMax)//fitAmplitudes(pWave,spec,tval, ovpMax,currentFolder,Eini,Efin,gRes,maxTransCluster,pass)
			perDiffWave[0] = perDiff
		endif
		SetDataFolder $ampDirectory
	endfor
		
	//Kill extra unnecessary folders	
	SetDataFolder $baseFolderName
	String fList = ReplaceString("\r",ReplaceString(",",ReplaceString(";FOLDERS:",SortList(dataFolderDir(1)),""),";"),"")
	Variable nf = ItemsInList(fList)
	for(i=0;i<nf;i+=1)
		String cf = StringFromList(i,fList)
		if(stringmatch(cf,"cluster*"))
			KillDataFolder $cf
		endif
	endfor
	SetDataFolder $ampDirectory
	
	return perDiff
End

Function make2DpWave(pw1d,pw2d,tval,ovpVal,sym)

	Wave pw1d	 //This is the 1d merged params wave that has the peak position, width and fitted amplitude 
	Wave pw2d //This is the 2d parameter wave that contains the information for the max clusters
	Variable tval,ovpVal
	String sym
	
	Variable n1D = numpnts(pw1d)
	Variable n2D = DimSize(pw2D,0)
	
	Variable nPeaks = n1D/3
	String iniPwaveName = "mergedParams_" + replacestring(".",num2str(tval),"p") + "OS_" + replacestring(".",num2str(ovpVal),"p") + "_OVP_" + sym + "_2"

	Duplicate/O  pw2D, $iniPwaveName
	Wave w = $iniPwaveName
	Variable i,j
	for(i=0,j=0;i<n2D && j<n1D;i+=1 , j+=3)
		w[i][0] = pw1D[j]	//Peak position
		w[i][1] = pw1D[j+2]	//Peak Width
		w[i][2] = pw1D[j+1]	//Peak Amplitude
	endfor
End

Function/S makeFitWaves(pw2d,npks,sym)

	Wave pw2d
	Variable npks
	String sym
	
	Variable i,j=0,eps0,eps1,eps2,eps3,eps4,eps5,eps6,cVal,symVal
	String pwName  = "fitPW" + sym
	String epsName = "eps"   + sym
	String conName = "constraints" + sym
	String hName = "hold"  + sym
	String H = "1"
	Make/O/N=(16*npks+1) $pwName,$epsName,$hName
	Make/O/T/N=(npks) $conName
	Wave   w   = $pwName
	Wave   eps = $epsName
	Wave/T con = $conName
	Wave hw = $hName
	String hold = "1"
	//Define Epsilon Wave Values
	eps0 = 1e-4
	cVal = 1
	H = "1011111111111111"
	symVal = 0

	w[0] = symVal
	eps[0] = 0
	//Make 1D PW,Epsilon Wave, and Constraint Wave
	for(i=1;i<16*npks+1;i+=16)
		w[i+0]  = pw2d[j][0]
		w[i+1]  = pw2d[j][1]
		w[i+2]  = pw2d[j][2]
		w[i+3]  = pw2d[j][3]
		w[i+4]  = pw2d[j][4]
		w[i+5]  = pw2d[j][5]
		w[i+6]  = pw2d[j][6]
		w[i+7]  = pw2d[j][7]
		w[i+8]  = pw2d[j][8]
		w[i+9]  = pw2d[j][9]
		w[i+10] = pw2d[j][10]
		w[i+11] = pw2d[j][11]
		w[i+12] = pw2d[j][12]
		w[i+13] = pw2d[j][13]
		w[i+14] = pw2d[j][14]
		w[i+15] = pw2d[j][15]
		
		eps[i+0]  = 0
		eps[i+1]  = 0
		eps[i+2]  = eps0
		eps[i+3]  = 0
		eps[i+4]  = 0
		eps[i+5]  = 0
		eps[i+6]  = 0
		eps[i+7]  = 0
		eps[i+8]  = 0
		eps[i+9]  = 0
		eps[i+10] = 0
		eps[i+11] = 0
		eps[i+12] = 0
		eps[i+13] = 0
		eps[i+14] = 0
		eps[i+15] = 0
		
		String lowConstraint   = "K" + num2str(i+cVal) + " >= 0" 
		con[j] = lowConstraint
		hold += H
		j+=1
	endfor
	
	Variable hlen = strlen(hold)
	for(i=0;i<16*npks+1;i+=1)
		hw[i] = str2num(hold[i]) 
	endfor
	
	return hold
End

Function fitAmplitudes2(pWave,spec2Fit,Eini,Efin,sym,tval,ovpMax)

	Wave pWave,spec2fit
	Variable Eini,Efin
	String sym
	Variable tval,ovpMax
	
	String fitSpecName = "fitSpec" + sym
	String merSpecName = "mergedOri" + sym
	
	Duplicate/O spec2fit,$fitSpecName,res,energy,$merSpecName
	Duplicate/O pwave, PW2DOriginal
	Wave fitSpec = $fitSpecName
	fitspec=0
	Redimension/D pWave
	Wave pWaveOriginal
	if(!WaveExists(pWaveOriginal))
		Duplicate/O pWave,pWaveOriginal
	endif
	
	Variable i,gres = abs((Efin-Eini)/2000)
	
	//Define Hold String,Epsilon wave, constraint wave. Make 1D parameter wave from 2D
	NVAR numPeaks = root:numPeaks
	numPeaks = DimSize(pWave,0)
	String H = makeFitWaves(pWave,numPeaks,sym)//Hold string for peaks. Keep position and width constant. Vary the amplitude.
	String pwName  = "fitPW" + sym
	String epsName = "eps"   + sym
	String conName = "constraints" + sym
	Wave   pw1D   = $pwName
	Wave   eps = $epsName
	Wave/T con = $conName
	Duplicate/O pw1d,pw1dOriginal

	//Make x wave
	for(i=0;i<2000;i+=1)
		energy[i] = gres*i + Eini
	endfor
	
	Variable V_FitTol = 0.00001
	Variable V_FitError = 0

	FuncFit/Q/H=H fitRandGauss , pw1D, spec2fit /X=energy /D=fitSpec /R=res /E=eps /C=con	
	placeFitAmpIn2DPW(pw1D,pWave)
	
	//Calculate percent difference between original NEXAFS and merged DFT
	//NVAR perDiff = root:perDiff
	//perDiff = calcPercentDiff(spec2fit,fitSpec)	
	
	//Make the merged Spectrum using the original parameter wave
	Wave mergeSpec = $merSpecName
	Wave mergeSpec2 = makeMergeSpec(pw1dOriginal,mergeSpec,energy)
	NVAR perDiff = root:perDiff
	perDiff = calcPercentDiff(spec2fit,mergeSpec2)	
	//Make the 1D fitted Oscillator Strength wave to add to the final plot
	makeFittedOSWave(pw1D,numPeaks)
	Wave OS, En	
	
	//Plot the target spectra, the fit, the original merged spectrum and the residuals
	//String graphName = "ampFit_OS" + replacestring(".",num2str(tval),"p") + "_OVP" + replacestring(".",num2str(ovpMax),"p") + "_" + sym// + "_" + num2str(pass)
	//String targetSpecName = NameOfWave(spec2fit)
	//DoWindow $graphName
	//if(!V_Flag)
	//	Display/N=$graphName/W=(0,0,400,500)/K=1 spec2fit,fitSpec,mergeSpec vs energy
	//	Label left sym + "Transition Intensity [a.u.]\\U"
	//	Label bottom "Transition Energy[eV]"
	//	NewFreeAxis/L residual
	//	NewFreeAxis/L oss
	//	ModifyGraph axisEnab(oss)={0,0.2},axisEnab(left)={0.25,0.75},axisEnab(residual)={0.8,1},freePos(residual)=0, freePos(oss)=0,lblPosMode=1
	//	AppendToGraph/L=residual res
	//	AppendToGraph/L=oss OS vs En
	//	ModifyGraph lsize=1.5,lstyle($fitSpecName)=3,rgb($fitSpecName)=(0,0,0),rgb(res)=(52428,34958,1),rgb($merSpecName)=(2,39321,1),grid=2,mirror=1,nticks=5,minor=1,fStyle=1,lsize($fitSpecName)= 2.5
	//	ModifyGraph mode(OS)=1,rgb(OS)=(1,9611,39321),mirror(bottom)=2,log(oss)=1
	//	WaveStats/Q res
	//	SetAxis residual V_min,V_max
	//	Label residual "Residual\\U"
	//	WaveStats/Q OS
	//	SetAxis oss 0.001,1
	//	Label oss "OS [a.u.]\\U"
	//	Legend/C/N=text0/J/A=RC/X=4.11/Y=-1.88 "OS = " + num2str(tval) + "% OVP=" + num2str(ovpMax) + "%\r\\s("+ targetSpecName + ") Unfiltered\r\\s("+fitSpecName+") Fit\r\\s("+merSpecName+") Merged\r\\s(res) Residual\r% Diff = " + num2str(perDiff) +"%\r\\s(OS) OS"
	//endif
	return perDiff
End

Function fitRandGauss(pWave,yw,xw) : FitFunc
	Wave pWave //2d parameter wave containing the positions, widths, and amplitudes for al of the merged peaks
	Wave yw	//DFT NEXAFS to fit
	Wave xw	//Energy wave	
	
	yw = 0 
	Make/O/N=3 gpw
	Variable i
	for(i=1; i<numpnts(pWave)+0; i+=16)
		gpw[0] = pWave[i+1] //amplitude
		gpw[1] = pWave[i+0] //position
		gpw[2] = pWave[i+2] //width
		yw+= sqrt(2*Pi)*gpw[0]*gpw[2]*gauss(xw,gpw[1],gpw[2])
	endfor
End