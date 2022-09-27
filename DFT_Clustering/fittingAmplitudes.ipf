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

Function/S makeFitWaves(pw2d,npks,sym,holdPos,holdWid,holdAmp,stage)

	Wave pw2d
	Variable npks,holdPos,holdWid,holdAmp,stage
	String sym
	
	Variable i,j=0,eps0,eps1,eps2,eps3,eps4,eps5,eps6,cVal,symVal
	String pwName  = "fitPW" + sym + "_" + num2str(stage)
	String epsName = "eps"   + sym+ "_" + num2str(stage)
	String conName = "constraints" + sym + "_" + num2str(stage)
	String hName = "hold"  + sym + "_" + num2str(stage)
	String H = "1"
	Make/O/N=(16*npks+1) $pwName,$epsName,$hName
	Variable openParams = 4 - 2*holdPos - holdWid - holdAmp,pLo,pHi,k=0
	Make/O/T/N=(openParams*npks) $conName
	Wave   w   = $pwName
	Wave   eps = $epsName
	Wave/T con = $conName
	Wave hw = $hName
	String hold = "1"
	//Define Epsilon Wave Values
	eps0 = 1e-6
	cVal = 1
	if(!holdPos && !holdWid && !holdAmp)
		H = "0001111111111111"
	elseif(holdPos && !holdWid && !holdAmp)
		H = "1001111111111111"
	elseif(holdPos && holdWid && !holdAmp)
		H = "1011111111111111"
	endif
	symVal = 0

	w[0] = symVal
	eps[0] = 0
	String lowAmpConstraint,lowWidConstraint,loPos,hiPos
	//Make 1D PW,Epsilon Wave, and Constraint Wave
	for(i=1;i<16*npks+1;i+=16)
		w[i+0]  = pw2d[k][0]
		w[i+1]  = pw2d[k][1]
		w[i+2]  = pw2d[k][2]
		w[i+3]  = pw2d[k][3]
		w[i+4]  = pw2d[k][4]
		w[i+5]  = pw2d[k][5]
		w[i+6]  = pw2d[k][6]
		w[i+7]  = pw2d[k][7]
		w[i+8]  = pw2d[k][8]
		w[i+9]  = pw2d[k][9]
		w[i+10] = pw2d[k][10]
		w[i+11] = pw2d[k][11]
		w[i+12] = pw2d[k][12]
		w[i+13] = pw2d[k][13]
		w[i+14] = pw2d[k][14]
		w[i+15] = pw2d[k][15]
		
		if(holdPos)
			eps[i+0]  = 0
		else
			eps[i+0]  = eps0
		endif
		
		if(holdAmp)
			eps[i+1]  = 0
		else
			eps[i+1]  = eps0
		endif
		
		if(holdWid)	
			eps[i+2]  = 0
		else
			eps[i+2]  = eps0
		endif
		
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
		
		if(!holdPos && !holdWid && !holdAmp)
			pLo = w[i+0] - w[i+2]
			pHi = w[i+0] + w[i+2]
			lowAmpConstraint = "K" + num2str(i+cVal) + " > 0"
			lowWidConstraint = "K" + num2str(i+cVal+1) + " > 0"
			loPos = "K" + num2str(i) + " > " + num2str(pLo)
			hiPos = "K" + num2str(i) + " < " + num2str(pHi)
			con[j]   = loPos
			con[j+1] = hiPos
			con[j+2] = lowAmpConstraint
			con[j+3] = lowWidConstraint
			j+=openParams
		elseif(holdPos && !holdWid && !holdAmp)
			lowAmpConstraint = "K" + num2str(i+cVal) + " > 0"
			lowWidConstraint = "K" + num2str(i+cVal+1) + " > 0"
			con[j]   = lowAmpConstraint
			con[j+1] = lowWidConstraint
			j+=openParams
		elseif(holdPos && holdWid && !holdAmp)
			lowAmpConstraint   = "K" + num2str(i+cVal) + " > 0"
			con[j] = lowAmpConstraint
			j+=openParams
		endif
		hold += H
		k+=1
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
	Duplicate/O pwave,pwave2,pwave3,PW2DOriginal,pw2DFit
	Wave fitSpec = $fitSpecName
	fitspec=0
	Redimension/D pWave
	Wave pWaveOriginal
	if(!WaveExists(pWaveOriginal))
		Duplicate/O pWave,pwave2,pwave3,pWaveOriginal
	endif
	
	Variable i,gres = abs((Efin-Eini)/2000)
	
	//Define Hold String,Epsilon wave, constraint wave. Make 1D parameter wave from 2D
	NVAR numPeaks = root:numPeaks
	numPeaks = DimSize(pWave,0)
	String H = makeFitWaves(pWave,numPeaks,sym,1,1,0,1)//Hold string for peaks. Keep position and width constant. Vary the amplitude.
	String pwName1  = "fitPW" + sym + "_1"
	String epsName1 = "eps"   + sym + "_1"
	String conName1 = "constraints" + sym + "_1"
	
	String pwName2  = "fitPW" + sym + "_2"
	String epsName2 = "eps"   + sym + "_2"
	String conName2 = "constraints" + sym + "_2"
	
	String pwName3  = "fitPW" + sym + "_3"
	String epsName3 = "eps"   + sym + "_3"
	String conName3 = "constraints" + sym + "_3"
	Wave   pw1D   = $pwName1
	Wave   eps = $epsName1
	Wave/T con = $conName1
	Duplicate/O pw1d,pw1dOriginal

	//Make x wave
	for(i=0;i<2000;i+=1)
		energy[i] = gres*i + Eini
	endfor
	
	Variable V_FitTol = 0.00001,V_FitError = 0

	FuncFit/Q/H=H fitRandGauss , pw1D, spec2fit /X=energy /D=fitSpec /R=res /E=eps /C=con
	Variable num=0,den=0,chisq=0
	NVAR chiSqBBCl = root:chiSqBBCl
	NVAR GFBBCl = root:GFBBCl
	
	NVAR chiSqBBCluref = root:chiSqBBCluref
	NVAR GFBBCluref = root:GFBBCuref
	
	GFBBCl = 0 
	GFBBCluref = 0
	//Calculte GF for refined cluster vs DFT
	for(i=0;i<2000;i+=1)
		num = (spec2fit[i]-fitSpec[i])^2
		den = spec2fit[i]
		chiSq += num/den
		GFBBCl +=num/(den^2)
	endfor
	chiSqBBCl = chiSq
	NVAR redchiSqBBCl = root:redchiSqBBCl
	NVAR redchiSqBBCluref = root:redchiSqBBCluref
	redchiSqBBCl = chiSq/(2000-numPeaks)
	print chiSqBBCl,redchiSqBBCl
	placeFitAmpIn2DPW(pw1D,pw2DFit)
	//Make the merged Spectrum using the original parameter wave
	Wave mergeSpec = $merSpecName
	Wave mergeSpec2 = makeMergeSpec(pw1dOriginal,mergeSpec,energy)
	//Calculte GF for refined cluster vs DFT
	num=0;den=0;chiSqBBCluref=0
	for(i=0;i<2000;i+=1)
		num = (spec2fit[i]-mergeSpec2[i])^2
		den = spec2fit[i]
		chiSqBBCluref += num/den
		GFBBCluref +=num/(den^2)
	endfor
	redchiSqBBCluref = chiSqBBCluref/(2000-numPeaks)
	NVAR perDiff = root:perDiff
	perDiff = calcPercentDiff(spec2fit,fitSpec)
	//perDiff = calcPercentDiff(spec2fit,mergeSpec2)
	//Make the 1D fitted Oscillator Strength wave to add to the final plot
	makeFittedOSWave(pw1D,numPeaks)
	Wave OS, En	
	//plotRawDFTvsClusterDFT(0,83,"CuPc")
	return perDiff
End

Function calcGF(expW,model,n)
	Wave expW,model
	Variable n
	
	Variable i,GF=0,num=0,den=0,den2=0,chiSq=0
	for(i=0;i<n;i+=1)
		num = (expW[i]-model[i])^2
		den = model[i]
		den2 = expW[i]^2
		chiSq += num/den
		GF +=num/den2
	endfor
	print GF
End

Function/WAVE placeFitAmpIn2DPW(pw1,pw2)

	Wave pw1,pw2
		
	Variable i,n= DimSize(pw2,0),j=2
	for(i=0;i<n;i+=1)
		pw2[i][1] = pw1[j]
		j+=16
	endfor
	
	return pw2
End

Function/WAVE placeFitWidIn2DPW(pw1,pw2)

	Wave pw1,pw2
		
	Variable i,n= DimSize(pw2,0),j=3
	for(i=0;i<n;i+=1)
		pw2[i][2] = pw1[j]
		j+=16
	endfor
	
	return pw2
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
