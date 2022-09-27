#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function dftPlots(baseFolderName,tval,ovpMax)
	
	String baseFolderName
	Variable tval,ovpMax
	
	//Define paths to relevant data folders containing the parameter waves for each filtering stage
	String iniFolder = GetDataFolder(1)
	String pathToOriginal    = baseFolderName + ":originalPeaks:"
	String pathToFirstFilter = baseFolderName + ":firstFilter:"
	String pathToAllClusters = baseFolderName + ":firstFilter:cluster_" + replaceString(".",num2str(tval),"p") + "_OS_" + replaceString(".",num2str(ovpMax),"p") + "_OVP_All:"
	
	//Find the parameter waves for each of the filtering/clustering steps
	SetDataFolder $pathToOriginal
	Wave allParamsSorted
	
	SetDataFolder $pathToFirstFilter
	Wave allParamsF1Sorted	
	
	String pathToMergedParams = pathToAllCLusters + "mergedParams:"
	String merged1DPwaveName = "mergedParamsWave_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p") + "_OVP_All"
	SetDataFolder $pathToMergedParams
	Wave mergedPwave1D = $merged1DPwaveName
	Variable nClusters = numpnts(mergedPWave1D)/3
	Make/O/N=(nClusters,3) mergedPWave2D
	Variable i,j=0
	
	for(i=0;i<3*nClusters;i+=3)
		mergedPWave2D[j][0] = mergedPWave1D[i]
		mergedPWave2D[j][1] = mergedPWave1D[i+1]
		mergedPWave2D[j][2] = mergedPWave1D[i+2]
		j+=1
	endfor
	
	Wave mergedPWave2D
	WaveStats/Q/RMD=[][1,1] mergedPWave2D
	SetDataFolder iniFolder
	
	//Plot the oscillator strengths for each of the filtering stages
	String osPlot = "OscillatorStrengths_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	
	DoWindow $osPlot
	if(!V_Flag)
		Display/N=$osPlot/K=1 allParamsSorted[*][1] vs allParamsSorted[*][0]
		NewFreeAxis/L f1
		NewFreeAxis/L f2
		AppendtoGraph/L=f1/W=$osPlot allParamsF1Sorted[*][1] vs allParamsF1Sorted[*][0]
		AppendtoGraph/L=f2/W=$osPlot mergedPWave2D[*][2] vs mergedPWave2D[*][0]
		Label left "ƒ\\U";DelayUpdate
		Label bottom "Energy[eV]";DelayUpdate
		Label f1 "ƒ\\U";DelayUpdate
		Label f2 "ƒ\\U";DelayUpdate
		ModifyGraph/W=$osPlot mode=1,tick(left)=2,tick(f1)=2,zero(f1)=1,mirror=1,axisEnab(left)={0,0.30},axisEnab(f1)={0.34,0.64},freePos(f1)=0,lblPosMode(f1)=2
		ModifyGraph/W=$osPlot zero(f2)=1,axisEnab(f2)={0.7,1},freePos(f2)=0,lblPosMode(f2)=2,lsize=1.5
		ModifyGraph/W=$osPlot rgb(allParamsF1Sorted)=(0,0,0)
		ModifyGraph/W=$osPlot rgb(mergedPWave2D)=(1,4,52428)

		SetAxis left *,0.02;DelayUpdate
		SetAxis f1 *,0.02
		SetAxis f2 *,0.2
		
		String legendText = "\\s(allParamsSorted) OS Filter\r\\s(allParamsF1Sorted) OVP Filter\r\\s(mergedPWave2D) Clustering"

		SetAxis bottom 284,320;DelayUpdate
		Legend/C/N=text0/A=MC/X=-36/Y=40.00 legendText
	endif
	
	//Paths to various TOTAL NEXAFS
	String pathToOrigNXFS   = pathToOriginal    + "Spectra:tdmCompSpecs:"
	String pathToFFNXFS     = pathToFirstFilter + "Spectra:tdmCompSpecs:"
	String pathToMergedNXFS = baseFolderName    + ":AmplitudeFitting:ALL:" 
	
	SetDataFolder $pathToOrigNXFS
	Wave Total_Specf0
	
	SetDataFolder $pathToFFNXFS
	Wave Total_Specf1
	
	SetDataFolder $pathToMergedNXFS
	Wave fitSpec
	
	//Plot the TOTAL NEXAFS spectrum obtained from each of the filtering stages
	String totalPlotWindow = "Total_NEXAFS_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	DoWindow $totalPlotWindow
	if(!V_Flag)
		Display/W=(0,0,500,300)/N=$totalPlotWindow/K=1 Total_Specf0,Total_Specf1,fitSpec
		ModifyGraph lsize=2,rgb(Total_Specf0)=(0,0,0),rgb(fitSpec)=(1,4,52428)
		ModifyGraph grid=2,mirror=1,nticks=10,minor=1,fStyle=1,fSize=16;DelayUpdate
		Label left "Transition Intensity [a.u.] \\U";DelayUpdate
		Label bottom "Transiton Energy [eV]";DelayUpdate
		SetAxis bottom 283,320
		Legend/C/N=text0/J/A=RT "\\s(Total_Specf0) Original DFT \r\\s(Total_Specf1) OS Filter\r\\s(fitSpec) Clustered"
	endif
	
	//Plot the overlap Matrix for the original DFT and the first filter
	SetDataFolder $pathToOriginal
	Wave ovpWaveAll,oriAllPeaks2D
	
	String osOVPWin = "OS_OVP_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	DoWindow $osOVPWin
	if(!V_Flag)
		Display/N=$osOVPWin/W=(0,0,730,600)/K=1
		AppendImage ovpWaveAll
		ModifyImage/W=$osOVPWin ovpWaveAll ctab={*,*,ColdWarm,0}
		ColorScale/C/N=text0/A=MC image=ovpWaveAll;DelayUpdate
		ColorScale/C/N=text0 axisRange={0,100}
		ColorScale/C/N=text0 "Peak Overlap[%]"
		ColorScale/C/N=text0/X=60.00/Y=5.00
		ModifyGraph fsize=16,fstyle=1;DelayUpdate
		Label/W=$osOVPWin left "Peak ID";DelayUpdate
		Label/W=$osOVPWin bottom "Peak ID";DelayUpdate
		ModifyGraph/W=$osOVPWin margin(left)=72,margin(bottom)=63,margin(right)=106,gfsize=16
		ModifyGraph nticks=10,minor=1
	endif
	
	SetDataFolder $pathToFirstFilter
	Wave ovpWaveAll,f1AllPeaks2D
	String symOVPWin = "SymOVP_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	DoWindow $symOVPWin
	if(!V_Flag)
		Display/N=$symOVPWin/W=(0,0,730,600)/K=1
		AppendImage ovpWaveAll
		ModifyImage/W=$symOVPWin ovpWaveAll ctab={*,*,ColdWarm,0}
		ColorScale/C/N=text0/A=MC image=ovpWaveAll;DelayUpdate
		ColorScale/C/N=text0 axisRange={0,100}
		ColorScale/C/N=text0 "Peak Overlap[%]"
		ColorScale/C/N=text0/X=60.00/Y=5.00
		ModifyGraph fsize=16,fstyle=1;DelayUpdate
		Label/W=$symOVPWin left "Peak ID";DelayUpdate
		Label/W=$symOVPWin bottom "Peak ID";DelayUpdate
		ModifyGraph/W=$symOVPWin margin(left)=63,margin(bottom)=72,margin(right)=106,gfsize=16
		ModifyGraph nticks=10,minor=1
	endif
	
	//Plot 2d Matrix of Transition peaks for oriignal DFT and first filter
	String ori2DGaussWin = "orig2DGauss_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	String f12DGaussWin  = "f12DGauss_"   + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	
	DoWindow $ori2DGaussWin
	if(!V_Flag)
		Display/N=$ori2DGaussWin/W=(0,0,730,300)/K=1
		AppendImage oriAllPeaks2D
		ModifyImage/W=$ori2DGaussWin oriAllPeaks2D ctab={*,0.001,SpectrumBlack,0}
		ColorScale/C/N=text0/A=MC image=oriAllPeaks2D;DelayUpdate
		ColorScale/C/N=text0 axisRange={0,0.0001}
		ColorScale/C/N=text0 "Peak Intensity\U"
		ColorScale/C/N=text0/X=60.00/Y=5.00
		ModifyGraph fsize=16,fstyle=1;DelayUpdate
		Label/W=$ori2DGaussWin left "Transition Energy[eV]";DelayUpdate
		Label/W=$ori2DGaussWin bottom "Peak ID";DelayUpdate
		ModifyGraph/W=$ori2DGaussWin margin(left)=63,margin(bottom)=63,margin(right)=106,gfsize=16
		ModifyGraph nticks=10,minor=1
	endif
	
	DoWindow $f12DGaussWin
	if(!V_Flag)
		Display/N=$f12DGaussWin/W=(0,0,730,300)/K=1
		AppendImage f1AllPeaks2D
		ModifyImage/W=$f12DGaussWin f1AllPeaks2D ctab={*,0.001,SpectrumBlack,0}
		ColorScale/C/N=text0/A=MC image=f1AllPeaks2D;DelayUpdate
		ColorScale/C/N=text0 axisRange={0,0.0001}
		ColorScale/C/N=text0 "Peak Intensity\U"
		ColorScale/C/N=text0/X=60.00/Y=5.00
		ModifyGraph fsize=16,fstyle=1;DelayUpdate
		Label/W=$f12DGaussWin left "Transition Energy[eV]";DelayUpdate
		Label/W=$f12DGaussWin bottom "Peak ID";DelayUpdate
		ModifyGraph/W=$f12DGaussWin margin(left)=63,margin(bottom)=63,margin(right)=106,gfsize=16
		ModifyGraph nticks=10,minor=1
	endif
	
	//Make peaks for different filtering stages from 2D Gauss
	Variable xOri = DimSize(oriAllPeaks2D,0)
	Variable xF1  = DimSize(f1AllPeaks2D,0)
	String peaksWin = "peaks_" + replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p")
	DoWindow $peaksWin
	if(!V_Flag)
		Display/N=$peaksWin/K=1 
		NewFreeAxis/L f1
		
		for(i=0;i<=xOri-1;i+=1)
			AppendToGraph/W=$peaksWin oriAllPeaks2D[i][*]	//Original Peaks
			String oriPkName = NameOfWave(oriAllPeaks2D) + "#" + num2str(i)
		endfor
		
		for(i=0;i<=xF1-1;i+=1)
			AppendToGraph/L=f1/W=$peaksWin f1AllPeaks2D[i][*]	//OS filtered peaks
			String f1PkName = NameOfWave(f1AllPeaks2D) + "#" + num2str(i)
			ModifyGraph/W=$peaksWin rgb($f1PkName)=(0,0,0)
		endfor
		
		ModifyGraph mirror(bottom)=1, axThick=2,nticks=10,fstyle=1,fsize=14,minor=1,lsize=1,mirror(left)=1,mirror(f1)=1
		ModifyGraph axisEnab(left)={0,0.49},axisEnab(f1)={0.51,1},freePos(f1)=0,lblPosMode(f1)=2,lsize=1.5
		SetAxis left 0,0.02
		SetAxis f1 0,0.02
		SetAxis bottom 280,320
		Label left "Intensity\\U";DelayUpdate
		Label bottom "Transition Energy[eV]";DelayUpdate
		Label f1 "Intensity\\U";DelayUpdate
	endif
	
	SetDataFolder iniFolder	
End

Function plotTDMcomps(aName,fnum)

	String aName
	Variable fnum
	
	String iniFolder = GetDataFolder(1)
	Variable i
	for(i=1;i<=fnum;i+=1)
		String fName = aName + num2str(i)
		SetDataFolder $fName	
		Wave TDMx_Ori
		Wave TDMy_Ori
		Wave TDMz_Ori
		Wave symmetryParams
		Wave eV_
		String gNameX = fName + "_TDMx"
		String gNameY = fName + "_TDMy"
		String gNameZ = fName + "_TDMz"
		DoWindow $gNameX
		if(!V_Flag)
			Display/N=$gNameX/K=1 TDMx_Ori vs eV_
			NewFreeAxis/L symRot
			AppendToGraph/L=symRot symmetryParams[*][0] vs eV_
			Label left "μ\\Bx [a.u.] \\U";DelayUpdate
			Label symRot "μ\\Bx [a.u.] \\U";DelayUpdate
			Label bottom "Transition Energy [eV]";DelayUpdate
			ModifyGraph tick=2,mirror=1,nticks=10,minor=1,fSize=14;DelayUpdate
			SetAxis bottom 280,360
			SetAxis left -0.05,0.05
			SetAxis symRot -1,1
			ModifyGraph mode=1,rgb(TDMx_Ori)=(0,0,0)
			ModifyGraph axisEnab(left)={0,0.47},axisEnab(symRot)={0.53,1},freePos(symRot)=0,lblPosMode(symRot)=1
			ModifyGraph mode=2,lsize=2
			Legend/C/N=text0/J/A=MC "\\JC"+fName+"\r\\s(TDMx_Ori) Original\r\\s(symmetryParams) Symmetry"
		endif
		
		DoWindow $gNameY
		if(!V_Flag)
			Display/N=$gNameY/K=1 TDMy_Ori vs eV_
			NewFreeAxis/L symRot
			AppendToGraph/L=symRot symmetryParams[*][1] vs eV_
			Label left "μ\\By [a.u.] \\U";DelayUpdate
			Label symRot "μ\\By [a.u.] \\U";DelayUpdate
			Label bottom "Transition Energy [eV]";DelayUpdate
			ModifyGraph tick=2,mirror=1,nticks=10,minor=1,fSize=14;DelayUpdate
			SetAxis bottom 280,360
			SetAxis left -0.05,0.05
			SetAxis symRot -1,1
			ModifyGraph mode=1,rgb(TDMy_Ori)=(0,0,0)
			ModifyGraph axisEnab(left)={0,0.47},axisEnab(symRot)={0.53,1},freePos(symRot)=0,lblPosMode(symRot)=1
			ModifyGraph mode=2,lsize=2
			Legend/C/N=text0/J/A=MC "\\JC"+fName+"\r\\s(TDMy_Ori) Original\r\\s(symmetryParams) Symmetry"
		endif
		
		DoWindow $gNameZ
		if(!V_Flag)
			Display/N=$gNameZ/K=1 TDMz_Ori vs eV_
			NewFreeAxis/L symRot
			AppendToGraph/L=symRot symmetryParams[*][2] vs eV_
			Label left "μ\\Bz [a.u.] \\U";DelayUpdate
			Label symRot "μ\\Bz [a.u.] \\U";DelayUpdate
			Label bottom "Transition Energy [eV]";DelayUpdate
			ModifyGraph tick=2,mirror=1,nticks=10,minor=1,fSize=14;DelayUpdate
			SetAxis bottom 280,360
			SetAxis left -0.05,0.05
			SetAxis symRot -1,1
			ModifyGraph mode=1,rgb(TDMz_Ori)=(0,0,0)
			ModifyGraph axisEnab(left)={0,0.47},axisEnab(symRot)={0.53,1},freePos(symRot)=0,lblPosMode(symRot)=1
			ModifyGraph mode=2,lsize=2
			Legend/C/N=text0/J/A=MC "\\JC"+fName+"\r\\s(TDMz_Ori) Original\r\\s(symmetryParams) Symmetry"
		endif
		SetDataFolder $iniFolder
	endfor
End

Function compareSymmetries(aName,fnum)
String aName
	Variable fnum
	
	String iniFolder = GetDataFolder(1)
	Variable i
	for(i=1;i<=fnum;i+=1)
		String fName = aName + num2str(i)
		SetDataFolder $fName	
		Wave originalSymmetries
		Wave symmetryParams
		Wave eV_
		String gName = fName + "_SymmetryComparison"
		DoWindow $gName
		if(!V_Flag)
			Display/N=$gName/K=1 originalSymmetries vs eV_
			NewFreeAxis/L symRot
			AppendToGraph/L=symRot symmetryParams[*][4] vs eV_
			Label left "Symmetry \\U";DelayUpdate
			Label symRot "Symmetry \\U";DelayUpdate
			Label bottom "Transition Energy [eV]";DelayUpdate
			ModifyGraph tick=2,mirror=1,nticks=10,minor=1,fSize=14;DelayUpdate
			SetAxis bottom 280,360
			SetAxis left -0.5,3.5
			SetAxis symRot -0.5,3.5
			ModifyGraph mode=1,rgb(originalSymmetries)=(0,0,0)
			ModifyGraph axisEnab(left)={0,0.47},axisEnab(symRot)={0.53,1},freePos(symRot)=0,lblPosMode(symRot)=1
			ModifyGraph mode=2,lsize=2
			ModifyGraph margin(left)=43,margin(bottom)=36
			Legend/C/N=text0/J/A=RT/B=1/X=5.00/Y=0.00 "\\JC"+fName+"\r\\s(originalSymmetries) Original\r\\s(symmetryParams) Rotated"
		endif	
		SetDataFolder $iniFolder
	endfor
End
