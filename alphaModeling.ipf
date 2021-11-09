#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//To Do: Add capability to model/fit a second orientation based on results of initial tensor fit
Function simDFT(tval,ovpVal,IPwave,mol,alpha,i0,phi,expSpecName,expEnergyName,expFolderPath,fit,rigidShift,thetaList,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,[d,holdAmps,holdWidths,holdPos,NEXAFStype,stepShift,startPre,endPre,startPost,endPost,refine,maskEnergy1,maskEnergy2,pkToRefine,refpw,holdmodTheta])
	
	Variable tval
	Variable ovpVal
	Wave IPwave
	String mol
	Variable alpha
	Variable i0
	Variable phi
	String expSpecName	//Base name of spectra
	String expEnergyName	//Base name of experimental energy waves
	String expFolderPath
	String fit	//If yes then try to fit the data, else just model it
	Variable rigidShift
	Variable stepShift
	Variable d
	String thetaList //Semicolon separated list of theta values to simulate/fit
	Variable holdAmps
	Variable holdWidths
	Variable holdPos
	String NEXAFStype
	Variable startPre 		//At what energy should the preedge start to be defined?
	Variable endPre 			//At what energy should the preedge stop to be defined?
	Variable startPost 		//At what energy should the postedge start to be defined?
	Variable endPost			//At what energy should the postedge stop to be defined?
	Variable anchorStep1,anchorStep2
	Variable anchorExp1,anchorExp2
	Variable stepWid1,stepWid2,stepE1,stepE2
	Variable refine
	Variable maskEnergy1
	Variable maskEnergy2
	Variable pkToRefine
	Wave refpw
	Variable holdmodTheta
	
	if(ParamIsDefault(holdAmps))
		holdAmps = 0
	endif
	
	if(ParamIsDefault(holdWidths))
		holdWidths = 1
	endif
	
	if(ParamIsDefault(holdPos))
		holdPos = 1
	endif
	
	if(ParamIsDefault(NEXAFStype))
		NEXAFStype = ""
	endif
	
	if(ParamIsDefault(refine))
		refine = 0
	endif
	
	if(ParamIsDefault(holdmodTheta))
		holdmodTheta = 1
	endif
	
	Variable openParams = 4 - holdAmps - holdWidths - holdPos - holdmodTheta
	Variable timerStart = startMSTimer
	
	if(!refine)	
		//Start by getting the peak parameters from the clustering algorithm
		Wave pw2d = getParams(tval,ovpVal)		
		//Combine TDM components that take place at same energy into 2D parameter wave, then make tensor for each resonance.
		makeTensor(pw2d)	
		//Normalize the tensor by dividing by the maximum TDM component of the original tensor
		normalizeTensor()		
		//Get the experimental NEXAFS and energies, make the step edges and concatenate them into a single wave
		prepSimDFTfit(expSpecName,expEnergyName,mol,IPwave,expFolderPath,anchorStep1,anchorStep2,stepWid1,stepWid2,stepE1,stepE2,stepShift = rigidShift)	
		//Make 1D parameter wave for fits
		make1DparamWave(IPwave,pw2d,alpha,i0,phi,rigidShift = rigidShift)
	endif	
	
	Wave pWave1D,paramNames,allExpSpec,allExpEnergy,allDFTSteps,allStepEnergy
	Variable n = numpnts(allExpSpec),nSpec = ItemsInList(thetaList)
	
	if(!refine)
		Duplicate/O pWave1D,pwCopy
	else
		Duplicate/O refpw,refitpWave2
	endif
	
	Wave valuesPolar = makePolarWave2(thetaList,n)	//Make wave containing theta values
	Wave pwAlpha = fitAlpha2(maskEnergy1,maskEnergy2,pWave1D,allExpEnergy,allExpSpec,nSpec,thetaList,alpha,i0,tval,ovpval,allDFTSteps)
	pwAlpha[2] = alpha
	Wave pwI0    = fitI0(maskEnergy1,maskEnergy2,pwAlpha,allExpEnergy,allExpSpec,valuesPolar,allDFTSteps,nSpec,tval,ovpval)//Fit i0
	
	if(!refine)
		Duplicate/O allExpSpec,res,results1,results2,results3,results
		results = 0;res = 0 
	endif
	
	Variable nPeaks = (numpnts(pWave1D) - 4 - pWave1D[0])/11
	Variable npnts = numpnts(allExpSpec)
	//Make Hold String for fit. Open amplitudes, hold position and width constant
	Variable pks2Fit = nPeaks//This value can be changed to troubleshoot singular matrix errors
	String H1 = makeHoldString(holdAmps,holdWidths,holdPos,pWave1D,holdModTheta,1)
	
	//Set up values for epsilon wave
	Wave eps1 = makeEpsiltonWave(pWave1D,holdAmps,holdWidths,holdPos,holdModTheta,1)

	//Make constraint wave	
	Wave/T constraint1 = makeConstraintWave(holdAmps,holdWidths,holdPos,nPeaks,pWave1D,holdmodTheta)
	
	//Do the fit or model the data using a set of DFT-BB peaks
	if(StringMatch(fit,"yes"))
		Variable V_FitTol = 0.00001,V_FitError = 0
		if(refine)
			fitSinglePeak(refpw,pkToRefine,holdPos,holdWidths,holdAmps,holdmodTheta,tval,ovpVal,alpha,mol)//adjust hold alpha
			String newpwName = "refine_pwCopy_pk" + num2str(pkToRefine)
			Wave refitpWave = $newpwName
			Abort
		else
			//Fit1 --> Open amp,position and modTheta
			Duplicate/O pwI0, pwFit
			Make/O/N=2 rcs,alphaFit
			FuncFit/Q/H=H1/M=2 simDFTfit2,pwFit, allExpSpec /X={allExpEnergy,valuesPolar,allDFTSteps} /R=res /E=eps1 /C=constraint1 /D=results1
			Wave pwFit
			NVAR redChiSq = root:redChiSq 
			redChiSq = calcRedChiSq(allExpSpec,results1)
			print "The Chi Squared is: ",V_Chisq
			print "The reduced Chi Squared is: ",redChiSq
			rcs[0] = redChiSq
			Wave M_Covar,W_Sigma,pwMolAdj
			//Fit2 --> Open everything
			Duplicate/O pwFit, pwFit2
			H1 = makeHoldString(0,0,0,pwFit2,0,0)
			Wave eps1 = makeEpsiltonWave(pwFit2,0,0,0,0,0)
			Wave/T constraint1 = makeConstraintWave(0,0,0,nPeaks,pwFit2,1)
			FuncFit/Q/H=H1/M=2 simDFTfit2,pwFit2, allExpSpec /X={allExpEnergy,valuesPolar,allDFTSteps} /R=res /E=eps1 /C=constraint1 /D=results
			NVAR redChiSq = root:redChiSq 
			redChiSq = calcRedChiSq(allExpSpec,results)
			print "The Chi Squared is: ",V_Chisq
			print "The reduced Chi Squared is: ",redChiSq
			rcs[1] = redChiSq
			alphaFit[0] = pwFit2[2]
			Wave M_Covar,W_Sigma,pwMolAdj
		endif
		Variable V_FitQuitReason
		//if(V_FitQuitReason == 0)
		//	makeCoVar(pwCopy,M_Covar,W_Sigma,paramNames,alpha)
		//endif
	else
		Make/O/N=18 alphaWave = {0,5,10,15,20,25,30,40,45,50,55,60,65,70,75,80,85,90}
		Variable i,nAlpha = numpnts(alphaWave)
		for(i=0;i<nAlpha;i+=1)
			pWave1D[2] = alphaWave[i]
			simDFTModel(pwCopy,allExpSpec,allExpEnergy,valuesPolar,allDFTSteps)
		endfor
	endif
	cleanUpTensorWaves2()
	
	//Generate the individual waves the concatenated waves
	Wave totalNEXAFS
	
	if(StringMatch(fit,"no"))
		NewDataFolder/O Modeling
		for(i=0;i<nAlpha;i+=1)
			Variable alphaVal = alphaWave[i]
			String tnxfsName = "totalNEXAFS_" + num2str(alphaVal)
			Wave totalNEXAFS = $tnxfsName
			alphaModelPlotting(fit,allExpEnergy,allExpSpec,allStepEnergy,allDFTSteps,results,totalNEXAFS,nSpec,nPeaks,pwCopy,tval,ovpVal,thetaList,d,NEXAFStype,alpha=alphaVal)		
			organizeTensorWaves(alphaVal,fit)
		endfor
	elseif(!refine)
		alphaModelPlotting(fit,allExpEnergy,allExpSpec,allStepEnergy,allDFTSteps,results,totalNEXAFS,nSpec,nPeaks,pwAlpha,tval,ovpVal,thetaList,d,NEXAFStype)		
		getAmpEnFromFit(pwFit2,pwI0,round(pwFit2[2]),tval,ovpVal,pwMolAdj)
		organizeTensorWaves(round(alpha),fit)
	else
		alphaModelPlotting(fit,allExpEnergy,allExpSpec,allStepEnergy,allDFTSteps,results,totalNEXAFS,nAlpha,nPeaks,refitpWave,tval,ovpVal,thetaList,d,NEXAFStype)		
		getAmpEnFromFit(refitpWave,pWave1D,refitpWave[2],tval,ovpVal,pwMolAdj)
		organizeTensorWaves(refitpWave[2],fit)
	endif
	Variable stopTimer = stopMSTimer(timerStart)/1000000
	print "Process ended after " + num2str(stopTimer) + " seconds."	
End

Function stepWrapper(allExpEnergy,anchorStep1,anchorStep2,mol,expSpecName,expEnergyName,IPcorr,stepWid1,stepWid2,stepE1,stepE2,nSpec)

	Wave allExpEnergy,IPcorr
	Variable anchorStep1,anchorStep2,stepWid1,stepWid2,stepE1,stepE2,nSpec
	String mol,expSpecName,expEnergyName
	
	Variable pntsPerSpec = numpnts(allExpEnergy)/nSpec
	WAVE allDFTSteps = makeStep(IPcorr,pntsPerSpec,stepWid1,stepWid2,stepE1,stepE2,mol,expSpecName,expEnergyName)
	Wave allExpSpec = scaleExptoDFTStep(allExpEnergy,allDFTSteps,anchorStep1,anchorStep2,mol,expSpecName)
	scaleDFTStepToBareAtom(allExpEnergy,allExpSpec,allDFTSteps,anchorStep1,anchorStep2,pntsPerSpec,mol)
End

Function/WAVE scaleExptoDFTStep(allExpEnergy,stepSum,anchorStep1,anchorStep2,mol,expSpecName)
	
	Wave allExpEnergy,stepSum
	Variable anchorStep1,anchorStep2
	String mol,expSpecName
	
	//Make the concatenated NEXAFS wave
	Wave allExpSpec = makeLongSpecWave(expSpecName)
	
	//Find the bare atom waves
	Wave mu_energy = findBAenergy()
	Wave mu = findBAMA(mol)
	WaveStats/Q stepSum
	
	//Find the experimental energies used for scaling   
	Variable exp280 = findWaveValAtEergy(allExpEnergy,allExpSpec,anchorStep1)
	Variable exp360 = findWaveValAtEergy(allExpEnergy,allExpSpec,anchorStep2)
	
	//Find the bare atom energies used for scaling   
	Variable bareStepLo = findWaveValAtEergy(mu_energy,mu,anchorStep1)
	Variable bareStepHi = findWaveValAtEergy(mu_energy,mu,anchorStep2)
	//		
	Variable expScale = (exp360-exp280) / (bareStepHi-bareStepLo)
	allExpSpec = allExpSpec / expScale//- V_min
	Variable exp280_2 = findWaveValAtEergy(allExpEnergy,allExpSpec,anchorStep1)
	if(exp280_2 > bareStepLo)
		allExpSpec -= V_min
	elseif(exp280 < bareStepLo)
		allExpSpec += V_min
	endif
	return allExpSpec
End

Function/WAVE scaleDFTStepToBareAtom(allExpEnergy,allExpSpec,stepSum,anchorStep1,anchorStep2,pntsPerSpec,mol)
	Wave stepSum,allExpEnergy,allExpSpec
	Variable anchorStep1,anchorStep2,pntsPerSpec
	String mol
	
	//Find the bare atom waves
	Wave mu_energy = findBAenergy()
	Wave mu = findBAMA(mol)
	
	//Find the bare atom energies used for scaling   	
	Variable bareStepLo = findWaveValAtEergy(mu_energy,mu,anchorStep1)
	Variable bareStepHi = findWaveValAtEergy(mu_energy,mu,anchorStep2)
	//Find the dft energies used for scaling 
	Variable dftStepLo = findWaveValAtEergy(allExpEnergy,stepSum,anchorStep1)//stepSum[0]
	Variable dftStepHi = findWaveValAtEergy(allExpEnergy,stepSum,anchorStep2)//stepSum[pntsPerSpec-1] 
	
	Variable scaleDFTstepToBAMA = (bareStepHi-bareStepLo) / (dftStepHi-dftStepLo)
	stepSum = (stepSum - dftStepLo) * scaleDFTstepToBAMA + bareStepLo 
	
	return stepSum
End

Function findWaveValAtEergy(xw,yw,val)
	
	Wave xw,yw
	Variable val
	
	Variable p1 = round(BinarySearchInterp(xw,val))
	Variable dftLo = yw[p1]
	
	return dftLo
End

Function/WAVE findBAenergy()
	
	String currentFolder=GetdataFolder(1)
	SetDataFolder root:Packages:NXA
	WAVE/Z mu_energy
	String newNameEn = currentFolder + "mu_energy"
	Duplicate/O mu_energy,$newNameEn
	SetDataFolder currentFolder
	
	return mu_energy
End

Function/WAVE findBAMA(mol)
	String mol
	
	String currentFolder=GetdataFolder(1)
	SetDataFolder root:Packages:NXA
	WAVE/Z mu=$(mol+"_mu")
	String newNameMu = currentFolder + "mu"
	Duplicate/O mu,$newNameMu
	
	If( !WaveExists(mu) )
		Abort "Couldn't find the indicated mass absorption wave! Aborting!"
	endif
	
	SetDataFolder currentFolder
	
	return mu
End

Function/WAVE makeLongSpecWave(expSpecName)
	
	String expSpecName
	String listOfExpSpec   = SortList(WaveList(expSpecName   + "*",";",""),";",24)
	Make/O/D/N=0 allExpSpec = 0
	Concatenate/O/NP listOfExpSpec,allExpSpec
	
	return allExpSpec
End

Function/WAVE makeStep(IPcorr,pntsPerSpec,stepWid1,stepWid2,stepE1,stepE2,mol,expSpecName,expEnergyName)

	Wave IPcorr
	Variable pntsPerSpec
	Variable stepWid1,stepWid2,stepE1,stepE2
	String mol
	String expSpecName		//Base name of spectra
	String expEnergyName	//Base name of experimental energy waves
	
	String listOfExpSpec   = SortList(WaveList(expSpecName   + "*",";",""),";",24)
	String listOfExpEnergy = SortList(WaveList(expEnergyName + "*",";",""),";",24)
	Variable nSpec = ItemsInList(listOfExpSpec)
	String dummyName = StringFromList(0,listOfExpEnergy)
	Wave dummyEnergy = $dummyName 
	
	Variable startPre= 240 		//At what energy should the preedge start to be defined?
	Variable endPre 	= 284.1		//At what energy should the preedge stop to be defined?
	Variable startPost = 284.3 		//At what energy should the postedge start to be defined?
	Variable endPost	 = 360		//At what energy should the postedge stop to be defined?
	
	//Determine min and max values for StepEdge
	WaveStats/Q IPcorr
		
	Variable stepWid, stepDecay = 0

	//Determine pre and post edges of bare atom absoprtion using a polynomial fit
	//Find the bare atom waves
	Wave mu_energy = findBAenergy()
	Wave mu = findBAMA(mol)
	FitPoly3(mu,mu_energy,startPre,endPre)
	WAVE W_Coef
	Duplicate/d/o W_Coef, PreEdge_Coef
	FitPoly3(mu,mu_energy,startPost,endPost)
	Duplicate/d/o W_Coef, PostEdge_Coef
	
	Variable StitchLow =poly(PreEdge_Coef,endPre)
	Variable StitchHigh=poly(PostEdge_Coef,startPost)
	//Make step edges for each atom
	Variable nAtoms = numpnts(IPcorr),i
	
	Variable dftPreEn2 = round(BinarySearchInterp(dummyEnergy,IPcorr[V_minloc]))
	Variable dftPosEn2 = round(BinarySearchInterp(dummyEnergy,IPcorr[V_maxloc]))
			
	Make/O/N=(pntsPerSpec) stepSum=0 
	for(i=0;i<nAtoms;i+=1)
		//Determine width for average value of IP
		Variable E0 = IPcorr[i]
			
		//stepWid =(stepWid1) +((stepWid2-stepWid1)/(stepE2-stepE1))*(E0-IPcorr[V_minloc])
		stepWid =(stepWid1) +((stepWid2-stepWid1)/(IPcorr[V_maxloc]-IPcorr[V_minloc]))*(E0-IPcorr[V_minloc])
		String stepName = "step_" + num2str(i)
		Make/o/n=(pntsPerSpec) $stepName	
		Wave w = $stepName
		w=0
		w[0,dftPreEn2-1]=poly(PreEdge_Coef,dummyEnergy)
		w[dftPreEn2,dftPosEn2]=((StitchHigh-StitchLow)*Gstep(dummyEnergy,E0,stepWid))*exp(-stepDecay*(dummyEnergy-E0-stepWid)) + StitchLow
		w[dftPosEn2+1,numpnts(w)-1]=poly(PostEdge_Coef,dummyEnergy)
		Smooth 10,w
		stepSum += w
	endfor
	
	//Get ready to make the bare atom absorption waves set in appropriate energy range and then concatenate
	String dftStepList = ""
	for(i=0;i<nSpec;i+=1)
		dftStepList = AddListItem(NameOfWave(stepSum),dftStepList)
	endfor
	
	//Make the generated DFT Step into a concatenated wave for simultaneous fitting
	Concatenate/O/NP dftStepList, allDFTSteps
	
	return allDFTSteps
End	

Function prepSimDFTfit(expSpecName,expEnergyName,mol,IPwave,expFolderPath,anchorStep1,anchorStep2,stepWid1,stepWid2,stepE1,stepE2,[startPre,endPre,startPost,endPost,stepShift])

	String expSpecName		//Base name of spectra
	String expEnergyName	//Base name of experimental energy waves
	String mol				//Name of molecule used to build the bare atom absorption
	Wave IPwave				//Wave containing the Ionization Potentials from DFT
	String expFolderPath	//Igor path to experimental data to fit to.
	Variable startPre 		//At what energy should the preedge start to be defined?
	Variable endPre 			//At what energy should the preedge stop to be defined?
	Variable startPost 		//At what energy should the postedge start to be defined?
	Variable endPost			//At what energy should the postedge stop to be defined?
	Variable stepShift		//What shift to apply to IP energies?
	Variable anchorStep1
	Variable anchorStep2
	Variable stepWid1,stepWid2,stepE1,stepE2
	
	if(ParamIsDefault(startPre))
		startPre = 240
	endif
	
	if(ParamIsDefault(endPre))
		endPre = 284.1
	endif
	
	if(ParamIsDefault(startPost))
		startPost = 284.7
	endif
	
	if(ParamIsDefault(endPost))
		endPost = 360
	endif
	
	if(ParamIsDefault(stepShift))
		stepShift = 0
	endif
	
	String currentFolder=GetdataFolder(1)
	
	Duplicate/O IPwave, IPcorr
	IPcorr = IPwave + stepShift
	
	SetDataFolder expFolderPath
	
	String listOfExpSpec   = SortList(WaveList(expSpecName   + "*",";",""),";",24)
	String listOfExpEnergy = SortList(WaveList(expEnergyName + "*",";",""),";",24)
	
	String dummyName = StringFromList(0,listOfExpEnergy)
	Wave dummyEnergy = $dummyName 
	
	Variable nSpec = ItemsInList(listOfExpSpec),i
	
	String newNameMu = currentFolder + "mu"
	String newNameEn = currentFolder + "mu_energy"
	
	SetDataFolder expFolderPath
	
	Wave mu = $newNameMu
	Wave mu_energy = $newNameEn
	
	//Concatenate the experimental angle-resolved NEXAFS and corresponding energies		
	Make/O/D/N=0 allExpSpec = 0 , allExpEnergy = 0 
	Concatenate/O/NP listOfExpSpec,allExpSpec
	Concatenate/O/NP listOfExpEnergy,allExpEnergy
	
	WaveStats/Q allExpEnergy
	
	Variable pntsPerSpec = numpnts(allExpSpec)/ItemsInList(listOfExpSpec)
	//Make the step based on the IP from DFT
	Wave allDFTSteps = makeStep(IPcorr,pntsPerSpec,stepWid1,stepWid2,stepE1,stepE2,mol,expSpecName,expEnergyName)
	//Anchor the intensity of the DFT step function to the intensity of the bare atom absorption. Means that DFT step is now in terms of MA[g/cm^2]
	Wave allDFTSteps2 = scaleDFTStepToBareAtom(allExpEnergy,allExpSpec,allDFTSteps,anchorStep1,anchorStep2,pntsPerSpec,mol)
	//Normalize the experimental NEXAFS to bare atom absorption
	Wave allExpSpec = scaleExptoDFTStep(allExpEnergy,allDFTSteps2,anchorStep1,anchorStep2,mol,expSpecName)
	//Wave allExpSpec = scaleExpToDFTStep2(allExpEnergy,allDFTSteps2,anchorStep1,anchorStep2,mol,expSpecName,nSpec,pntsPerSpec)
	//Get ready to make the bare atom absorption waves set in appropriate energy range and then concatenate
	String newNameallExpSpec = currentFolder + "allExpSpec"
	String newNameallExpEnergy = currentFolder + "allExpEnergy"
	String newNameallDFTSteps = currentFolder + "allDFTSteps"
	Duplicate/O allExpSpec,$newNameallExpSpec
	Duplicate/O allExpEnergy,$newNameallExpEnergy
	Duplicate/O allDFTSteps,$newNameallDFTSteps
	
	SetDataFolder currentFolder
End

Function make1DparamWave(IPwave,tensorPWave,alpha,i0,phi,[atomName,rigidShift])

	Wave IPwave			//1D wave containing the ionization potentials from DFT
	Wave tensorPWave		//2D wave containing peak positions,widths and amplitudes
	Variable alpha		//Molecular tilt angle
	Variable i0			//Global intensity scaling factor
	Variable phi			//Azimuthal angle of E-field. Usually 0 for NEXAFS
	String atomName		//What element are we processing?
	Variable rigidShift	//Global energy shift to be applied to DFT energies to match experiment
	
	if(ParamIsDefault(atomName))
		atomName = "C"
	endif
	
	if(ParamIsDefault(rigidShift))
		rigidShift = 0
	endif
	
	Variable nTrans = DimSize(tensorPWave,0)
	Make/O/D/N=(nTrans,3) pWave2D
	
	Variable nAtoms = numpnts(IPwave)
	
	Variable i,j=0

	for(i=0;i<=nTrans-1;i+=1)
		pWave2D[i][0] = tensorPWave[i][0]	//Peak position
		pWave2D[i][1] = tensorPWave[i][2]	//Peak width
		pWave2D[i][2] = tensorPWave[i][1]	//Peak amplitude
	endfor
	
	Make/O/D/N=(11*nTrans + nAtoms + 4) pWave1D
	
	//Make text wave with parameter names. Easier to browse parameters.
	Make/O/T/N=(11*nTrans + nAtoms + 4) paramNames 
	
	//Convert 2D wave to 1D
	for(i=nAtoms+4;i<=11*nTrans + nAtoms;i+=11)
		pWave1D[i+0] = pWave2D[j][0]	+ rigidShift//Peak position
		pWave1D[i+1] = pWave2D[j][1]	//Peak Width
		pWave1D[i+2] = pWave2D[j][2]	//Peak Max amplitude
		pWave1D[i+3] = abs(tensorPWave[j][8])	//XX
		pWave1D[i+4] = abs(tensorPWave[j][11])	//XY
		pWave1D[i+5] = abs(tensorPWave[j][12])	//XZ
		pWave1D[i+6] = abs(tensorPWave[j][9])	//YY
		pWave1D[i+7] = abs(tensorPWave[j][13])	//YZ
		pWave1D[i+8] = abs(tensorPWave[j][10])	//ZZ
		pWave1D[i+9] = alpha						//Local orientation
		pWave1D[i+10] = 0	//Modified amplitude value to be used for peak refinements
		
		paramNames[i+0]  = "peakEnergy_"     + num2str(j)
		paramNames[i+1]  = "peakWidth_"      + num2str(j)
		paramNames[i+2]  = "maxAmplitude_"   + num2str(j)
		paramNames[i+3]  = "OSxx_"           + num2str(j)
		paramNames[i+4]  = "OSxy_"           + num2str(j)
		paramNames[i+5]  = "OSxz_"           + num2str(j)
		paramNames[i+6]  = "OSyy_"           + num2str(j)
		paramNames[i+7]  = "OSyz_"           + num2str(j)
		paramNames[i+8]  = "OSzz_"           + num2str(j)
		paramNames[i+9]  = "LocalAlpha_"     + num2str(j)
		paramNames[i+10] = "MolTensorTheta_" + num2str(j)
		j+=1
	endfor
	
	//Check that there no transitions with zero energy
	
	//Add the ionization potentials to make the step edge
	for(i=4;i<nAtoms+4;i+=1)
		pWave1D[i] = IPwave[i-4]
		String ipAtom = atomName + num2str(i-3) + "_IP"
		paramNames[i] = ipAtom
	endfor
	
	//Add the number of atoms to the parameter wave
	pWave1D[0] = nAtoms
	paramNames[0] = "nAtoms"

	//Add the i0 to the parameter wave
	pWave1D[1] = i0
	paramNames[1] = "i0"
	
	//Add the alpha value to the parameter wave
	pWave1D[2] = alpha//*(Pi/180)	
	paramNames[2] = "alpha"
	
	//Add the value of phi of the Electric Field to the parameter wave
	pWave1D[3] = phi
	paramNames[3] = "phi"
End

Function makeTensor1D(pw)

	Wave pw
	
	Variable nAtoms = pw[0],i,j=0
	Variable nTransitions = (numpnts(pw) - 4 - nAtoms)/11
	for(i=0;i<nTransitions;i+=1)
		String tensorName =  "resonance_" + num2str(i)
		Make/O/D/N=(3,3) $tensorName
		Wave w = $tensorName
		w = 0 
	endfor
	
	for(i=4 + nAtoms;i<11*nTransitions + nAtoms + 4;i+=11)
		tensorName =  "resonance_" + num2str(j)
		Wave w = $tensorName
		w = 0
		w[0][0] = pw[i+3]; w[0][1] = pw[i+4]; w[0][2] = pw[i+5]
		w[1][0] = pw[i+4]; w[1][1] = pw[i+6]; w[1][2] = pw[i+7]
		w[2][0] = pw[i+5]; w[2][1] = pw[i+7]; w[2][2] = pw[i+8] 
		w = abs(w)
		j+=1
	endfor

End

Function/WAVE make3DnormWave()
		
	String normTensors = WaveList("norm_*",";","")
	Variable nTrans = ItemsInList(normTensors)
	Make/O/D/N=(3,3,nTrans) norm3D
	Concatenate/O/NP=2 normTensors, norm3D	
	
	return norm3D
End

Function/WAVE make3DResonance()
		
	String tensors = WaveList("resonance_*",";","")
	Variable nTrans = ItemsInList(tensors)
	Make/O/D/N=(3,3,nTrans) tensor3D
	Concatenate/O/NP=2 tensors, tensor3D	
	
	return tensor3D
End

Function simDFTfit2(pw,yw,ew,thw,sw) :FitFunc

	Wave pw	//1D Parameter wave. Has peak positions,widths, and amplitudes. has initial guess for alpha and IPs required to build step edge
	Wave yw	//Wave containing the different spectra to be fit
	Wave ew	//Wave containing the energies of the various spectra to be fit (i.e. the x wave)
	Wave thw	//Wave containing the various values of theta
	Wave sw	//Wave containing the concatenated steps from the DFT IPs
		
	Variable i,j=0,k,m
	Variable targetDP = 10
	//Make the absorption tensor for each transition, normalize it, and place each transition into a layer of a 3D wave
	makeTensor1D(pw)
	normalizeTensor()	
	Wave norm3D   = make3DnormWave()
	Wave tensor3D = make3DResonance()
	
	//Add the step edge
	yw = sw
	
	//Reference hold Wave. To be used to determine how the peak amplitudes will be defined
	Wave holdWave
	
	Variable nTransitions = DimSize(norm3D,2)	
	Variable pos,wid,amp,alpha,i0,phi,localalpha,amp2	
	Variable p0,pf,nSpec=howManySpec(thw),cSpec=1,thetaVal	
	
	Make/O/N=(3,3,nTransitions) filmTensor3D,r0x3D,r90x3D,r180x3D,r270x3D,tilt3D,correctedMolTensor3D
	
	Variable nAtoms = pw[0]
		
	Make/O/D/N=(3,3) currentTensor = 0 
	
	Variable nPnts = numpnts(yw)/nSpec
	
	Make/O/N=(nTransitions) MA = 0,vecAngle = 0,vecAngle2 = 0
	 
	//Wave that will populate the XX/YY and ZZ tensor elements into 2D vector
	Make/O/N=2 vector
	//Make rotation matrices for rotation of 90,180,and 270 degrees around z-axis
	Make/O/D/N=(3,3) rotMat0   = {{cos(0*(Pi/180))  ,sin(0*(Pi/180)),0}  ,{-sin(0*(Pi/180))  ,cos(0*(Pi/180)),0}  ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat90  = {{cos(90*(Pi/180)) ,sin(90*(Pi/180)),0} ,{-sin(90*(Pi/180)) ,cos(90*(Pi/180)),0} ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat180 = {{cos(180*(Pi/180)),sin(180*(Pi/180)),0},{-sin(180*(Pi/180)),cos(180*(Pi/180)),0},{0,0,1}}
	Make/O/D/N=(3,3) rotMat270 = {{cos(270*(Pi/180)),sin(270*(Pi/180)),0},{-sin(270*(Pi/180)),cos(270*(Pi/180)),0},{0,0,1}}
	
	//Duplicated parameter wave that will contain the results of rotating the molecular tensor
	Duplicate/O pw,pwMolAdj,pwFilm
	Wave pwMolAdj,pwFilm
	//*****POTENTIALLY CHANGE CLUSTERING ALGORITHM****
	//Generate clusters based on overall/total OS instead of component OS
	//How to open up individual tensor elements for fitting:
	//1. Ratio of in-plane vs out-of-plane components?
	//2. Need a second parameter to describe amplitude? Possibly related to a local alpha?
		//*********Will need to change pWave from 10 parameters per peak to 11!*******
	//3. Visualize tensor elements for each cluster (Tensor Element vs Cluster ID) -->Log plot	
	//4. Additive oscillator strength... How to accomplish this programatically???
	i0         = pw[1] 
	alpha      = pw[2] * (Pi/180)
	phi        = pw[3] * (Pi/180)
	Variable count = 0
	for(i=nAtoms + 4 ;i<11*nTransitions + nAtoms;i+=11)
		pos        = pw[i + 0]
		wid        = ABS(pw[i + 1])
		amp        = ABS(pw[i + 2])
		localalpha = pw[i + 9]
		amp2       = pw[i + 10]
		
	//	if(count == 14)
	//		Debugger
	//	endif
	//	count +=1
		
		currentTensor = tensor3D[p][q][j]// norm3D[p][q][j]
		//print atan(currentTensor[2][2]/currentTensor[0][0]) * (180/pi)
		vector[0] = currentTensor[0][0]
		vector[1] = currentTensor[2][2]
		Wave rv = vecOpsWrap(vector,amp2,"Z")
		currentTensor[0][0] =  abs(rv[0])
		currentTensor[1][1] =  abs(rv[0])
		currentTensor[2][2] =  abs(rv[1])
		
		if(WaveExists(pwMolAdj))
			pwMolAdj[i+3] = currentTensor[0][0]
			pwMolAdj[i+4] = currentTensor[0][1]
			pwMolAdj[i+5] = currentTensor[0][2]
			pwMolAdj[i+6] = currentTensor[1][1]
			pwMolAdj[i+7] = currentTensor[1][2]
			pwMolAdj[i+8] = currentTensor[2][2]
		endif
		
		Variable ang2 = atan(0.5*(tensor3D[2][2][j]/tensor3D[0][0][j])) * (180/pi)//atan(0.5*(rv[1]/rv[0])) * (180/pi)
		Variable ang3 = atan((0.5*rv[1])/rv[0]) * (180/pi)
		vecAngle[j] = ang2
		vecAngle2[j] = ang3
		currentTensor *= i0 * amp// * currentTensor
		Make/O/D/N=(3,3) rotMatAlignX = {{1,0,0},{0,cos(alpha),sin(alpha)},{0,-sin(alpha),cos(alpha)}}
			
		//Tilt the tensor by angle alpha
		MatrixOP/O  tiltedTensorX = rotMatAlignX x currentTensor x rotMatAlignX^t			 
		//Construct the film tensor by adding the azimuthal rotations 
		MatrixOP/O  RotTensor0X   =  rotMat0   x tiltedTensorX x rotMat0^t
		MatrixOP/O  RotTensor90X  =  rotMat90  x tiltedTensorX x rotMat90^t
		MatrixOP/O  RotTensor180X =  rotMat180 x tiltedTensorX x rotMat180^t
		MatrixOP/O  RotTensor270X =  rotMat270 x tiltedTensorX x rotMat270^t
		MatrixOP/O filmTensor = RotTensor0X + RotTensor90X + RotTensor180X + RotTensor270X
		
		
		Wave filmTensorT = truncate2D(filmTensor,targetDP)
		
		if(WaveExists(pwFilm))
			pwFilm[i+3] = filmTensorT[0][0]
			pwFilm[i+4] = filmTensorT[0][1]
			pwFilm[i+5] = filmTensorT[0][2]
			pwFilm[i+6] = filmTensorT[1][1]
			pwFilm[i+7] = filmTensorT[1][2]
			pwFilm[i+8] = filmTensorT[2][2]
		endif
		
		correctedMolTensor3D[][][j] = currentTensor[p][q]	//Molecular tensor after rotating ZZ/XX elements by ModTheta
		filmTensor3D[][][j]         = filmTensorT[p][q]		//Film Tensor after tilting molecular tensor by alpha and 4-fold addition
		r0x3D[][][j]                = RotTensor0X[p][q]
		r90x3D[][][j]               = RotTensor90X[p][q]
		r180x3D[][][j]              = RotTensor180X[p][q]
		r270x3D[][][j]              = RotTensor270X[p][q]
		tilt3D[][][j]               = tiltedTensorX[p][q]
		
		if(j<nTransitions)
			if(cSpec == 1)
				p0 = 0
				pf = nPnts-1
				thetaVal = thw[p0] * (Pi/180)
			elseif(cSpec == nSpec)
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0] * (Pi/180)
			else
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0] * (Pi/180)
			endif
		endif
		//Need to extract xx and zz tensor elements and create NEXAFS waves from them
		
		//Apply electric field to current tensor to obtain the mass absorbance
		Make/O/D/N=(3) eField = {sin(thetaVal)*cos(phi),sin(thetaVal)*sin(phi),cos(thetaVal)}	
		MatrixOP/O tempMA = eField^t x filmTensor x eField
		
		if(numtype(tempMA[0]) !=0)
			print "Transition " + num2str(j) + " with a theta value of " + num2str(thetaVal) + " has a problem"
			//MA[j] = 0
		else
			MA[j] = tempMA[0] //abs(tempMA[0])
		endif
		//Make the BB peaks for display
		String pkName = "pk" + num2str(j) + "_spec" + num2str(cSpec)
		Make/O/N=2000 $pkName
		Wave z = $pkName
		ControlInfo/W=ClusteringAlgorithm Eini
		Variable Eini = V_Value
			
		ControlInfo/W=ClusteringAlgorithm Efin
		Variable Efin = V_Value
			
		SetScale/i x,Eini,Efin,z
		z = 0
		z = (MA[j]) * gauss(x,pos,wid)

		//yw[p0,pf] += sqrt(2*Pi) * wid * MA[j] * gauss(ew,pos,wid)
		yw[p0,pf] += MA[j] * gauss(ew,pos,wid)
		j+=1
		
		if(j<=nTransitions)
			if(cSpec == 1)
				if(j==(nTransitions) && cSpec < nSpec)
					i=nAtoms + 4 - 11
					cSpec+=1
					j=0
				endif 			
			else
				if(j==(nTransitions) && cSpec < nSpec)
					i=nAtoms + 4 - 11
					cSpec+=1
					j=0
				endif 
			endif
		endif
	endfor
	
	KillWaves/Z rotMat0,rotMat90,rotMat180,rotMat270,RotTensor0,RotTensor90,RotTensor180,RotTensor270,TEMPma
	
End

Function simDFTModel(pw,yw,ew,thw,sw) :FitFunc

	Wave pw	//1D Parameter wave. Has peak positions,widths, and amplitudes. has initial guess for alpha and IPs required to build step edge
	Wave yw	//Wave containing the different spectra to be fit
	Wave ew	//Wave containing the energies of the various spectra to be fit (i.e. the x wave)
	Wave thw	//Wave containing the various values of theta
	Wave sw	//Step Wave
		
	Variable i,j=0,k,m
	
	//Make the absorption tensor for each transition, normalize it, and place each transition into a layer of a 3D wave
	makeTensor1D(pw)
	normalizeTensor()	
	make3DnormWave()
	Wave norm3D
	
	Variable nTransitions = DimSize(norm3D,2)
	Variable pos,wid,amp,alpha,i0,phi,localalpha,amp2	
	Variable p0,pf,nSpec=howManySpec(thw),cSpec=1,thetaVal	
	Variable nAtoms = pw[0]
	Variable targetDP = 10

	Make/O/D/N=(3,3) currentTensor = 0
	//This is the wave that will contain the DFT tensor model NEXAFS
	String tNXFSname = "totalNEXAFS_" + num2str(pw[2]) 
	Duplicate/D/O yw,$tNXFSname
	Wave totalNEXAFS = $tNXFSname 
	totalNEXAFS = sw
	Variable nPnts = numpnts(yw)/nSpec
	Make/O/D/N=(nTransitions) MA = 0 
	
	//Make rotation matrices for rotation of 90,180,and 270 degrees around z-axis
	Make/O/D/N=(3,3) rotMat0   = {{cos(0*(Pi/180))  ,sin(0*(Pi/180)),0}  ,{-sin(0*(Pi/180))  ,cos(0*(Pi/180)),0}  ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat90  = {{cos(90*(Pi/180)) ,sin(90*(Pi/180)),0} ,{-sin(90*(Pi/180)) ,cos(90*(Pi/180)),0} ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat180 = {{cos(180*(Pi/180)),sin(180*(Pi/180)),0},{-sin(180*(Pi/180)),cos(180*(Pi/180)),0},{0,0,1}}
	Make/O/D/N=(3,3) rotMat270 = {{cos(270*(Pi/180)),sin(270*(Pi/180)),0},{-sin(270*(Pi/180)),cos(270*(Pi/180)),0},{0,0,1}}
	
	for(i=nAtoms + 4 ;i<11*nTransitions + nAtoms;i+=11)
		i0    = pw[1]
		alpha = pw[2]//90-pw[i+9]
		phi   = pw[3]
		pos   = pw[i + 0]
		wid   = pw[i + 1]
		amp   = pw[i + 2]
		localalpha = pw[i + 9]
		amp2       = pw[i + 10]
		
		currentTensor = i0 * amp * norm3D[p][q][j]
		
		Make/O/D/N=(3,3) rotMatAlignX = {{1,0,0},{0,cos(alpha*(Pi/180)),sin(alpha*(Pi/180))},{0,-sin(alpha*(Pi/180)),cos(alpha*(Pi/180))}}
		
		MatrixOP/O tiltedTensor = rotMatAlignX x currentTensor x rotMatAlignX^t
		
		MatrixOP/O  RotTensor0   =  rotMat0   x tiltedTensor x rotMat0^t
		MatrixOP/O  RotTensor90  =  rotMat90  x tiltedTensor x rotMat90^t
		MatrixOP/O  RotTensor180 =  rotMat180 x tiltedTensor x rotMat180^t
		MatrixOP/O  RotTensor270 =  rotMat270 x tiltedTensor x rotMat270^t
		
		MatrixOP/O filmTensor = RotTensor0 + RotTensor90 + RotTensor180 + RotTensor270
		
		for(k=0;k<=2;k+=1)
			for(m=0;m<=2;m+=1)
				Variable currComp = filmTensor[k][m]
				targetDP = round(targetDP)
				currComp = round(currComp * (10^targetDP)) / (10^targetDP)
				if(currComp == -0)
					currComp = 0
				endif
				filmTensor[k][m] = currComp
			endfor
		endfor
		
		if(j<nTransitions)
			if(cSpec == 1)
				p0 = 0
				pf = nPnts-1
				thetaVal = thw[p0]
			elseif(cSpec == nSpec)
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0]
			else
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0]
			endif
		endif
		
		Make/O/D/N=(3) eField = {sin(thetaVal*(Pi/180))*cos(phi*(Pi/180)),sin(thetaVal*(Pi/180))*sin(phi*(Pi/180)),cos(thetaVal*(Pi/180))}	
		MatrixOP/O tempMA = amp * (eField^t x filmTensor x eField)
		
		if(numtype(tempMA[0]) !=0)
			MA[j] = 0
		else
			MA[j] = tempMA[0]
		endif
		
		String pkName = "pk" + num2str(j) + "_spec" + num2str(cSpec) + "_Alpha" + num2str(pw[2]) 
		Make/O/N=2000 $pkName
		Wave z = $pkName
		SetScale/i x,280,360,z
		z =  MA[j] * gauss(x,pos,wid)		
		totalNEXAFS[p0,pf] += MA[j] * gauss(ew,pos,wid)
		j+=1
		
		if(j<=nTransitions)
			if(cSpec == 1)
				if(j==(nTransitions) && cSpec < nSpec)
					i=nAtoms + 4 - 11
					cSpec+=1
					j=0
				endif 			
			else
				if(j==(nTransitions) && cSpec < nSpec)
					i=nAtoms + 4 - 11
					cSpec+=1
					j=0
				endif 
			endif
		endif		
	endfor
	
End

Function modifyWaves(w,nSpec,thetaList)

	Wave w
	Variable nSpec
	String thetaList
	
	Variable i,j=0,n = DimSize(w,0),pps=n/nSpec,newSpec=ItemsInList(thetaList),newSize = newSpec*pps
	Duplicate/O w, test
	Redimension/N=(newSize,3) test
	//Copy step values and energies into redimensioned wave
	for(i=n;i<newSize;i+=1)
		test[i][0] = w[j][0]
		test[i][2] = w[j][2]
		j+=1
		if(j>=pps)
			j=0
		endif
	endfor
	//Fill the redimensioned wave with new theta values
	j=0
	for(i=0;i<newSize;i+=1)
		Variable cTheta = str2num(StringFromList(j,thetaList))
		test[i][1] = cTheta
		if((mod(i+1,pps)==0) && i!=0)
			j+=1
		endif
	endfor
End

Function modelWrapper(thetaList,pw,eval,pwOri,d,[fit,fetchWaves,moveStepY])
	Wave pw,pwOri
	Variable eVal,fit,fetchWaves,moveStepY,d
	String thetaList
	
	//Look for the waves used in the fit and to construct the ETS wave?
	if(fetchWaves)
		Wave holdWave     = findHoldWave() 
		Wave valuesPolar  = findThetaWave() 
		Wave allExpEnergy = findExpEnergyWave() 
		Wave allDFTSteps  = findStepWave() 
		Wave allExpSpec   = findExpSpecWave() 
		Wave eps          = findEpsWave()
	else
		Wave holdWave,valuesPolar,allExpEnergy,allDFTSteps,allExpSpec,eps
	endif
	Duplicate/O allDFTSteps,DFTStepNew
	DFTStepNew = allDFTSteps + moveStepY 
	Duplicate/O allExpEnergy,totalNEXAFS
	Duplicate/O pw,pwMolAdj,pwFilm,pwRes
	totalNEXAFS = 0
	simDFTfit2(pw,totalNEXAFS,allExpEnergy,valuesPolar,DFTStepNew)
	Wave totalNEXAFS,pwMolAdj
	Variable nSpec = ItemsInList(thetaList)
	splitSpec(totalNEXAFS,nSpec,"NXFS")
	splitSpec(DFTStepNew,nSpec,"dftStep")
	//plotAlphaIntensities2(eVal,enw,"NXFS",nSpec,thetaList)
	
	//Make function that plots DFT Model against Experimental NEXAFS
	if(!fit)
		plotModel(thetaList)	//Plot Change in parameters 
		getAmpEnFromFit(pw,pwOri,20,1,50,pwMolAdj,d=d)
		//gaugeChange(pwOri,pw)	
	else 						   //Fit the model
		DFTfitNew(pw,holdWave,valuesPolar,allExpEnergy,allDFTSteps,allExpSpec,eps)
	endif
End

Function/WAVE findHoldWave()

	String iniFolder = GetDataFolder(1)
	String miscFolder = iniFolder + "TensorMisc"
	SetDataFolder $miscFolder
	Wave holdWave
	SetDataFolder $iniFolder
	String nName = iniFolder + "holdWave"
	Duplicate/O holdWave,$nName
	return holdWave
End

Function/WAVE findThetaWave()

	String iniFolder = GetDataFolder(1)
	String miscFolder = iniFolder + "TensorMisc"
	SetDataFolder $miscFolder
	Wave valuesPolar
	SetDataFolder $iniFolder
	String nName = iniFolder + "valuesPolar"
	Duplicate/O valuesPolar,$nName
	return valuesPolar
End

Function/WAVE findExpEnergyWave()

	String iniFolder = GetDataFolder(1)
	String miscFolder = iniFolder + "TensorMisc"
	SetDataFolder $miscFolder
	Wave allExpEnergy
	SetDataFolder $iniFolder
	String nName = iniFolder + "allExpEnergy"
	Duplicate/O allExpEnergy,$nName
	
	return allExpEnergy
End

Function/WAVE findStepWave()

	String iniFolder = GetDataFolder(1)
	String miscFolder = iniFolder + "TensorMisc"
	SetDataFolder $miscFolder
	Wave allDFTSteps
	SetDataFolder $iniFolder
	String nName = iniFolder + "allDFTSteps"
	Duplicate/O allDFTSteps,$nName
	return allDFTSteps
End

Function/WAVE findExpSpecWave()

	String iniFolder = GetDataFolder(1)
	String miscFolder = iniFolder + "TensorMisc"
	SetDataFolder $miscFolder
	Wave allExpSpec
	SetDataFolder $iniFolder
	String nName = iniFolder + "allExpSpec"
	Duplicate/O allExpSpec,$nName
	return allExpSpec
End

Function/WAVE findEpsWave()

	String iniFolder = GetDataFolder(1)
	String miscFolder = iniFolder + "TensorMisc"
	SetDataFolder $miscFolder
	Wave eps
	SetDataFolder $iniFolder
	String nName = iniFolder + "eps"
	Duplicate/O eps,$nName
	return eps
End

Function plotModel(thetaList,[fit])

	String thetaList
	Variable fit
	
	Variable n = ItemsInList(thetaList),i,pDiff
	Make/O/N=(n) percentDiff
	String iniFolder = GetDataFolder(1)
	String DFTNXFS = WaveList("NXFS*",";",""),fitNXFS = WaveList("fitResult*",";","")
	String pks = WaveList("pk*"+num2str(2),";","")
	Variable npks = ItemsInList(pks)
	String fitFolder  = iniFolder + "TensorMisc"
	SetDataFolder $fitFolder
	String energies = WaveList("expEnergy*",";",""),enWaveName = StringFromList(0,energies)
	Wave enWave = $enWaveName 
	String expNXFS = WaveList("expSpec*",";","")	
	String steps = WaveList("dftStep*",";","")
	Wave step = $StringFromList(0,steps)
	DoWindow ModelvsExpPlot
	
	String Rlist = "0;1;39321;65535",Glist = "3204;34817;1;43690",Blist = "13107;52428;1;0"
	DoWindow ModelvsExpPlot
	if(!V_Flag)
		Display/N=ModelvsExpPlot/K=1/W=(0,0,1500,400)
		NewFreeAxis/L residuals
		NewFreeAxis/L peaks
		
		for(i=0;i<n;i+=1)
			Variable r = str2num(StringFromList(i,Rlist))	,g = str2num(StringFromList(i,Glist)),b = str2num(StringFromList(i,Blist))	
			SetDataFolder $fitFolder
			String expName = StringFromList(i,expNXFS)
			Wave nxfs = $expName
			AppendToGraph/W=ModelvsExpPlot nxfs vs enWave
			
			Variable x = numpnts(nxfs)
			String resName = "res" + num2str(i)
			Make/O/N=(x) $resName
			Wave res = $resName
			
			SetDataFolder $iniFolder
			if(!fit)
				String dftName = StringFromList(i,DFTNXFS)
				Wave dft  = $dftName
				AppendToGraph/W=ModelvsExpPlot dft vs enWave
				ModifyGraph lstyle($dftName)=3
				res = ((nxfs - dft)/nxfs)*100
				pDiff = calcPercentDiff(nxfs,dft)
				ModifyGraph rgb($expName)=(r,g,b),rgb($dftName)=(r,g,b),lsize($expName)=3,lsize($dftName)=3
			else
				String fitName = StringFromList(i,fitNXFS)
				Wave fitW  = $fitName
				AppendToGraph/W=ModelvsExpPlot fitW vs enWave
				ModifyGraph lstyle($fitName)=3
				res =(( nxfs - fitW)/nxfs)*100	
				pDiff = calcPercentDiff(nxfs,fitW)
				ModifyGraph rgb($expName)=(r,g,b),rgb($fitName)=(r,g,b),lsize($expName)=3,lsize($fitName)=3
			endif
			
			AppendToGraph/W=ModelvsExpPlot/L=residuals res vs enWave
			ModifyGraph rgb($resName)=(r,g,b)
			percentDiff[i] = pDiff
		endfor
		
		SetDataFolder $iniFolder
		
		//Append the DFT step edge onto graph
		AppendToGraph/W=ModelvsExpPlot step vs enWave
		ModifyGraph lsize($StringFromList(0,steps))=2,rgb($StringFromList(0,steps))=(0,0,0)
		
		for(i=0;i<npks;i+=1) 	
			String pkName = StringFromList(i,pks)
			Wave pk  = $pkName
			AppendToGraph/W=ModelvsExpPlot/L=peaks pk
			ModifyGraph rgb($pkName)=(0,0,0),lsize($pkName)=1.5
		endfor
		
		ModifyGraph grid=2,mirror=1,minor=1,fStyle=1,fSize=12,axisEnab(left)={0.2,0.8}
		ModifyGraph lblPosMode(residuals)=1,axisEnab(residuals)={0.82,1},freePos(residuals)=0
		ModifyGraph lblPosMode(peaks)=1,axisEnab(peaks)={0,0.18},freePos(peaks)=0,lblPosMode(left)=1
		Label left "Mass Absorbance [cm\S2\M/g]\\U"
		Label peaks "Peaks [a.u.]\\U"
		Label residuals "Residuals[%]\\U"
		Label bottom "Transition Energy[eV]"
		SetAxis bottom 283,310
		SetAxis residuals -100,100
		SetAxis peaks *,70000
		//Make the legend
		String totalLegend = "θ[°]  EXP  DFT  %Diff\r",legendPortion = "",whichList
		if(fit)
			whichList = fitNXFS
		else
			whichList = DFTNXFS
		endif
		for(i=0;i<n;i+=1)
			String alpha = StringFromList(i,thetaList)
			String dftSpec = StringFromList(i,whichList)
			String expSpec = StringFromList(i,expNXFS)
			if(i != (n-1))
				legendPortion += alpha + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) +"\r\n"
			else
				legendPortion += alpha + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) 
			endif
		endfor
		Legend/C/A=RT/N=text0/J totalLegend + legendPortion
	else
		for(i=0;i<n;i+=1)	
			SetDataFolder $fitFolder
			expName = StringFromList(i,expNXFS)
			Wave nxfs = $expName
			
			x = numpnts(nxfs)
			resName = "res" + num2str(i)
			Make/O/N=(x) $resName
			Wave res = $resName
			
			SetDataFolder $iniFolder
			if(!fit)
				dftName = StringFromList(i,DFTNXFS)
				Wave dft  = $dftName
				res = ((nxfs - dft)/nxfs)*100
				pDiff = calcPercentDiff(nxfs,dft)
			else
				fitName = StringFromList(i,fitNXFS)
				Wave fitW  = $fitName
				res =(( nxfs - fitW)/nxfs)*100	
				pDiff = calcPercentDiff(nxfs,fitW)
			endif
			percentDiff[i] = pDiff
		endfor
		
		totalLegend = "θ[°]  EXP  DFT  %Diff\r"
		legendPortion = ""
		if(fit)
			whichList = fitNXFS
		else
			whichList = DFTNXFS
		endif
		for(i=0;i<n;i+=1)
			alpha = StringFromList(i,thetaList)
			dftSpec = StringFromList(i,whichList)
			expSpec = StringFromList(i,expNXFS)
			if(i != (n-1))
				legendPortion += alpha + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) +"\r\n"
			else
				legendPortion += alpha + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) 
			endif
		endfor
		Legend/C/A=RT/N=text0/J/W=ModelvsExpPlot totalLegend + legendPortion	
	endif
	SetDataFolder $iniFolder
End

Function gaugeChange(pwOri,pwNew)

	Wave pwOri,pwNew
	Variable nPeaks = (numpnts(pwOri)-pwOri[0]-4)/11,i,j=0
	Make/O/N=(nPeaks) difPos,difWid,difAmp,difAmp2
	for(i=pwOri[0]+4;i<numpnts(pwOri);i+=11)
		Variable oriPos   = pwOri[i+0]
		Variable newPos   = pwNew[i+0]
		Variable oriWid   = pwOri[i+1]
		Variable newWid   = pwNew[i+1]
		Variable oriAmp   = pwOri[i+2]
		Variable newAmp   = pwNew[i+2]
		Variable oriAmp2  = pwOri[i+10]
		Variable newAmp2  = pwNew[i+10]
		//Use log to track change
		//log(new/old)
		difPos[j]  = (newPos - oriPos)///oriPos)*100
		difWid[j]  = ((newWid - oriWid)/oriWid)*100
		difAmp[j]  = ((newAmp - oriAmp)/oriAmp)*100
		difAmp2[j] = ((newAmp2 - oriAmp2)/oriAmp2)*100
		j+=1
	endfor
	
	DoWindow ParameterChanges
	if(!V_Flag)
		NewFreeAxis/L ampchange
		NewFreeAxis/L modampchange
		Display/N=ParameterChanges/K=1 difPos
		AppendToGraph/L=ampchange difAmp
		AppendToGraph/L=modampchange difAmp2
		ModifyGraph grid=2,mirror=1,nticks=10,fStyle=1,fSize=12
		Label left "Position Diff [eV]"
		Label ampchange "Amplitude Diff"
		Label modampchange "Ratio Diff"
		Label bottom "Peak ID"
		ModifyGraph mode=4,marker=19,lsize=1.5
		ModifyGraph lblPosMode(left)=1,lblPosMode(ampchange)=1,lblPosMode(modampchange)=1
		ModifyGraph axisEnab(left)={0,0.3},axisEnab(ampchange)={0.34,0.64},axisEnab(modampchange)={0.68,1}
		ModifyGraph freePos(ampchange)=0,freePos(modampchange)=0
	endif
End

Function DFTfitNew(pw,hw,tw,xw,sw,ew,epsw)

	Wave pw,hw,tw,xw,sw,ew,epsw //parameter,hold,theta,energy,step,experimental,and epsilon waves
	
	//Wave makeSingleConstraintWave(holdAmps,holdWidths,holdPos,holdAlpha,holdSecAmp,pWave1D,pkID)
	String H = holdWaveToStr(hw)
	Duplicate/O ew,fitRes,results
	Duplicate/O pw,refineFitpw
	Variable V_FitError = 0.000001
	FuncFit/H=H/M=2/Q simDFTfit2,refineFitpw, ew /X={xw,tw,sw} /R=fitRes /E=epsw /D=results// /C=constraint 
	Wave totalNEXAFS 
	splitSpec(totalNEXAFS,4,"fitResult")
	plotModel("40;55;70;90",fit=1)
End

Function/S holdWaveToStr(hw)

	Wave hw
	
	Variable n = numpnts(hw),i
	String H =""
	for(i=0;i<n;i+=1)
		H += num2str(hw[i])
	endfor
//	print H
	return H
End

Function cleanUpTensorWaves2()

	String resonanceList  = WaveList("resonance_*",";","")	
	String normResonanceList  = WaveList("norm_resonance_*",";","")
	String pWaveList = SortList(WaveList("pWave_*",";",""))
	
	Variable nResonances  = ItemsInList(resonanceList)
	Variable nPwaves = ItemsInList(pWaveList)
	
	Variable k
	
	//Place the tensors for each step into a 3D wave. Each layer corresponds to a resonance tensor
	Make/O/N=(3,3,nResonances) oritensor3D,normTensor3D
	
	Concatenate/NP=2/O resonanceList,oritensor3D 	
	Concatenate/NP=2/O normResonanceList,normTensor3D
	
	//Remove the 2D tensor waves that clutter the workspace
	for(k=0;k<=nResonances-1;k+=1)
		String currentResWave  = StringFromList(k,resonanceList)
		String currentNormResWave  = StringFromList(k,normResonanceList)
		
		Wave w = $currentResWave
		Wave z = $currentNormResWave
		KillWaves/Z $currentResWave,$currentNormResWave
	endfor
	
	for(k=0;k<=nPwaves-1;k+=1)
		String currentPwave   = StringFromList(k,pWaveList)
		Wave b = $currentPwave
		KillWaves/Z $currentPwave
	endfor
	
	KillWaves/Z rotMatAlignX,rotMatAlignY,rotMatAlignZ,rotXY,rotYX
	KillWaves/Z totAzi,tensor3D,TransitionInstances
End

Function/WAVE makePolarWave2(thetaList,nPnts)

	String thetaList
	Variable nPnts
	
	Make/O/D/N=(nPnts) valuesPolar
	Variable nThetas = ItemsInList(thetaList)
	Variable pointsPerTheta = nPnts/nThetas
	
	Variable i,j=0
	for(i=0;i<nPnts;i+=1)
		valuesPolar[i] =str2num(StringFromList(j,thetaList))
		if(i==(j+1)*pointsPerTheta-1)
			j+=1
		endif
	endfor
	
	return valuesPolar
End

Function splitSpec(w,n,name,[alpha,reScale])

	Wave w	//Wave containing concatenated spectra
	Variable n	//How many spectra are there in the concatenated wave?
	String name	//What should be the base names of the split spectra?
	Variable alpha	//What's the alpha value for the generated spectra?
	Variable reScale	//SetScale of wave?
	
	Variable pnts = numpnts(w)/n
	String cSpec
	Variable i,j=0,k=0
	
	for(i=1;i<=n;i+=1)
		if(ParamIsDefault(alpha))
			cSpec = name + num2str(i)// + "_alpha" + num2str(alpha)
		else
			cSpec = name + num2str(i) + "_alpha" + num2str(alpha)
		endif
		Make/O/N=(pnts) $cSpec
		Wave y = $cSpec
		if(reScale)
			SetScale/i x,280,360,y	//Update this so that the x scale is obtained from the panel
		endif
		
		for(j=j;j<=i*pnts-1;j+=1)
			
			if(n==1)
				y[k] = w[j]
			else
				y[k] = w[j]
			endif
			k+=1
		endfor
		k=0
	endfor 
End

Function plotPeaks(nSpec,nPeaks,alpha,tval,ovpVal)

	Variable nSpec,nPeaks,alpha,tval,ovpVal
	
	Variable i,j
	for(j=1;j<=nSpec;j+=1)
		String graphName = "Spec_" + num2str(j) +"_Alpha" + replaceString(".",num2str(alpha),"p") + "_OS" + replaceString(".",num2str(tval),"p") + "_OVP" + replaceString(".",num2str(ovpVal),"p")
		DoWindow $graphName
		if(!V_Flag)
			Display/N=$graphName
			for(i=0;i<nPeaks;i+=1)
				String cPeakName = "pk" + num2str(i) + "_spec" + num2str(j)  
				Wave z = $cPeakName
				AppendToGraph/W=$graphName z
			endfor
			Label left "Mass Absorbance [cm2/g] \\U";DelayUpdate
			Label bottom "Transition Energy[eV]";DelayUpdate
			ModifyGraph grid=2,mirror=1,nticks=10,minor=1,fStyle=1,lsize=2
			ApplyColorTableToTopGraph("ColdWarm")
		endif
	endfor
End

Function plotResults(baseNameList,alpha,eName,thetaList,tval,ovpVal)
	
	String baseNameList
	Variable alpha
	String eName
	String thetaList
	Variable tval
	Variable ovpVal
	
	Variable nItems = ItemsInList(baseNameList,";")

	Variable i,j,k

	for(i=0;i<nItems;i+=1)
	
		String currentWaveList = StringFromList(i,baseNameList,";")
		String currentSetList = WaveList(currentWaveList+"*alpha"+num2str(alpha),";","")
		String energyList = WaveList(eName +"*alpha"+num2str(alpha),";","")
		Variable nSet = ItemsInList(currentSetList)
		
		String plotName = currentWaveList +"_Alpha" + replaceString(".",num2str(alpha),"p") + "_OS" + replaceString(".",num2str(tval),"p") + "_OVP" + replaceString(".",num2str(ovpVal),"p")
		DoWindow $plotName
			
		if(!V_Flag)			
			String totalLegend = ""
		
			for(j=0;j<nSet;j+=1)
				String currentWave = StringFromList(j,currentSetList)
				String currentEnergy = StringFromList(j,energyList)
				Wave w = $currentWave
				Wave x = $currentEnergy
				
				if(j==0)
				Display/N=$plotName/K=1 w vs x
				elseif(j<nSet)
					AppendToGraph/W=$plotName w vs x
				endif
				
				String currentTheta = StringFromList(j,thetaList)
				String legendPortion 	
				if(StringMatch(currentWaveList,"dftSpec"))
				if(j<nSet-1)
					legendPortion = "\\s(" + currentWave + ") " + currentTheta + "\r" 
					else
						legendPortion = "\\s(" + currentWave + ") " + currentTheta
					endif
				elseif(StringMatch(currentWaveList,"dftStep"))
					if(j<nSet-1)
						legendPortion = "\\s(" + currentWave + ") " + currentTheta + "\r"
					else
						legendPortion = "\\s(" + currentWave + ") " + currentTheta
					endif
				elseif(StringMatch(currentWaveList,"expSpec"))		
					if(j<nSet-1)
						legendPortion = "\\s(" + currentWave + ") " + currentTheta + "\r"
					else
						legendPortion = "\\s(" + currentWave + ") " + currentTheta
					endif
				endif
		
				totalLegend += legendPortion 	
			endfor
			
			if(StringMatch(currentWaveList,"dftSpec"))	
				Legend/C/N=text0/J "DFT \r \\JCθ[°]\r \\JL" + totalLegend	
			elseif(StringMatch(currentWaveList,"dftStep"))
				Legend/C/N=text0/J "Step \r \\JCθ[°]\r \\JL" + totalLegend	
			elseif(StringMatch(currentWaveList,"expSpec"))
				Legend/C/N=text0/J "Experiment \r\\JCθ[°]\r \\JL" + totalLegend	
			endif
			
			Label left "Mass Absorbance [cm\\S2\\M/g] \\U";DelayUpdate
			Label bottom "Transition Energy [eV]";DelayUpdate
			ModifyGraph grid=2,mirror=1,nticks=10,minor=1,fStyle=1,fSize=16,lsize=1.5
			SetAxis bottom 280,320
			ApplyColorTableToTopGraph("ColdWarm")
		endif
	endfor
	
End

Function plotResultsFit2(baseNameList,alpha,eName,thetaList,tval,ovpMax,nPeaks,NEXAFStype)
	
	String baseNameList
	Variable alpha
	String eName
	String thetaList
	Variable tval
	Variable ovpMax
	Variable nPeaks
	String NEXAFStype
	
	Variable nItems = ItemsInList(baseNameList,";")

	Variable i,j
	
	alpha = round(alpha)
	for(i=0;i<nItems;i+=1)
		String currentWaveList = StringFromList(i,baseNameList,";")
		String currentSetList = WaveList(currentWaveList+"*alpha"+num2str(alpha),";","")
		String energyList = WaveList(eName +"*alpha"+num2str(alpha),";","")
		Variable nSet = ItemsInList(currentSetList)
	endfor

	String summaryPlot = "Final_Comparison" + replaceString(".",num2str(alpha),"p") + "_OS" + replaceString(".",num2str(tval),"p") + "_OVP" + replaceString(".",num2str(ovpMax),"p") + "_" + NEXAFStype 
	DoWindow $summaryPlot	
	if(!V_Flag)
		Display/N=$summaryPlot/K=1
		
		for(i=0;i<nSet;i+=1)
			String currentExp  = StringFromList(0,baseNameList,";")
			String currentDFT  = StringFromList(1,baseNameList,";")
			
			String currentExpList  = WaveList(currentExp  +"*alpha"+num2str(alpha),";","")
			String currentDFTList  = WaveList(currentDFT  +"*alpha"+num2str(alpha),";","")
			String enList          = WaveList(eName +"*alpha"+num2str(alpha),";","")
		
			String currentTheta = StringFromList(i,thetaList)	
			
			String currentExpWave  = StringFromList(i,currentExpList ,";")
			String currentDFTWave  = StringFromList(i,currentDFTList ,";")
			String currentEnerWave = StringFromList(i,enList,";")
		
			Wave w = $currentExpWave
			Wave y = $currentDFTWave		
			Wave x = $currentEnerWave
			
			AppendToGraph w,y vs x
			ModifyGraph lstyle($currentDFTWave)=3
		endfor
		
		ApplyColorTableToTopGraph("ColdWarm")
		
		//Append Step
		String currentStepWave = "dftStep1_alpha" + num2str(alpha)
		Wave z = $currentStepWave
		AppendToGraph z vs x
		ModifyGraph rgb($currentStepWave)=(39321,39321,39321)
		
		Variable magicAngleSpec = findMagicAngle(thetaList)
		NewFreeAxis/L BBPeaks					
		//Append BB peaks to plot
		for(i=0;i<nPeaks;i+=1)
			String cPeakName = "pk" + num2str(i) + "_spec" + num2str(magicAngleSpec)  
			Wave pk = $cPeakName
			AppendToGraph/L=BBPeaks/W=$summaryPlot pk
			ModifyGraph rgb($cPeakName)=(0,0,0)
		endfor
		
		Label left "Mass Absorbance [cm\\S2\\M/g] \\U";DelayUpdate
		Label bottom "Transition Energy [eV]";DelayUpdate
		ModifyGraph grid=2,mirror=1,nticks=20,minor=1,fStyle=1,fSize=16,lsize=1.5
		ModifyGraph axisEnab(left)={0.35,1},axisEnab(BBPeaks)={0,0.3},freePos(BBPeaks)=0,lblPosMode=1
		SetAxis bottom 283,320
		SetAxis BBPeaks 0,100000//*
		Label BBPeaks "BB Peaks\\u"
		String legendText = "\\s(" + currentExpWave +") Experiment \\s(" + currentDFTWave + ") DFT\r\\s(" + currentStepWave +") Step"
		Legend/C/N=text0/J "\\JCα = " + num2str(alpha) + "°\r" + legendText
	endif
	
End

Function/WAVE getParams(tval,ovpVal)

	Variable tval,ovpVal
	
	String iniFolder = GetDataFolder(1)	
	String dftFitFolder = "All"
	//Look for folder containing results from amplitude fitting DFT clusters of each component
	String iniPwaveName
	Variable fit = 0
	if(DataFolderExists(dftFitFolder))
		SetDataFolder	$dftFitFolder
		if(fit)
			iniPwaveName = "combClusterPWAll"//Use this if using the pwave before fitting to dft
		else	
			iniPwaveName = "pw2dOriginal"//Use this using pwave after fitting to dft
		endif
		Wave pWave = $iniPwaveName
		SetDataFolder	iniFolder
		Duplicate/O pWave, $iniPwaveName
	else	
		print "Data folder for " + dftFitFolder + " not found."
	endif
	Wave pw2D = $iniPwaveName
	
	return pw2d
End

Function/WAVE make2Dfrom1DPwave(pw1d)

	Wave pw1d
	
	Variable i=0,j,nClusters = (numpnts(pw1d)-4-pw1d[0])/11
	Make/O/N=(nClusters,3) pw2D
	for(j=4+pw1d[0];j<4+pw1d[0]+11*nClusters;j+=11)
		pw2D[i][0]  = pw1D[j]	//Energy
		pw2D[i][1]  = pw1D[j+1]//Width
		pw2D[i][2]  = pw1D[j+2]//Amp
		i+=1
	endfor	
	
	return pw2d
End

//This should be at the end of the clustering process.
Function makeTensor(pWave)

	Wave pWave
	
	Variable nClusters = DimSize(pWave,0),i=0
	Variable tol = 1/100	//What is the minimum threshold of values to consider for the TDM tensor? If component is less than the V_max*tol then make it 0
	Variable targetDP = 10
	
	for(i=0;i<nClusters;i+=1)
		String tensorName =  "resonance_" + num2str(i)
		Make/O/N=(3,3) $tensorName
		Wave w = $tensorName
		w = 0
		w[0][0] = pWave[i][8] ; w[0][1] = pWave[i][11]; w[0][2] = pWave[i][12]
		w[1][0] = pWave[i][11]; w[1][1] = pWave[i][9] ; w[1][2] = pWave[i][13]
		w[2][0] = pWave[i][12]; w[2][1] = pWave[i][13]; w[2][2] = pWave[i][10] 
		truncateSym(w,targetDP,tol)
	endfor	
End

Function normalizeTensor()

	String resonanceList = WaveList("resonance_*",";","")
	Variable nResonances = ItemsInList(resonanceList)
	Variable k
	
	//Add a wave that will contain the max amplitudes
	Make/O/N=(nResonances) maxAmplitudes
		
	for(k=0;k<=nResonances-1;k+=1)
		String currentTensorName = StringFromList(k,resonanceList)
		Wave currentTensor  = $currentTensorName
		WaveStats/Q $currentTensorName
		String normTensorName = "norm_resonance_" + num2str(k)
		Duplicate/O currentTensor,$normTensorName
		Wave normTensor = $normTensorName
		//normTensor  = normTensor/V_max
		maxAmplitudes[k] = V_max
	endfor
End

Function getAmpEnFromFit(pwFit,pwOri,alpha,tval,ovpVal,pwMolAdj,[d])

	Wave pwFit
	Wave pwOri
	Wave pwMolAdj
	Variable alpha
	Variable tval
	Variable ovpVal
	Variable d
	
	Variable nPeaks = (numpnts(pwOri) - pwOri[0] - 4)/11
	Variable i,j=0,k=0
	Variable i0dIF = pwFit[1]
	Make/o/n=(nPeaks) iniAmps,fitAmps,iniEns,fitEns,iniTheta,fitTheta
	Make/O/N=(nPeaks) ampChange,enChange,modThetaChange,tdmThetaChange,widChange	
	Make/O/N=(nPeaks) iniTDMTheta,fitTDMTheta,iniWidth,fitWidth
	//Get the fitted and original Peak Amplitudes,Widths, and Positions
	for(i=pwOri[0] + 4;i<11*nPeaks+pwOri[0]+4;i+=11)
		fitEns[j]    = pwFit[i]
		fitAmps[j]   = pwFit[i+2]
		fitTheta[j]  = pwFit[i+10]
		fitTDMTheta[j] = atan(pwMolAdj[i+8]/pwMolAdj[i+6])*(180/pi)
		fitWidth[j] = pwFit[i+1]*2.355
		
		iniEns[j]   = pwOri[i]
		iniAmps[j]  = pwOri[i+2]
		iniTheta[j] = pwOri[i+10]
		iniTDMTheta[j] = atan(pwOri[i+8]/pwOri[i+6])*(180/pi)
		iniWidth[j] = pwOri[i+1]*2.355
		Variable wid = pwFit[i+1]*2.355
		enChange[j]    = (fitEns[j]   - iniEns[j])/(wid)	//Relative energy change
		modThetaChange[j] = fitTheta[j] - iniTheta[j]
		ampChange[j]   = log(fitAmps[j]/iniAmps[j])//Include i0 into iniAmps to fix this ratio!
		tdmThetaChange[j] =  fitTDMTheta[j] - iniTDMTheta[j] 
		widChange[j] = fitWidth[j] - iniWidth[j]
		j+=1
	endfor
	
	if(d)
		String gName = "FitParamChange_OS" +num2str(tval) + "_OVP" + num2str(ovpval)
		DoWindow $gName
		if(!V_Flag)
			Display/N=$gName/W=(0,0,300,500)/K=1 ampChange
			NewFreeAxis/L deltaEn 
			NewFreeAxis/L deltaTh
			NewFreeAxis/L deltaTDM
			AppendToGraph/W=$gName/L=deltaEn enChange
			AppendToGraph/W=$gName/L=deltaTh modThetaChange
			AppendToGraph/W=$gName/L=deltaTDM fitTDMTheta,iniTDMTheta//tdmThetaChange
			ModifyGraph mirror=1,minor=1,fStyle=1,fSize=14,lblPosMode(left)=1,lblPosMode(deltaEn)=1,lblPosMode(deltaTh)=1,lblPosMode(deltaTDM)=1
			ModifyGraph axisEnab(left)={0,0.22},axisEnab(deltaEn)={0.25,0.47},axisEnab(deltaTh)={0.5,0.72},axisEnab(deltaTDM)={0.75,1},freePos(deltaEn)=0,freePos(deltaTh)=0,freePos(deltaTDM)=0
			ModifyGraph mode=4,marker=19,rgb(enChange)=(1,34817,52428),rgb(modThetaChange)=(1,39321,19939)//,rgb(tdmThetaChange)=(19729,1,39321)
			Label left "log(A\\BF\\M/A\\BI\\M)"
			Label deltaEn "∆E/FWHM"
			Label deltaTh "∆θ[°]"
			Label deltaTDM "Final,Initial[°]"
			Label bottom "Cluster ID"
			SetAxis deltaEn -1,1
			SetAxis deltaTh -90,90
			ModifyGraph manTick(deltaTh)={0,45,0,0},manMinor(deltaTh)={2,0}
			SetAxis deltaTDM 0,90
			ModifyGraph manTick(deltaTDM)={0,45,0,0},manMinor(deltaTDM)={2,0}
			SetAxis left -1,1
			ModifyGraph zero(left)=1,zero(deltaEn)=1,zero(deltaTh)=1,zero(deltaTDM)=1,zeroThick(left)=2,zeroThick(deltaEn)=2,zeroThick(deltaTh)=2,zeroThick(deltaTDM)=2
			ModifyGraph rgb(fitTDMTheta)=(19729,1,39321),rgb(iniTDMTheta)=(0,0,0)
		endif
	endif
//	String gName2 = "FitAmplitudes_Alpha_" + replaceString(".",num2str(alpha),"p") + "OS_" + replaceString(".",num2str(tval),"p") + "OVP_" + replaceString(".",num2str(ovpVal),"p")
//	DoWindow $gName2
//	if(!V_Flag)
//		Display/N=$gName2/K=1 fitAmps,iniAmplitudes VS fitEnergies
//		Label left "Transition Amplitude [a.u.] \U"
//		Label bottom "Peak Energy[eV]"
//		ModifyGraph grid=2,zero(left)=1,mirror=1,nticks(bottom)=nPeaks,nTicks(left)=10,minor=1,fStyle=1,zeroThick(left)=2
//		ModifyGraph lsize=1.5
//		ModifyGraph marker=16,mode(fitAmps)=8,toMode(fitAmps)=1,mode(iniAmplitudes)=3
//		ModifyGraph rgb(iniAmplitudes)=(0,0,65535)
//		ModifyGraph useNegRGB(fitAmps)=1,usePlusRGB(fitAmps)=1,plusRGB(fitAmps)=(65535,43690,0),negRGB(fitAmps)=(29524,1,58982)
//		Legend/C/N=text0/J/A=MC "\\JCAmplitudes\r\\s(fitAmps) Fit\r\\s(iniAmplitudes) Orignal"
//	endif
End

//This function fits the molecular tilt angle to a given energy range defined by E1-E2. The pi-manifold is a good choice.
Function/WAVE fitAlpha(E1,E2,pw,xw,ew,tw,sw,nSpec,os,ovp)

	Wave pw,xw,ew,sw,tw	//Parameter wave, concatenated energy wave,concatenated experimental NEXAFS,concatenated sample theta wave
	Variable E1,E2	//Energy range that we want to fit alpha in
	Variable nSpec,os,ovp//How many spectra are we fitting?
	
	Wave mw = makeMaskWave(E1,E2,xw)
	String H = makeHoldStrAlphaFit(pw,E2)
	Wave eps = makeEpsiltonWaveAlphaFit(pw,E2)
	Wave con = makeConstraintAlphaFit(pw,E2)
	Duplicate/O ew,alphaFitResults,res
	Duplicate/O pw, pwAlphaFit
	Variable V_FitError=0
	FuncFit/H=H simDFTfit2,pwAlphaFit, ew /X={xw,tw,sw} /R=res /E=eps /D=alphaFitResults /M=mw /C=con
	splitSpec(alphaFitResults,nSpec,"alphaFit",alpha=pwAlphaFit[2])
	Variable alpha = pwAlphaFit[2],i0 = pwAlphaFit[1]
//	plotFitResults("alphaFit",alpha,i0,os,ovp)
	
	return pwAlphaFit
End

Function/WAVE fitAlpha2(E1,E2,pw,xw,ew,nSpec,tl,alpha,i0,os,ovp,sw)
	Wave pw,xw,ew,sw	//Parameter wave, concatenated energy wave,concatenated experimental NEXAFS,concatenated sample theta wave
	Variable E1,E2	//Energy range that we want to fit alpha in
	Variable nSpec //How many spectra are we fitting?
	Variable alpha,i0//Tilt angle guess
	Variable os,ovp
	String tl
	Duplicate/O pw, pwAlphaFit
	Variable nEnergies = detEnPnts(xw,E1,E2),i
	Make/O/N=(nEnergies,(nSpec+1)) alphaFitVals
	Make/O/N=(nEnergies) energies,alphas
	Duplicate/O ew,ew_Sub
	ew_Sub = ew - sw
	Variable newAlpha = popAlphaPWave(xw,E1,E2,ew_Sub,alphaFitVals,nSpec,tl,alpha,i0,os,ovp)
	pwAlphaFit[2] = newAlpha 
	return pwAlphaFit
End

Function detEnPnts(w,E1,E2)
	Wave w
	Variable E1,E2
	
	Variable p1 = round(BinarySearchInterp(w,E1))
	Variable p2 = round(BinarySearchInterp(w,E2))
	Variable dif = p2-p1
	//print p1,p2,dif
	return dif
End

Function popAlphaPWave(xw,E1,E2,ew,pw,nSpec,tl,alpha,i0,os,ovp)
	Wave xw//Energy wave
	Wave ew//Experiment wave
	Wave pw//2D wave that will be populated with energies and intensities to fit
	String tl//List containing sample theta values
	Variable E1,E2,nSpec,alpha,i0,os,ovp
	
	Variable i,j,pps=numpnts(xw)/nSpec//,n = DimSize(pw,0)
	//Variable p1 = round(BinarySearchInterp(xw,E1))
	//Variable p2 = round(BinarySearchInterp(xw,E2))
	FindLevel/Q/R=[0,pps-1] xw,E1
	Variable p1 = round(V_levelX)
	FindLevel/Q/R=[0,pps-1] xw,E2
	Variable p2 = round(V_levelX)
	Variable n = p2-p1
	Make/O/N=(nSpec) thetas
	Make/O/N=(n) alphas,energies,sigmas,avgAlpha
	//This loop will populate the energies for desired energy range
	for(i=0;i<n;i+=1)
		pw[i][0] = xw[p1 + i]
		energies[i] = xw[p1 + i]		
	endfor
	
	//This loop will populate the intensities for desired energy range for each sample theta into pw
	for(j=0;j<nSpec;j+=1)
		for(i=0;i<n;i+=1)
			pw[i][j] = ew[p1 + j*pps + i]
		endfor
	endfor
	
	//This loop will populate the theta wave based on the theta list
	for(i=0;i<nSpec;i+=1)
		thetas[i] = str2num(StringFromList(i,tl))
	endfor
	
	//This loop will make and populate the 1d waves containing:
	//1. the intensites for each sample theta for each energy
	//2. The pw for the alpha fits
	//3. The result waves of the fits
	for(i=0;i<n;i+=1)
		String name = "aFit_" + num2str(i)
		String name2 = "aval_" + num2str(i)
		String name3 = "alphaResultFit_" + num2str(i)
		String name4 = "alphaEpsFit_" + num2str(i)
		String name5 = "alphaResiFit_" + num2str(i)
		String name6 = "alphaConsFit_" + num2str(i)
		Make/O/N=(nSpec) $name,$name3,$name5
		Make/D/O/N=(2) $name2,$name4
		Make/O/T/N=1 $name6
		Wave w = $name
		Wave fitpw = $name2
		Wave w3 = $name3
		Wave eps = $name4
		eps[0] = 1e-4
		eps[1] = 1e-6
		Wave res = $name5
		Wave/T con = $name6
		con[0] = {"K0>=0","K0<=90","K1>=0"}
		fitpw[0] = alpha
		fitpw[1] = i0
		for(j=0;j<nSpec;j+=1)
			w[j] = pw[i][j]
		endfor
		
		//Fit alpha
		FuncFit/Q fitStohrAlpha , fitpw, w /X=thetas /D=w3 /R=res /E=eps /C=con
		Wave W_Sigma
		alphas[i] = fitpw[0]
		sigmas[i] = W_sigma[0]
		KillWaves $name,$name3,$name2,$name4,$name5,$name6
	endfor
	
	WaveStats/Q alphas
	avgAlpha = V_avg
	Variable alphaMin = V_min
	WaveStats/Q sigmas
	Variable avgError = V_avg
	print avgAlpha
	//Display the results of the alpha fit
//	String gName = "alphaFit_OS" + num2str(os) + "_OVP" + num2str(ovp)
//	DoWindow $gName
//	if(!V_Flag)
//		Display/N=$gName/K=1 alphas,avgAlpha vs energies
//		ErrorBars/W=$gName alphas Y,wave=(sigmas,sigmas)
//		ModifyGraph mirror=1,minor=1,fStyle=1,fSize=14
//		ModifyGraph mode(alphas)=4,marker(alphas)=19,lstyle(avgAlpha)=3,lsize(avgAlpha)=2,rgb(avgAlpha)=(0,0,0)
//		Label left "α[°]"
//		Label bottom "Transition Energy[eV]"
//		Legend/C/N=text0/J/A=MC "\\s(alphas) alphas\r\\s(avgAlpha) avgAlpha\rAVG α = "+num2str(avgAlpha[0])+" ± "+num2str(avgError)+"°"
//		Legend/C/N=text0/J/A=MC "\\s(alphas) alphas\r\\s(avgAlpha) avgAlpha\rAVG α = "+num2str(alphaMin)+" ± "+num2str(avgError)+"°"
//	endif
	
	return alphaMin//avgAlpha[0]//
End

Function fitStohrAlpha(pw,yw,xw):FitFunc
	Wave pw //This parameter wave contains the value of alpha to be fitted
	Wave yw//This is the intensity wave at a given energy
	Wave xw//This is the theta wave
	pw[0] *= (pi/180)
	xw *= (pi/180)
	Variable i,n=numpnts(xw)
	//if(pw[0] > 55)
	//	Debugger
	//endif
	
	for(i=0;i<n;i+=1)
		if(pw[0]*(180/pi) < 55)
			yw[i] = pw[1]*(1/3)*(1+0.50*(3*cos(pw[0])^2-1)*(3*cos(xw[i])^2-1))//Vector Orbital
		else
			yw[i] = pw[1]*(2/3)*(1-0.25*(3*cos(pw[0])^2-1)*(3*cos(xw[i])^2-1))//Planar Orbital
		endif
	endfor
	pw[0] *= (180/pi)//Convert alpha back to degrees
	xw *= (180/pi)
End

Function fitStohrAlpha2(pw,yw,xw)
	Wave pw //This parameter wave contains the value of alpha to be fitted
	Wave yw//This is the intensity wave at a given energy
	Wave xw//This is the theta wave
	pw[0] *= (pi/180)
	xw *= (pi/180)
	Variable i,n=numpnts(xw)
	Duplicate/O yw,test
	for(i=0;i<n;i+=1)
		test[i] = pw[1]*(1/3)*(1+0.5*(3*cos(pw[0])^2-1)*(3*cos(xw[i])^2-1))
	endfor
	pw[0] *= (180/pi)
	xw *= (180/pi)
End

Function plotFitResults(fitStage,alpha,i0,os,ovp)
	String fitStage
	Variable alpha
	Variable i0,os,ovp
	
	Variable i,j=0
	
	Wave M_colors = makeColorWave(0.5)
	Variable nColors = DimSize(M_colors,0)
	
	String wList = WaveList(fitStage+"*alpha"+num2str(round(alpha)),";","")
	String xList = WaveList("energy*",";","")
	String eList = WaveList("exp*",";","")
	Variable n = itemsInList(wList)
	String pltName = fitStage + "_Results" + "OS" +num2str(os) + "_OVP" + num2str(ovp)
	DoWindow $pltName
	if(!V_Flag)
		Display/N=$pltName/K=1 
		for(i=0;i<n;i+=1)
			String cwn = StringFromList(i,wList)
			String xwn = StringFromList(i,xList)
			String ewn = StringFromList(i,eList)
			Wave w = $cwn, e = $ewn ,x = $xwn
			AppendToGraph/W=$pltName w,e vs x
			ModifyGraph lstyle($cwn)=3,lsize($ewn)=1.5
			ModifyGraph rgb($cwn)=(M_colors[j][0],M_colors[j][1],M_colors[j][2],M_colors[j][3])
			ModifyGraph rgb($ewn)=(M_colors[j][0],M_colors[j][1],M_colors[j][2],M_colors[j][3])
			j+=1
			if(j>=nColors)
				j = 0
			endif
		endfor
		Label left "Mass Absorbance [cm\S2\M/g] \U"
		Label bottom "Transition Energy [eV]"
		ModifyGraph mirror=1,nticks=10,minor=1,fStyle=1,fSize=14
		if(StringMatch(fitStage,"alphaFit"))
			TextBox/C/N=text0/A=MC "α="+num2str(alpha)+"°"
		elseif(StringMatch(fitStage,"i0Fit"))
			TextBox/C/N=text0/A=MC "i0="+num2str(i0)
		endif
	endif
End

Function/WAVE fitI0(E1,E2,pw,xw,ew,tw,sw,nSpec,tval,ovpval)

	Wave pw,xw,ew,sw,tw	//Parameter wave, concatenated energy wave,concatenated experimental NEXAFS,concatenated sample theta wave
	Variable E1,E2	//Energy range that we want to fit alpha in
	Variable nSpec //How many spectra are we fitting?
	Variable tval,ovpval
	Wave mw = makeMaskWave(E1,E2,xw)
	String H = makeHoldStrI0Fit(pw)
	Wave eps = makeEpsiltonWaveI0Fit(pw)
	Duplicate/O ew,i0FitResults,res
	Duplicate/O pw, pwi0Fit
	Variable V_FitError=0
	FuncFit/Q/H=H simDFTfit2,pwi0Fit, ew /X={xw,tw,sw} /R=res /E=eps /D=i0FitResults/M=mw
	splitSpec(i0FitResults,nSpec,"i0Fit",alpha=round(pwi0Fit[2]))
	splitSpec(xw,nSpec,"energy")
	splitSpec(ew,nSpec,"exp")
	Variable alpha = pwi0Fit[2],i0 = pwi0Fit[1]
	//plotFitResults("i0Fit",alpha,i0,tval,ovpval)
	return pwi0Fit
End

Function organizeTensorWaves(alpha,fit)
	
	Variable alpha
	String fit

	String iniFolder = GetDataFolder(1)
	Variable i,j
	
	NewDataFolder/O TensorMisc
	String misc = WaveList("!totalNEXAFS*",";","")
	Variable nMisc = ItemsInList(misc)
	
	for(j=0;j<nMisc;j+=1)
		String currentWave = StringFromList(j,misc)
		Wave w = $currentWave
		String destFolder = iniFolder + "TensorMisc:"
		MoveWave w,$destFolder
	endfor	
	
	
	String finalFolderName = "Alpha_" + replaceString(".",num2str(alpha),"p")
	NewDataFolder/O $finalFolderName
	
	MoveDataFolder/O=3 TensorMisc,$finalFolderName
	
	if(StringMatch(fit,"no"))
		MoveDataFolder/O=3 $finalFolderName,Modeling
	endif
End

Function howManySpec(w)
	
	Wave w
	
	Variable dim = WaveDims(w)
	Variable n,i,j=1,val
	if(dim == 1)
		n = numpnts(w)
	else
		n = DimSize(w,0)
	endif
	
	for(i=1;i<n;i+=1)
		if(dim == 1)
			val = w[i]
			if(val != w[i-1])
				j+=1
			endif
		else
			val = w[i][1]
			if(val != w[i-1][1])
				j+=1
			endif
		endif
	endfor

	return j

End

Function/S makeThetaList(w)
	
	Wave w
	
	Variable dim = WaveDims(w)
	Variable n,i,val
	String thetaList = ""
	if(dim == 1)
		n = numpnts(w)
	else
		n = DimSize(w,0)
	endif
	
	for(i=1;i<n;i+=1)
		if(dim == 1)
			val = w[i]
			if((i-1)==0)
				thetaList = AddListItem(num2str(val),thetaList)
			elseif(val != w[i-1])
				thetaList = AddListItem(num2str(val),thetaList)
			endif
		else
			val = w[i][1]
			if((i-1)==0)
				thetaList = AddListItem(num2str(val),thetaList)
			elseif(val != w[i-1][1])
				thetaList = AddListItem(num2str(val),thetaList)
			endif
		endif
	endfor
	thetaList = SortList(thetaList,";",2)
	return thetaList
End

Function makeCoVar(pw,cvMat,sd,pNames,alpha)

	Wave pw,cvMat,sd,pNames
	Variable alpha
	
	Variable i,j=0
	Variable n=numpnts(pw)
	Variable nPeaks = (numpnts(pw) - 4 - pw[0])/11
	
	Duplicate/O cvMat,cvMadAdj,corrMat
	Duplicate/O pw,pw2
	Duplicate/O pNames,pNames2
	
	for(i=0;i<n;i+=1)
		Variable sdev = sd[i]
		Variable corr = cvMadAdj[i][j]/sdev
		if(numtype(corr) != 0 )
			corrMat[i][j] = 0
		else
			corrMat[i][j] = corr
		endif
		
		if(j>=n-1)
			break
		endif
		
		if(i==n-1)
			i=0
			j+=1
		endif
		
	endfor
	
	//Remove Normalized Values of Oscillator Strength. Not used for fitting.
	for(i=11*nPeaks+4+pw[0]-6;i>4 + pw[0];i-=11)
		DeletePoints i,6,cvMadAdj,corrMat
		DeletePoints/M=1 i,6,cvMadAdj,corrMat
		DeletePoints i,6,pw2,pNames2
	endfor
	
	//Remove Ionization Potentials. Not used for fitting.
	for(i=pw[0] + 3;i>3;i-=1)
		DeletePoints	 i,1,cvMadAdj,corrMat
		DeletePoints/M=1	 i,1,cvMadAdj,corrMat
		DeletePoints i,1,pw2,pNames2
	endfor
	
	////Remove phi and number of Atoms. Not used for fitting.
	
	String corrGName = "CORRELATION_" + replaceString(".",num2str(alpha),"p")
	DoWindow $corrGName//CORRELATION
	if(!V_Flag)
		WaveStats/Q corrMat
		NewImage/N=$corrGName/K=1 corrMat
		ModifyImage corrMat ctab= {*,*,RedWhiteBlue256,0}
		ModifyGraph nticks=10,fStyle=1,gridStyle=2,gridRGB=(0,0,0),margin(left)=25,margin(top)=25,margin(right)=90,fsize=10
		Label left "Fit Parameters";DelayUpdate
		Label top "Fit Parameters"
		ColorScale/C/N=text0/A=RC/E  ctab={-0.2,0.2,RedWhiteBlue256,0},axisRange={V_min,V_max}
		ColorScale/C/N=text0 "Correlation \\U"
	endif
End

Function FitPoly3(yw,xw,StartX, EndX)
	WAVE yw, xw
	Variable startX, EndX
	FindLevel/P/Q xw, startX
	Variable startP=floor(V_LevelX)
	Findlevel/P/Q xw, endX
	Variable endP=floor(V_LevelX)
	CurveFit/W=2/Q poly 3, yw[startP,endP]/X=xw[startP,endP]
end

Function Lstep(x, x0,width)
	Variable x, x0, width
	return 1/2 + 1/pi * atan((x-x0)/(width/2))
end

Function Gstep(x,x0,width)
	Variable x, x0, width
	Variable c=2*sqrt(2)
	return 1/2 + 1/2*erf((x-x0)/(width/c))
end

Function/S makeHoldString(holdAmps,holdWidths,holdPos,pWave1D,holdModTheta,holdAlpha)
	Variable holdAmps,holdWidths,holdPos,holdModTheta,holdAlpha
	Wave pWave1D
	
	Variable i,nPeaks = (numpnts(pwave1d)-pwave1d[0]-4)/11
	//Make Hold String for fit. Open amplitudes, hold position and width constant
	if(holdAlpha)
		String H = "1111"		//H1 = Number of atoms, i0, alpha, phi
	else
		H = "1101"
	endif
	for(i=1;i<=pWave1D[0];i+=1)
		H +="1" //For Ionization Potentials from DFT. Used to build the Step Edge. Hold them constant.
	endfor
	
	for(i=0;i<=nPeaks-1;i+=1)
		//Hold string for peak parameters. Position, Width, Max Amplitude, xxNoRM, xyNorm,xzNorm,yyNorm,yzNorm,zzNorm
		//Hold cases:
		if(!holdPos && !holdWidths && !holdAmps && !holdModTheta)   
			H +="00011111110"//"00011111110"
		elseif(holdPos && !holdWidths && holdAmps && !holdModTheta)  
			H +="10111111110"
		elseif(holdPos && !holdWidths && !holdAmps && !holdModTheta)  
			H +="10011111110"//"10011111110"
		elseif(holdPos && holdWidths && !holdAmps && !holdModTheta)  
			H +="11011111110"//"11011111110"
		elseif(!holdPos && holdWidths && !holdAmps && !holdModTheta) 
			H +="01011111110"
		elseif(!holdPos && !holdWidths && holdAmps && !holdModTheta)
			H +="00111111110"
		elseif(!holdPos && holdWidths && holdAmps && !holdModTheta)
			H +="01111111110"
		elseif(holdAmps && holdWidths && holdAmps && !holdModTheta) 
			H +="11111111110"
		elseif(!holdPos && !holdWidths && !holdAmps && holdModTheta)   
			H +="00011111111"//"00011111110"
		elseif(holdPos && !holdWidths && holdAmps && holdModTheta)  
			H +="10111111111"
		elseif(holdPos && !holdWidths && !holdAmps && holdModTheta)  
			H +="10011111111"//"10011111110"
		elseif(holdPos && holdWidths && !holdAmps && holdModTheta)  
			H +="11011111111"//"11011111110"
		elseif(!holdPos && holdWidths && !holdAmps && holdModTheta) 
			H +="01011111111"
		elseif(!holdPos && !holdWidths && holdAmps && holdModTheta)
			H +="00111111111"
		elseif(!holdPos && holdWidths && holdAmps && holdModTheta)
			H +="01111111111"
		elseif(holdAmps && holdWidths && holdAmps && holdModTheta) 
			H +="11111111111"	
		endif
	endfor
	    
	Make/O/N=(numpnts(pWave1D)) holdWave
	for(i=0;i<numpnts(pWave1D);i+=1)
		holdWave[i] = str2num(H[i,i])
	endfor
	
	//print H
	return H
End

Function/S makeHoldStrAlphaFit(pw,Emax)
	Wave pw
	Variable Emax
	
	Variable i,nPeaks = (numpnts(pw)-pw[0]-4)/11
	//Make Hold String for fit. Open amplitudes, hold position and width constant
	String H	= "1101"	//H1 = Number of atoms, i0, alpha, phi	
	
	for(i=1;i<=pw[0];i+=1)
		H +="1" //For Ionization Potentials from DFT. Used to build the Step Edge. Hold them constant.
	endfor
	
	for(i=0;i<=nPeaks-1;i+=1)  
		Variable E = pw[i*11+4+pw[0]]
		if( E < Emax)
			H +="11011111111"//If peak position is less than max mask energy then open amplitudes
		else
			H +="11111111111"//If peak position is greater than max mask energy then close amplitudes
		endif
	endfor
	
	Make/O/N=(numpnts(pw)) holdWave_alphaFit
	for(i=0;i<numpnts(pw);i+=1)
		holdWave_alphaFit[i] = str2num(H[i,i])
	endfor
	return H
End

Function/S makeHoldStrI0Fit(pw)
	Wave pw
	
	Variable i,nPeaks = (numpnts(pw)-pw[0]-4)/11
	//Make Hold String for fit. Open amplitudes, hold position and width constant
	String H	= "1011"	//Number of atoms, i0, alpha, phi. Open only i0
	
	for(i=1;i<=pw[0];i+=1)
		H +="1" //For Ionization Potentials from DFT. Used to build the Step Edge. Hold them constant.
	endfor
	
	for(i=0;i<=nPeaks-1;i+=1)  
		H +="11111111111"
	endfor
	
	Make/O/N=(numpnts(pw)) holdWave_I0Fit
	for(i=0;i<numpnts(pw);i+=1)
		holdWave_I0Fit[i] = str2num(H[i,i])
	endfor
	return H
End

Function/WAVE makeConstraintAlphaFit(pw,Emax)
	Wave pw
	Variable Emax	
	
	Variable i,nPeaks = (numpnts(pw)-pw[0]-4)/11,pos,j=0
	String c
	Make/O/T/N=2 constraintAlpha	
	constraintAlpha[0] = {"K2 >= 0","K2 <= 90"}
	for(i=4 + pw[0];i < 4 + pw[0] + 11*(nPeaks);i+=11)
		pos   = pw[i+0]	
		if(pos <= Emax)
			c = "K" + num2str(i+2)  + " >= 0"
			Redimension/N=(2+j) constraintAlpha
			constraintAlpha[j+1]  = c
			j+=1
		else
			break
		endif
	endfor
	return constraintAlpha
End

Function/WAVE makeEpsiltonWave(pWave1D,holdAmps,holdWidths,holdPos,holdModTheta,holdAlpha)
	
	Wave pWave1D
	Variable holdAmps,holdWidths,holdPos,holdModTheta,holdAlpha
	
	Variable i,nPeaks = (numpnts(pwave1d)-pwave1d[0]-4)/11
	//Set up values for epsilon wave
	Make/O/N=(numpnts(pWave1D)) eps
	
	Variable epsPos  = 1e-3
	Variable epsWid  = 1e-6
	Variable epsAmp  = 1e-6
	Variable epsMT   = 1e-2
	eps[0] = 0
	eps[1] = 0
	if(holdAlpha)
		eps[2] = 0	
	else
		eps[2] = 1e-2	
	endif
	eps[3] = 0
	
	for(i=4;i<pWave1D[0] + 4;i+=1)
		eps[i]   = 0	//For Ionization Potentials from DFT. Used to build the Step Edge
	endfor
	
	for(i=4 + pWave1D[0];i<11*nPeaks+4+pWave1D[0];i+=11)
		
		if(!holdPos && !holdWidths && !holdAmps && holdModTheta)   
			eps[i + 0]    = epsPos		
			eps[i + 1]    = epsWid
			eps[i + 2]    = epsAmp
			eps[i + 10]   = 0
		elseif(holdPos && !holdWidths && holdAmps && holdModTheta)  
			eps[i + 0]   = 0	
			eps[i + 1]   = epsWid
			eps[i + 2]   = 0
			eps[i + 10]   = 0
		elseif(holdPos && !holdWidths && !holdAmps && holdModTheta)  
			eps[i + 0]    = 0		
			eps[i + 1]    = epsWid
			eps[i + 2]    = epsAmp
			eps[i + 10]   = 0
		elseif(holdPos && holdWidths && !holdAmps && holdModTheta)  
			eps[i + 0]    = 0	
			eps[i + 1]    = 0
			eps[i + 2]    = epsAmp
			eps[i + 10]   = 0
		elseif(!holdPos && holdWidths && !holdAmps && holdModTheta) 
			eps[i + 0]    = epsPos		
			eps[i + 1]    = 0
			eps[i + 2]    = epsAmp
			eps[i + 10]   = 0
		elseif(!holdPos && !holdWidths && holdAmps && holdModTheta)
			eps[i + 0]   = epsPos		
			eps[i + 1]   = epsWid
			eps[i + 2]   = 0
			eps[i + 10]   = 0
		elseif(!holdPos && holdWidths && holdAmps && holdModTheta)
			eps[i + 0]   = epsPos		
			eps[i + 1]   = 0
			eps[i + 2]   = 0
			eps[i + 10]   = 0
		elseif(holdAmps && holdWidths && holdAmps && holdModTheta) 
			eps[i + 0]   = 0		
			eps[i + 1]   = 0
			eps[i + 2]   = 0
			eps[i + 10]   = 0
		elseif(!holdPos && !holdWidths && !holdAmps && !holdModTheta)   //
			eps[i + 0]    = epsPos		
			eps[i + 1]    = epsWid
			eps[i + 2]    = epsAmp
			eps[i + 10]    = epsMT
		elseif(holdPos && !holdWidths && holdAmps && !holdModTheta)  
			eps[i + 0]   = 0	
			eps[i + 1]   = epsWid
			eps[i + 2]   = 0
			eps[i + 10]   = epsMT
		elseif(holdPos && !holdWidths && !holdAmps && !holdModTheta)  
			eps[i + 0]    = 0		
			eps[i + 1]    = epsWid
			eps[i + 2]    = epsAmp
			eps[i + 10]    = epsMT
		elseif(holdPos && holdWidths && !holdAmps && !holdModTheta)  
			eps[i + 0]    = 0	
			eps[i + 1]    = 0
			eps[i + 2]    = epsAmp
			eps[i + 10]    = epsMT
		elseif(!holdPos && holdWidths && !holdAmps && !holdModTheta) 
			eps[i + 0]    = epsPos		
			eps[i + 1]    = 0
			eps[i + 2]    = epsAmp
			eps[i + 10]    = epsMT
		elseif(!holdPos && !holdWidths && holdAmps && !holdModTheta)
			eps[i + 0]   = epsPos		
			eps[i + 1]   = epsWid
			eps[i + 2]   = 0
			eps[i + 10]   = epsMT
		elseif(!holdPos && holdWidths && holdAmps && !holdModTheta)
			eps[i + 0]   = epsPos		
			eps[i + 1]   = 0
			eps[i + 2]   = 0
			eps[i + 10]   = epsMT
		elseif(holdAmps && holdWidths && holdAmps && !holdModTheta) 
			eps[i + 0]   = 0		
			eps[i + 1]   = 0
			eps[i + 2]   = 0
			eps[i + 10]   = epsMT
		endif
		
		eps[i + 3]   = 0
		eps[i + 4]   = 0
		eps[i + 5]   = 0
		eps[i + 6]   = 0
		eps[i + 7]   = 0
		eps[i + 8]   = 0
		eps[i + 9]   = 0
	endfor
	
	return eps

End

Function/WAVE makeEpsiltonWaveAlphaFit(pw,Emax)
	
	Wave pw
	Variable Emax
	
	Variable i,nPeaks = (numpnts(pw)-pw[0]-4)/11
	//Set up values for epsilon wave
	Make/O/N=(numpnts(pw)) eps_alphaFit
	Wave eps = eps_alphaFit
	
	Variable epsAmp  = 1e-4
	Variable epsAmp2 = 1e-4
	eps[0] = 0
	eps[1] = 0//5e-9
	eps[2] = 1e-2
	eps[3] = 0
	
	for(i=4;i<pw[0] + 4;i+=1)
		eps[i]   = 0	//For Ionization Potentials from DFT. Used to build the Step Edge
	endfor
	
	for(i=4 + pw[0];i<11*nPeaks+4+pw[0];i+=11)
		 
		eps[i + 0]    = 0	
		eps[i + 1]    = 0
		if(pw[i] < Emax) 
			eps[i + 2]    = epsAmp
		else
			eps[i + 2]    = 0
		endif
		eps[i + 3]   = 0
		eps[i + 4]   = 0
		eps[i + 5]   = 0
		eps[i + 6]   = 0
		eps[i + 7]   = 0
		eps[i + 8]   = 0
		eps[i + 9]   = 0
		eps[i + 10]  = 0
	endfor
	
	return eps
End

Function/WAVE makeEpsiltonWaveI0Fit(pw)
	
	Wave pw
	
	Variable i,nPeaks = (numpnts(pw)-pw[0]-4)/11
	//Set up values for epsilon wave
	Make/O/N=(numpnts(pw)) eps_I0Fit
	Wave eps = eps_I0Fit
	
	eps[0] = 0
	eps[1] = 1e-2
	eps[2] = 0
	eps[3] = 0
	
	for(i=4;i<pw[0] + 4;i+=1)
		eps[i]   = 0	//For Ionization Potentials from DFT. Used to build the Step Edge
	endfor
	
	for(i=4 + pw[0];i<11*nPeaks+4+pw[0];i+=11)		 
		eps[i + 0]   = 0	
		eps[i + 1]   = 0
		eps[i + 2]   = 0
		eps[i + 3]   = 0
		eps[i + 4]   = 0
		eps[i + 5]   = 0
		eps[i + 6]   = 0
		eps[i + 7]   = 0
		eps[i + 8]   = 0
		eps[i + 9]   = 0
		eps[i + 10]  = 0
	endfor	
	return eps
End

Function/WAVE makeConstraintWave(holdAmps,holdWidths,holdPos,nPeaks,pWave1D,holdModTheta)
	
	Variable holdAmps,holdWidths,holdPos,holdModTheta
	Variable nPeaks
	Wave pWave1D
	
	Variable i,j=0,pos,pLow,pHigh,wid,n=1
	//Variable openParams = holdAmps + holdWidths + holdPos
	Make/O/T/N=(nPeaks) constraints
	
	String lc0,lc1,lc2,lc3,hc0,lowConstraint2,lowConstraint3
	
		if(!holdPos && !holdWidths && !holdAmps)
			Redimension/N=(4*nPeaks) constraints//(5*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid 
				lc0   = "K" + num2str(i+0)  + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0)  + " <=" + num2str(pHigh)  
				lc1   = "K" + num2str(i+1)  + " >= 0"
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc0
				constraints[j+1]  = hc0
				constraints[j+2]  = lc1
				constraints[j+3]  = lc2
		//		constraints[j+4] = lc3
				j+=4//5
			endfor
		elseif(holdPos && !holdWidths && holdAmps)  
			Redimension/N=(1*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				lc1   = "K" + num2str(i+1) + " >= 0"
				constraints[j+0]   = lc1
				j+=1
			endfor
		elseif(holdPos && !holdWidths && !holdAmps)  
			Redimension/N=(2*nPeaks) constraints//(3*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				lc1   = "K" + num2str(i+1)  + " >= 0"
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc1
				constraints[j+1]  = lc2
		//		constraints[j+2] = lc3
				j+=2//3
			endfor
		elseif(holdPos && holdWidths && !holdAmps)  
			Redimension/N=(1*nPeaks) constraints//(2*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1] = lc3
				j+=1//2
			endfor
		elseif(!holdPos && holdWidths && !holdAmps) 
			Redimension/N=(3*nPeaks) constraints//(4*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0)  + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0)  + " <=" + num2str(pHigh) 
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc0
				constraints[j+1]  = hc0
				constraints[j+2]  = lc2
			//	constraints[j+3] = lc3
				j+=3//4
			endfor
		elseif(!holdPos && !holdWidths && holdAmps)
			Redimension/N=(3*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0) + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0) + " <=" + num2str(pHigh)
				lc1   = "K" + num2str(i+1) + " >= 0"
				constraints[j+0] = lc0
				constraints[j+1] = hc0
				constraints[j+2] = lc1
				j+=3
			endfor
		elseif(!holdPos && holdWidths && holdAmps)
			Redimension/N=(2*nPeaks) constraints
			for(i=4 + pWave1D[0];i<11*nPeaks + 4 + pWave1D[0];i+=11)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0) + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0) + " <=" + num2str(pHigh)
				constraints[j+0] = lc0
				constraints[j+1] = hc0
				j+=2
			endfor			
		elseif(holdAmps && holdWidths && holdAmps) 
			print "Constraint wave not built due to all parameters being held constant"
		endif		
	
	return constraints
End

Function alphaModelPlotting(fit,allExpEnergy,allExpSpec,allStepEnergy,allDFTSteps,results,totalNEXAFS,nSpec,nPeaks,pwCopy,tval,ovpVal,thetaList,d,NEXAFStype,[alpha])

	String fit,thetaList,NEXAFStype
	Wave allExpEnergy,allExpSpec,allStepEnergy,allDFTSteps,results,totalNEXAFS,pwCopy
	Variable nSpec,nPeaks,tval,ovpVal,d,alpha
	
	if(StringMatch(fit,"yes"))
		splitSpec(allExpEnergy ,nSpec,"expEnergy" ,alpha=round(pwCopy[2]))
		splitSpec(allExpSpec   ,nSpec,"expSpec"   ,alpha=round(pwCopy[2]))
		splitSpec(allStepEnergy,nSpec,"stepEnergy",alpha=round(pwCopy[2]))
		splitSpec(allDFTSteps  ,nSpec,"dftStep"   ,alpha=round(pwCopy[2]))
		splitSpec(results  ,nSpec,"fitresults"    ,alpha=round(pwCopy[2]))
	//	if(d)
	//		print "Generating comparison plots from fits." 
	//		plotResultsFit2("expSpec;fitresults;dftStep;",pwCopy[2],"expEnergy",thetaList,tval,ovpVal,nPeaks,NEXAFStype)
	//	endif
	elseif(StringMatch(fit,"no"))
		splitSpec(allExpEnergy ,nSpec,"expEnergy" ,alpha=alpha)
		splitSpec(allExpSpec   ,nSpec,"expSpec"   ,alpha=alpha)
		splitSpec(allStepEnergy,nSpec,"stepEnergy",alpha=alpha)
		splitSpec(allDFTSteps  ,nSpec,"dftStep"   ,alpha=alpha)//,reScale=1)
		splitSpec(totalNEXAFS  ,nSpec,"dftSpec"   ,alpha=alpha)
		if(d)
			print "Generating comparison plots from modeling."
		//	plotPeaks(nSpec,nPeaks,pwCopy[2],tval,ovpVal)
			plotResults("dftSpec;",alpha,"expEnergy",thetaList,tval,ovpVal)
		endif	
	endif
End

Function findMagicAngle(thetaList)

	String thetaList
	
	Variable nThetas = ItemsInList(thetaList),i
	for(i=0;i<nThetas;i+=1)
		Variable thetaVal = str2num(StringFromList(i,thetaList))
		if((54 <= thetaVal) && (thetaVal <=55))
			break 
		endif
	endfor
	
	return i+1
End

Function/WAVE makeMaskWave(E1,E2,xw)
	Variable E1,E2	//Energies defining region to be fit. Energies less/greater than E1/E2 will be masked
	Wave xw	//Energy wave
	//Fit only region between energy 1 and energy 2
	//If maskwave point equals 0 then point won't be fitted.
	Duplicate/O xw,maskWave
	Variable i,x=numpnts(xw)
	for(i=0;i<x;i+=1)
		Variable en = xw[i]
		if(en >= E1 && en <= E2)
			maskWave[i] = 1
		else
			maskWave[i] = 0
		endif
	endfor
	
	return maskWave
End

Function simpleDFTModel(pw,norm3D,thetaValues)

	Wave pw	//1D Parameter wave. Has peak positions,widths, and amplitudes. has initial guess for alpha and IPs required to build step edge
	Wave norm3D //Absorption tensor	
	Wave thetaValues
	Variable i,j=0,k,m
	//Make the absorption tensor for each transition, normalize it, and place each transition into a layer of a 3D wave
	
	Variable nTransitions = DimSize(norm3D,2)
	Variable pos,wid,amp,alpha,i0,phi	
	Variable p0,pf,nSpec=numpnts(thetaValues),cSpec=1,thetaVal	
	Variable nAtoms = pw[0]
	Variable targetDP = 10

	Make/O/D/N=(3,3) currentTensor = 0
	//This is the wave that will contain the DFT tensor model NEXAFS
	Variable specpnts = 300*nSpec
	Variable ini = 280,fin=290,res=2000
	Wave xWave = makeX(ini,fin,specpnts,nSpec)
	Wave thw = makeThetaWave(ini,fin,specpnts,nSpec,thetaValues)
	Variable nPnts = specpnts/nSpec
	Make/O/D/N=(nTransitions) MA = 0 
	
	//Make rotation matrices for rotation of 90,180,and 270 degrees around z-axis
	Make/O/D/N=(3,3) rotMat0   = {{cos(0*(Pi/180))  ,sin(0*(Pi/180)),0}  ,{-sin(0*(Pi/180))  ,cos(0*(Pi/180)),0}  ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat90  = {{cos(90*(Pi/180)) ,sin(90*(Pi/180)),0} ,{-sin(90*(Pi/180)) ,cos(90*(Pi/180)),0} ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat180 = {{cos(180*(Pi/180)),sin(180*(Pi/180)),0},{-sin(180*(Pi/180)),cos(180*(Pi/180)),0},{0,0,1}}
	Make/O/D/N=(3,3) rotMat270 = {{cos(270*(Pi/180)),sin(270*(Pi/180)),0},{-sin(270*(Pi/180)),cos(270*(Pi/180)),0},{0,0,1}}
	
	for(i=3;i<3+3*nTransitions;i+=3)
		i0    = pw[0]
		alpha = pw[1]//90-pw[i+9]
		phi   = pw[2]
		pos   = pw[i + 0]
		wid   = pw[i + 1]
		amp   = pw[i + 2]
		
		currentTensor = i0 * amp * norm3D[p][q][j]
		
		Make/O/D/N=(3,3) rotMatAlignX = {{1,0,0},{0,cos(alpha*(Pi/180)),sin(alpha*(Pi/180))},{0,-sin(alpha*(Pi/180)),cos(alpha*(Pi/180))}}
		
		MatrixOP/O tiltedTensor = rotMatAlignX x currentTensor x rotMatAlignX^t
		
		MatrixOP/O  RotTensor0   =  rotMat0   x tiltedTensor x rotMat0^t
		MatrixOP/O  RotTensor90  =  rotMat90  x tiltedTensor x rotMat90^t
		MatrixOP/O  RotTensor180 =  rotMat180 x tiltedTensor x rotMat180^t
		MatrixOP/O  RotTensor270 =  rotMat270 x tiltedTensor x rotMat270^t
		
		MatrixOP/O filmTensor = RotTensor0 + RotTensor90 + RotTensor180 + RotTensor270
		
		for(k=0;k<=2;k+=1)
			for(m=0;m<=2;m+=1)
				Variable currComp = filmTensor[k][m]
				targetDP = round(targetDP)
				currComp = round(currComp * (10^targetDP)) / (10^targetDP)
				if(currComp == -0)
					currComp = 0
				endif
				filmTensor[k][m] = currComp
			endfor
		endfor
		
		if(j<nTransitions)
			if(cSpec == 1)
				p0 = 0
				pf = nPnts-1
				thetaVal = thw[p0]
			elseif(cSpec == nSpec)
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0]
			else
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0]
			endif
		endif
		
		String eName = "eField_" + num2str(thetaVal)
		Make/O/D/N=(3) $eName = {sin(thetaVal*(Pi/180))*cos(phi*(Pi/180)),sin(thetaVal*(Pi/180))*sin(phi*(Pi/180)),cos(thetaVal*(Pi/180))}	
		Wave eField = $eName
		MatrixOP/O tempMA = amp * (eField^t x filmTensor x eField)
		
		if(numtype(tempMA[0]) !=0)
			MA[j] = 0
		else
			MA[j] = tempMA[0]
		endif
		//Make totalNEXAFS for current Theta
		String tNXFSname = "totalNEXAFS_" + num2str(thetaVal) 
		Make/O/N=(specpnts) $tNXFSname
		Wave totalNEXAFS = $tNXFSname 
		//Make peaks associated with current theta
		String pkName = "pk" + num2str(j) + "_spec" + num2str(cSpec) + "_Alpha" + num2str(alpha) + "_Theta" + num2str(thetaVal) 
		Make/O/N=(res) $pkName
		Wave z = $pkName
		SetScale/i x,ini,fin,z
		z =  MA[j] * gauss(x,pos,wid)		
		totalNEXAFS += MA[j] * gauss(xWave,pos,wid)
		j+=1
		
		if(j<=nTransitions)
			if(cSpec == 1)
				if(j==(nTransitions) && cSpec < nSpec)
					i=0
					cSpec+=1
					j=0
				endif 			
			else
				if(j==(nTransitions) && cSpec < nSpec)
					i=0
					cSpec+=1
					j=0
				endif 
			endif
		endif		
	endfor
	
	KillWaves rotMat0,rotMat90,rotMat180,rotMat270,rotTensor0,rotTensor90,rotTensor180,rotTensor270,rotMatAlignX
	organizeSimpleModel(alpha)
	plotSimpleTensor(pw)
	Wave enWave = makeX2(ini,fin,res)
	plotAlphaIntensities(pw,thetaValues,enWave)
	concatenateEField()
End


Function organizeSimpleModel(alpha)

	Variable alpha
	
	String iniFolder = GetDataFolder(1)
	String fName = "Alpha_" + num2str(alpha)
	NewDataFolder/O $fName
	
	String pkWaves = WaveList("pk*",";","")
	Variable npks = ItemsInList(pkWaves),i
	for(i=0;i<npks;i+=1)
		String name = StringFromList(i,pkWaves)
		Wave w = $name
		String nName = iniFolder + fName + ":" + name
		Duplicate/O w,$nName
		KillWaves w
	endfor
End

Function plotSimpleTensor(pWave)
	
	Wave pWave
	
	Variable alpha = pWave[1]
	Variable energy = pWave[3],width = pWave[4],i
	String iniFolder = GetDataFolder(1)
	String fName = "Alpha_" + num2str(alpha)
	SetDataFolder $fName
	String wList = WaveList("pk*",";","")
	Variable n = ItemsInList(wList)
	String plotName = "Theta_Peaks"//"Alpha_" + num2str(alpha) + "_PLOT"
	DoWindow $plotName
	if(!V_Flag)
		Display/N=$plotName/K=1 
		for(i=0;i<n;i+=1)
			String wName = StringFromList(i,wList)
			Wave w = $wName
			AppendToGraph/W=$plotName w 
		endfor
		SetAxis bottom energy-3*width,energy+3*width
		ModifyGraph grid=2,mirror=1,nticks(left)=10,minor=1,fStyle=1,fSize=14,lsize=1.5
		Label left "Transition Intensity [a.u.]"
		Label bottom "Transition Energy [eV]"
		ApplyColorTableToTopGraph("ColdWarm")
		TextBox/C/N=text0 "Alpha " + num2str(alpha)
	endif 
	SetDataFolder $iniFolder
End

Function plotAlphaIntensities(pWave,tw,ew)
	
	Wave pWave,tw,ew
	
	Variable pt = BinarySearch(ew,pWave[3])
	Variable alpha = pWave[1]
	String iniFolder = GetDataFolder(1)
	String fName = "Alpha_" + num2str(alpha)
	SetDataFolder $fName
	String wList = WaveList("pk*",";","")
	Variable n = ItemsInList(wList),i
	String intName = "intensity_Alpha" + num2str(alpha)
	Make/O/N=(n) $intName
	Wave intensity = $intName
	for(i=0;i<n;i+=1)
		String name = StringFromList(i,wList)
		Wave w = $name
		intensity[i] = w[pt]
	endfor
	
	String plotName = "Theta_Intensities"//"alphaIntensities_" + num2str(alpha)
	DoWindow $plotName
	if(!V_Flag)
		Display/N=$plotName/K=1 intensity vs tw 
		ModifyGraph grid=2,mirror=1,nticks(left)=10,minor=1,fStyle=1,fSize=14,lsize=1.5,mode=4,marker=19,msize=3
		Label left "Transition Intensity [a.u.]"
		Label bottom "Sample θ [°]"
	endif
	
	SetDataFolder $iniFolder
End
 
Function plotAlphaIntensities2(en,ew,nxname,nSpec,thetaList)
	
	Wave ew
	Variable en,nSpec
	String nxName,thetaList
	
	Variable pt = BinarySearch(ew,en),i
	Make/O/N=(nSpec) intensity,tw
	String nxList = WaveList(nxname+"*",";","")
	for(i=0;i<nSpec;i+=1)
		String name = StringFromList(i,nxList)
		Wave w = $name
		intensity[i] = w[pt]
		tw[i] = str2num(StringFromList(i,thetaList))
	endfor
	
	String plotName = "Theta_Intensities"//"alphaIntensities_" + num2str(alpha)
	DoWindow $plotName
	if(!V_Flag)
		Display/N=$plotName/K=1 intensity vs tw 
		ModifyGraph grid=2,mirror=1,nticks(left)=10,minor=1,fStyle=1,fSize=14,lsize=1.5,mode=4,marker=19,msize=3
		Label left "Transition Intensity [a.u.]"
		Label bottom "Sample θ [°]"
	endif
	
End

Function/WAVE makeX(ini,fin,res,nSpec)

	Variable ini,fin,res,nSpec
	
	Make/O/N=(res) xWave = 0
	Variable i,j=1,c=1
	for(i=0;i<res;i+=1)
		xWave[i] = ini + j*(fin-ini)/(res/nSpec)
		j+=1
		if(mod(i+1,(res/nSpec))==0 && i!=0)
			j=0
			c+=1
			if(c<=nSpec)
				continue
			else
				break
			endif
		endif
	endfor
	
	return xWave
End

Function/WAVE makeThetaWave(ini,fin,res,nSpec,thetaValues)

	Variable ini,fin,res,nSpec
	Wave thetaValues
	
	Make/O/N=(res) thetaWave = 0
	Variable i,c=0
	for(i=0;i<res;i+=1)
		Variable thetaVal = thetaValues[c]
		thetaWave[i] = thetaVal
		if(mod(i+1,(res/nSpec))==0 && i!=0)
			c+=1
			if(c<=nSpec)
				continue
			else
				break
			endif
		endif
	endfor
	
	return thetaWave
End

Function/WAVE makeX2(ini,fin,res)

	Variable ini,fin,res
	
	Make/O/N=(res) enWave = 0
	Variable i
	for(i=0;i<res;i+=1)
		enWave[i] = ini + i*(fin-ini)/res
	endfor
	
	return enWave
End

Function concatenateEField()

	String eFieldList = WaveList("eField*",";","")
	Concatenate eFieldList,AllEFIelds
	Variable i,n=ItemsInList(eFieldList)
	for(i=0;i<n;i+=1)
		String name = StringFromList(i,eFieldList)
		Wave w = $name
		KillWaves w
	endfor

End

Function fitSinglePeak(pw,pkID,holdPos,holdWidths,holdAmps,holdSecAmp,tval,ovpMax,alpha,mol)
	
	Wave pw
	Variable pkId
	Variable holdPos,holdWidths,holdAmps,holdSecAmp
	Variable tval,ovpMax
	Variable alpha 
	String mol
	
	String baseFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":TransitionFiltering_"+ replaceString(".",num2str(tval),"p") + "OS_" + replaceString(".",num2str(ovpMax),"p") +"OVP:AmplitudeFitting1:" 
	String iniFolder = baseFolder + "Alpha_" + replaceString(".",num2str(alpha),"p")
	String tensorFolder = iniFolder + ":TensorMisc:"
	SetDataFolder $tensorFolder
	
	//TO DO:
	//1. Find the parameter wave for the fits
	String newpwName = "refine_pwCopy_pk" + num2str(pkID)
	Duplicate/O pw,$newpwName
	Wave newpw = $newpwName
	
	//2. Find the experimental NEXAFS to be fit
	String specName = tensorFolder + "allExpSpec"
	Wave spec = $specName
	String resultsName = "refine_results_pk" + num2str(pkID) 
	String residualName = "refine_res_pk" + num2str(pkID)
	Duplicate/O spec,$resultsName,$residualName
	Wave results = $resultsName
	Wave res = $residualName
	
	//3. Find the experimental energies
	String enName = tensorFolder + "allExpEnergy"
	Wave en = $enName
	
	//4.Get wave with sample theta values
	Wave valuesPolar
	
	//5.Get wave with DFT step edge
	Wave allDFTSteps
	
	//6. Modify hold string ,constraint wave and epsilon wave so that only the parameters of the current peak are modified
	Wave constraints = makeSingleConstraintWave(holdAmps,holdWidths,holdPos,holdSecAmp,pw,pkID)
	Wave epsilon = makeSingleEpsiltonWave(pw,pkID,holdAmps,holdWidths,holdPos,holdSecAmp)
	String H = makeSingleHoldString(holdAmps,holdWidths,holdPos,holdSecAmp,pw,pkID)
	
	//Variable V_FitError = 0
	Variable V_FitTol = 0.00001
	if(holdAmps && holdWidths && holdPos)
		FuncFit/H=H/M=2 simDFTfit2,newpw, spec /X={en,valuesPolar,allDFTSteps} /R=res /E=epsilon /D=results
	else
		FuncFit/H=H/M=2 simDFTfit2,newpw, spec /X={en,valuesPolar,allDFTSteps} /R=res /E=epsilon /C=constraints /D=results
	endif
	
	splitSpec(results,4,"refine_results_pk" + num2str(pkID) +"_",alpha=alpha)
End

//This function makes the constraint wave for the fit refinement.
//The amplitudes of all peaks remain open but the position,width, amplitude, local alpha, and tensor element ratio can be opened for a desired peak.
Function/WAVE makeSingleConstraintWave(holdAmps,holdWidths,holdPos,holdSecAmp,pWave1D,pkID)
	
	Variable holdAmps,holdWidths,holdPos,holdSecAmp
	Variable pkID
	Wave pWave1D
	
	Variable i=4 + pWave1D[0] + pkID*11,j=2,pos,pLow,pHigh,wid,n=2
	Variable nPeaks = ((numpnts(pWave1D) - 4 -  pWave1D[0])/11)-1,cpk=0
	
	Variable totalConstraints = 0
	if(!holdAmps)
		totalConstraints +=1	
	endif
	
	if(!holdWidths)
		totalConstraints +=1	
	endif
	
	if(!holdPos)
		totalConstraints +=2	
	endif
	
	if(!holdSecAmp)
		totalConstraints +=1	
	endif
	
	String constraintWName = "constraints_Pk" + num2str(pkID)
	Make/O/T/N=(nPeaks + totalConstraints) $constraintWName
	Wave/T constraints = $constraintWName
	//Change amplitude such that it's the absolute value.
	String lc0,lc1,lc2,lc3,hc0
	//lc0 --> Controls the lower bound for the peak position
	//hc0 --> Controls the upper bound for the peak position
	//lc1 --> Controls the lower bound for the peak width (always be greater than 0)
	//lc2 --> Controls the lower bound on the peak amplitude (Always greater than 0)
	//lc3 --> Controls the lower bound for the ratio between the tensor elements (always greater than zero)
	//lc4 --> Controls the lower bound for the tilt angle (always greater than or equal to 0 degrees)
	//lc5 --> Controls the upper bound on the tilt angle (always less than or equal to 90 degrees)
	if(!holdPos && !holdWidths && !holdAmps && !holdSecAmp)
		Redimension/N=(nPeaks + totalConstraints) constraints//(5*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid 
				lc0   = "K" + num2str(i+0)  + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0)  + " <=" + num2str(pHigh)  
				lc1   = "K" + num2str(i+1)  + " >= 0"
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc0
				constraints[j+1]  = hc0
				constraints[j+2]  = lc1
				constraints[j+3]  = lc2
				constraints[j+4]  = lc3
				j+=5
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
				constraints[j+0]  = lc2
				j+=1
				cpk+=1
			endif
		endfor
	elseif(holdPos && !holdWidths && holdAmps && !holdSecAmp)  
		Redimension/N=(nPeaks + totalConstraints) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				lc1   = "K" + num2str(i+1) + " >= 0"
				constraints[j+0]   = lc1
				j+=1
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor	
	elseif(holdPos && !holdWidths && !holdAmps && !holdSecAmp)  
		Redimension/N=(nPeaks + totalConstraints) constraints//(3*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				lc1   = "K" + num2str(i+1)  + " >= 0"
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc1
				constraints[j+1]  = lc2
			//	constraints[j+2] = lc3
				j+=2
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor	
	elseif(holdPos && holdWidths && !holdAmps && !holdSecAmp)  
		Redimension/N=(nPeaks + totalConstraints) constraints//(2*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1] = lc3
				j+=1
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && holdWidths && !holdAmps && !holdSecAmp) 
		Redimension/N=(nPeaks + totalConstraints) constraints//(4*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0)  + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0)  + " <=" + num2str(pHigh) 
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc0
				constraints[j+1]  = hc0
				constraints[j+2]  = lc2
				constraints[j+3] = lc3
				j+=4
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && !holdWidths && holdAmps && !holdSecAmp)
		Redimension/N=(nPeaks + totalConstraints) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0) + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0) + " <=" + num2str(pHigh)
				lc1   = "K" + num2str(i+1) + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0] = lc0
				constraints[j+1] = hc0
				constraints[j+2] = lc1
				constraints[j+3] = lc3
				j+=4
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && holdWidths && holdAmps && !holdSecAmp)
		Redimension/N=(nPeaks + totalConstraints) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0) + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0) + " <=" + num2str(pHigh)
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0] = lc0
				constraints[j+1] = hc0
				constraints[j+2] = lc3
				j+=3
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && !holdWidths && !holdAmps && holdSecAmp)
		Redimension/N=(nPeaks + totalConstraints) constraints//(5*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid 
				lc0   = "K" + num2str(i+0)  + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0)  + " <=" + num2str(pHigh)  
				lc1   = "K" + num2str(i+1)  + " >= 0"
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc0
				constraints[j+1]  = hc0
				constraints[j+2]  = lc1
				constraints[j+3]  = lc2
			//	constraints[j+4] = lc3
				j+=4
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(holdPos && !holdWidths && holdAmps && holdSecAmp)  
		Redimension/N=(nPeaks + totalConstraints) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				lc1   = "K" + num2str(i+1) + " >= 0"
				constraints[j+0]   = lc1
				j+=1
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor	
	elseif(holdPos && !holdWidths && !holdAmps && holdSecAmp)  
		Redimension/N=(nPeaks + totalConstraints) constraints//(3*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				lc1   = "K" + num2str(i+1)  + " >= 0"
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc1
				constraints[j+1]  = lc2
			//	constraints[j+2] = lc3
				j+=2
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor	
	elseif(holdPos && holdWidths && !holdAmps && holdSecAmp)  
		Redimension/N=(nPeaks + totalConstraints) constraints//(2*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1] = lc3
				j+=1
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && holdWidths && !holdAmps && holdSecAmp) 
		Redimension/N=(nPeaks + totalConstraints) constraints//(4*nPeaks) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0)  + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0)  + " <=" + num2str(pHigh) 
				lc2   = "K" + num2str(i+2)  + " >= 0"
				lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc0
				constraints[j+1]  = hc0
				constraints[j+2]  = lc2
			//	constraints[j+3] = lc3
				j+=3
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && !holdWidths && holdAmps && holdSecAmp)
		Redimension/N=(nPeaks + totalConstraints) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0) + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0) + " <=" + num2str(pHigh)
				lc1   = "K" + num2str(i+1) + " >= 0"
				constraints[j+0] = lc0
				constraints[j+1] = hc0
				constraints[j+2] = lc1
				j+=3
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(!holdPos && holdWidths && holdAmps && holdSecAmp)
		Redimension/N=(nPeaks + totalConstraints) constraints
		j=0
		for(i=4 + pWave1D[0];i < 4 + pWave1D[0] + 11*(nPeaks+1);i+=11)
			if(cpk == pkID)
				pos   = pWave1D[i+0]
				wid   = pWave1D[i+1]
				pLow  = pos - n*wid
				pHigh = pos + n*wid						
				lc0   = "K" + num2str(i+0) + " >=" + num2str(pLow)
				hc0   = "K" + num2str(i+0) + " <=" + num2str(pHigh)
				constraints[j+0] = lc0
				constraints[j+1] = hc0
				j+=2
				cpk+=1
			else
				lc2   = "K" + num2str(i+2)  + " >= 0"
			//	lc3   = "K" + num2str(i+10) + " >= 0"
				constraints[j+0]  = lc2
			//	constraints[j+1]  = lc3
				j+=1
				cpk+=1
			endif
		endfor
	elseif(holdAmps && holdWidths && holdAmps && holdSecAmp) 
		print "Constraint wave not built due to all parameters being held constant"
	endif		
	
	return constraints		
End

Function/WAVE makeSingleEpsiltonWave(pWave1D,pkID,holdAmps,holdWidths,holdPos,holdSecAmp)
	
	Wave pWave1D
	Variable pkID
	Variable holdAmps,holdWidths,holdPos,holdSecAmp
	
	Variable nPeaks = (numpnts(pWave1D) - 4 - pWave1D[0])/11,cPk=0
	Variable i
	//Set up values for epsilon wave
	String epsWaveName = "eps_Pk" + num2str(pkID)
	Make/O/N=(numpnts(pWave1D)) $epsWaveName
	Wave eps = $epsWaveName
	
	Variable epsPos  = 1e-7
	Variable epsWid  = 1e-7
	Variable epsAmp  = 1e-4
	Variable epsProp = 1e-4
	eps[0] = 0
	eps[1] = 0//5e-9
	
	eps[3] = 0
	
	for(i=4;i<pWave1D[0] + 4;i+=1)
		eps[i]   = 0	//For Ionization Potentials from DFT. Used to build the Step Edge
	endfor
	
	for(i=4 + pWave1D[0];i<11*nPeaks+4+pWave1D[0];i+=11)
		
		if(cPk == pkID)//i==(11*pkID + 4 + pWave1D[0]))
			if(!holdPos && !holdWidths && !holdAmps && holdSecAmp)   
				eps[i + 0]   = epsPos		
				eps[i + 1]   = epsWid
				eps[i + 2]   = epsAmp
				eps[i + 10]   = 0
			elseif(holdPos && !holdWidths && holdAmps && holdSecAmp)  
				eps[i + 0]   = 0	
				eps[i + 1]   = epsWid
				eps[i + 2]   = 0
				eps[i + 10]   = 0
			elseif(holdPos && !holdWidths && !holdAmps && holdSecAmp)  
				eps[i + 0]   = 0		
				eps[i + 1]   = epsWid
				eps[i + 2]   = epsAmp
				eps[i + 10]   = 0
			elseif(holdPos && holdWidths && !holdAmps && holdSecAmp)  
				eps[i + 0]   = 0	
				eps[i + 1]   = 0
				eps[i + 2]   = epsAmp
				eps[i + 10]   = 0
			elseif(!holdPos && holdWidths && !holdAmps && holdSecAmp) 
				eps[i + 0]   = epsPos		
				eps[i + 1]   = 0
				eps[i + 2]   = epsAmp
				eps[i + 10]   = 0
			elseif(!holdPos && !holdWidths && holdAmps && holdSecAmp)
				eps[i + 0]   = epsPos		
				eps[i + 1]   = epsWid
				eps[i + 2]   = 0
				eps[i + 10]   = 0
			elseif(!holdPos && holdWidths && holdAmps && holdSecAmp)
				eps[i + 0]   = epsPos		
				eps[i + 1]   = 0
				eps[i + 2]   = 0
				eps[i + 10]   = 0
			elseif(holdPos && holdWidths && holdAmps && holdSecAmp) 
				eps[i + 0]   = 0		
				eps[i + 1]   = 0
				eps[i + 2]   = 0
				eps[i + 10]   = 0
			elseif(!holdPos && !holdWidths && !holdAmps && !holdSecAmp)   
				eps[i + 0]   = epsPos		
				eps[i + 1]   = epsWid
				eps[i + 2]   = epsAmp
				eps[i + 10]   = epsProp
			elseif(holdPos && !holdWidths && holdAmps && !holdSecAmp)  
				eps[i + 0]   = 0	
				eps[i + 1]   = epsWid
				eps[i + 2]   = 0
				eps[i + 10]   = epsProp
			elseif(holdPos && !holdWidths && !holdAmps && !holdSecAmp)  
				eps[i + 0]   = 0		
				eps[i + 1]   = epsWid
				eps[i + 2]   = epsAmp
				eps[i + 10]   = epsProp
			elseif(holdPos && holdWidths && !holdAmps && !holdSecAmp)  
				eps[i + 0]   = 0	
				eps[i + 1]   = 0
				eps[i + 2]   = epsAmp
				eps[i + 10]   = epsProp
			elseif(!holdPos && holdWidths && !holdAmps && !holdSecAmp) 
				eps[i + 0]   = epsPos		
				eps[i + 1]   = 0
				eps[i + 2]   = epsAmp
				eps[i + 10]   = epsProp
			elseif(!holdPos && !holdWidths && holdAmps && !holdSecAmp)
				eps[i + 0]   = epsPos		
				eps[i + 1]   = epsWid
				eps[i + 2]   = 0
				eps[i + 10]   = epsProp
			elseif(!holdPos && holdWidths && holdAmps && !holdSecAmp)
				eps[i + 0]   = epsPos		
				eps[i + 1]   = 0
				eps[i + 2]   = 0
				eps[i + 10]   = epsProp
			elseif(holdPos && holdWidths && holdAmps && !holdSecAmp) 
				eps[i + 0]    = 0		
				eps[i + 1]    = 0
				eps[i + 2]    = 0
				eps[i + 10]   = epsProp
			endif
		
			eps[i + 3]    = 0
			eps[i + 4]    = 0
			eps[i + 5]    = 0
			eps[i + 6]    = 0
			eps[i + 7]    = 0
			eps[i + 8]    = 0
			eps[i + 9]    = 0
		else
			eps[i + 0]    = 0		
			eps[i + 1]    = 0
			eps[i + 2]    = 0//epsAmp
			eps[i + 3]    = 0
			eps[i + 4]    = 0
			eps[i + 5]    = 0
			eps[i + 6]    = 0
			eps[i + 7]    = 0
			eps[i + 8]    = 0
			eps[i + 9]    = 0
			eps[i + 10]   = 0
		endif
		cPk+=1
	endfor
	
	return eps

End

Function/S makeSingleHoldString(holdAmps,holdWidths,holdPos,holdSecAmp,pWave1D,pkID)
	Variable holdAmps,holdWidths,holdPos,holdSecAmp
	Wave pWave1D
	Variable pkID
	
	Variable i
	Variable nPeaks = (numpnts(pWave1D) - 4 - pWave1D[0])/11
	//Make Hold String for fit. Open amplitudes, hold position and width constant
	String H = "1111"//H1 = Number of atoms, i0, alpha, phi	
	
	for(i=1;i<=pWave1D[0];i+=1)
		H +="1" //For Ionization Potentials from DFT. Used to build the Step Edge. Hold them constant.
	endfor
	
	for(i=0;i<=nPeaks-1;i+=1)
	//Hold string for peak parameters. Position, Width, Max Amplitude, xxNoRM, xyNorm,xzNorm,yyNorm,yzNorm,zzNorm
	//Hold cases:
		if(i==pkID)
			if(!holdPos && !holdWidths && !holdAmps && holdSecAmp)   
				H +="00011111111"
			elseif(holdPos && !holdWidths && holdAmps && holdSecAmp)  
				H +="10111111111"
			elseif(holdPos && !holdWidths && !holdAmps && holdSecAmp)  
				H +="10011111111"
			elseif(holdPos && holdWidths && !holdAmps && holdSecAmp)  
				H +="11011111111"
			elseif(!holdPos && holdWidths && !holdAmps && holdSecAmp) 
				H +="01011111111"
			elseif(!holdPos && !holdWidths && holdAmps && holdSecAmp)
				H +="00111111111"
			elseif(!holdPos && holdWidths && holdAmps && holdSecAmp)
				H +="01111111111"
			elseif(holdAmps && holdWidths && holdAmps && holdSecAmp) 
				H +="11111111111"
			//
			elseif(!holdPos && !holdWidths && !holdAmps && !holdSecAmp)   
				H +="00011111110"
			elseif(holdPos && !holdWidths && holdAmps && !holdSecAmp)  
				H +="10111111110"
			elseif(holdPos && !holdWidths && !holdAmps && !holdSecAmp)  
				H +="10011111110"
			elseif(holdPos && holdWidths && !holdAmps && !holdSecAmp)  
				H +="11011111110"
			elseif(!holdPos && holdWidths && !holdAmps && !holdSecAmp) 
				H +="01011111110"
			elseif(!holdPos && !holdWidths && holdAmps && !holdSecAmp)
				H +="00111111110"
			elseif(!holdPos && holdWidths && holdAmps && !holdSecAmp)
				H +="01111111110"
			elseif(holdAmps && holdWidths && holdAmps && !holdSecAmp) 
				H +="11111111110"	
			//	
			endif
		else
			H +="11111111111"//"11011111111"
		endif
	endfor
	
	String holdWaveName = "holdWave_pk" + num2str(pkID)
	Make/O/N=(numpnts(pWave1D)) $holdWaveName
	Wave holdWave = $holdWaveName
	for(i=0;i<numpnts(pWave1D);i+=1)
		holdWave[i] = str2num(H[i,i])
	endfor
	
	print H
	return H
End

Function compareTransitionPositions(pw1,pw2)

	Wave pw1,pw2
	
	Variable i,nPeaks = (numpnts(pw1)-4-pw1[0])/11,en1,en2,amp1,amp2,j=0
	Make/O/N=(nPeaks) iniEnergy,finEnergy,diffEnergy,iniAmp,finAmp,diffAmp
	print nPeaks
	for(i=4+pw1[0];i<11*nPeaks+4+pw1[0];i+=11)
		en1 = pw1[i] 
		en2 = pw2[i]
		amp1 = pw1[i+2]
		amp2 = pw2[i+2]
		iniEnergy[j] = en1
		finEnergy[j] = en2
		iniAmp[j] = amp1
		finAmp[j] = amp2
		diffEnergy[j] = en2-en1
		diffAmp[j] = amp2 - amp1
		j+=1
	endfor
End

Function/WAVE modelVectorOrbital(tw,a)

	Wave tw
	Variable a
	
	Variable i,x=numpnts(tw),alphaR = a * (Pi/180)
	Make/O/N=(x) vecOrb,vecParaInt,vecPerpInt	
	
	Variable para,perp
	for(i=0;i<x;i+=1)
		Variable t_angle = tw[i] * (Pi/180)
		para = (1/3) * (1 + (1/2)*((3 * cos( t_angle)^2 - 1) * (3 * cos(alphaR)^2 - 1)))
		perp = (1/2) * sin(alphaR)^2
		vecParaInt[i] = para
		vecPerpInt[i] = perp	
		vecOrb[i] = para + perp
	endfor

	return vecOrb
End

Function/WAVE modelPlanarOrbital(tw,a)

	Wave tw
	Variable a
	
	Variable i,x=numpnts(tw),alphaR = a * (Pi/180)
	Make/O/N=(x) plOrb,plParaInt,plPerpInt	
	
	Variable para,perp
	for(i=0;i<x;i+=1)
		Variable t_angle = tw[i] * (Pi/180)
		para = (2/3) * (1 - (1/4)*((3 * cos( t_angle)^2 - 1) * (3 * cos(alphaR)^2 - 1)))
		perp = (1/2) * ( 1 + cos(alphaR)^2)
		plParaInt[i] = para
		plPerpInt[i] = perp	
		plOrb[i] = para + perp
	endfor
	
	return plOrb
End

Function orbitalWrapper(tw,a)

	Wave tw
	Variable a
	
	Wave vecOrb = modelVectorOrbital(tw,a)
	Wave plOrb = modelPlanarOrbital(tw,a)
	
	Wave plParaInt,vecParaInt,plPerpInt,vecPerpInt
	WaveStats/Q plParaInt
	Variable plPaMin = V_min,plPaMax = V_max
	WaveStats/Q vecParaInt
	Variable vPaMin = V_min,vPaMax = V_max
	WaveStats/Q plPerpInt
	Variable plPeMin = V_min,plPeMax = V_max
	WaveStats/Q vecPerpInt
	Variable vPeMin = V_min,vPeMax = V_max
	
	Variable perpMin,perpMax,paraMin,paraMax
	if(plPaMin < vPaMin)
		paraMin = plPaMin
	else
		paraMin = vPaMin
	endif
	
	if(plPaMax > vPaMax)
		paraMax = plPaMax
	else
		paraMax = vPaMax
	endif
	
	if(plPeMin < vPeMin)
		perpMin = plPeMin
	else
		perpMin = vPeMin
	endif
	
	if(plPeMax > vPeMax)
		perpMax = plPeMax
	else
		perpMax = vPeMax
	endif
	//Compare Total intensities of vector and planar orbitals
	DoWindow plotVecvsPlanarOrbs
	if(!V_Flag)
		Display/N=plotVecvsPlanarOrbs/K=1 vecOrb,plOrb vs tw
		NewFreeAxis/L para
		NewFreeAxis/L perp
		AppendToGraph/W=plotVecvsPlanarOrbs/L=para plParaInt,vecParaInt vs tw
		AppendToGraph/W=plotVecvsPlanarOrbs/L=perp plPerpInt,vecPerpInt vs tw
		ModifyGraph mode=4,marker=19,lsize=1.5, freePos(para)=0,freePos(perp)=0
		Label left "Total Intensity[a.u.]"
		Label bottom "θ[°]"
		Label para "Para Intensity[a.u.]"
		Label perp "Perp Intensity[a.u.]"
		SetAxis para paraMin,paraMax
		SetAxis perp perpMin,perpMax
		ModifyGraph grid=2,mirror=1,minor=1
		ModifyGraph rgb(plOrb)=(1,4,52428),rgb(plParaInt)=(1,4,52428),rgb(plPerpInt)=(1,4,52428)
		ModifyGraph lblPosMode(left)=1,axisEnab(left)={0.69,1}
		ModifyGraph lblPosMode(para)=1,axisEnab(para)={0.34,0.66}
		ModifyGraph lblPosMode(perp)=1,axisEnab(perp)={0,0.32}
		Legend/C/N=text0/J/A=RB "\\JCα="+num2str(a)+"°\r\\s(plOrb) Planar\r\\s(vecOrb) Vector"
	else
		Legend/C/N=text0/J/W=plotVecvsPlanarOrbs "\\JCα="+num2str(a)+"°\r\\s(plOrb) Planar\r\\s(vecOrb) Vector"
	endif
End

Function/S detOrbType(tensor)
	Wave tensor
	
	String type
	if((tensor[0][0]==tensor[1][1]==tensor[2][2]))
		type = "I"
	elseif((tensor[0][0]==0) && (tensor[1][1]==0) && (tensor[2][2]!=0))
		type = "V"
	elseif((tensor[0][0]!=0) && (tensor[1][1]!=0) && (tensor[2][2]==0))
		type = "P"
	elseif((tensor[0][0]!=0) && (tensor[1][1]!=0) && (tensor[2][2]!=0))
		type = "M"
	endif
	
	print type
	return type
End

Function modelTensorInt(pw,thw,tensor) 

	Wave pw	//1D Parameter wave. Has peak positions,widths, and amplitudes. has initial guess for alpha and IPs required to build step edge
	Wave thw	//Wave containing the various values of theta
	Wave tensor
		
	Variable i,j=0,k,m,targetDP = 10,n=numpnts(thw)
	Variable pos,wid,amp,alpha,i0,phi,localalpha,amp2	
	Duplicate/O tensor,currentTensor
	Duplicate/O pw,npw
	Make/O/N=(n) tensorInt,tempMA
	//Make rotation matrices for rotation of 90,180,and 270 degrees around z-axis
	Make/O/D/N=(3,3) rotMat0   = {{cos(0*(Pi/180))  ,sin(0*(Pi/180)),0}  ,{-sin(0*(Pi/180))  ,cos(0*(Pi/180)),0}  ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat90  = {{cos(90*(Pi/180)) ,sin(90*(Pi/180)),0} ,{-sin(90*(Pi/180)) ,cos(90*(Pi/180)),0} ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat180 = {{cos(180*(Pi/180)),sin(180*(Pi/180)),0},{-sin(180*(Pi/180)),cos(180*(Pi/180)),0},{0,0,1}}
	Make/O/D/N=(3,3) rotMat270 = {{cos(270*(Pi/180)),sin(270*(Pi/180)),0},{-sin(270*(Pi/180)),cos(270*(Pi/180)),0},{0,0,1}}
	
	String vecType = detOrbType(tensor)
	//*****POTENTIALLY CHANGE CLUSTERING ALGORITHM****
	//Generate clusters based on overall/total OS instead of component OS
	//How to open up individual tensor elements for fitting:
	//1. Ratio of in-plane vs out-of-plane components?
	//2. Need a second parameter to describe amplitude? Possibly related to a local alpha?
		//*********Will need to change pWave from 10 parameters per peak to 11!*******
	//3. Visualize tensor elements for each cluster (Tensor Element vs Cluster ID) -->Log plot	
	//4. Additive oscillator strength... How to accomplish this programatically???
	alpha      = npw[0]*(Pi/180)
	phi        = npw[1]*(Pi/180)
	pos        = npw[2]
	wid        = ABS(npw[3])
	amp        = ABS(npw[4])
	amp2       = npw[5]
	
	currentTensor = amp * tensor[p][q]

	Make/O/N=2 vector
	vector[0] = sqrt(tensor[0][0])
	vector[1] = sqrt(0.5*tensor[2][2])
	//vector[2] = tensor[2][2]
	Wave rv = vecOpsWrap(vector,amp2,"Z")
	currentTensor[0][0] = abs(rv[0])^2
	currentTensor[1][1] = abs(rv[0])^2
	currentTensor[2][2] = 2*abs(rv[1])^2
	Make/O/N=(3,3) modMolTensor = 0
	modMolTensor[0][0] = abs(rv[0])^2
	modMolTensor[1][1] = abs(rv[0])^2
	modMolTensor[2][2] = 2*abs(rv[1])^2
	//print sqrt(currentTensor[0][0]^2+currentTensor[1][1]^2+currentTensor[2][2]^2)
	
	Make/O/D/N=(3,3) rotMatAlignX = {{1,0,0},{0,cos(alpha),sin(alpha)},{0,-sin(alpha),cos(alpha)}}
		
	//Tilt the tensor by angle alpha
	MatrixOP/O  tiltedTensorX = rotMatAlignX x currentTensor x rotMatAlignX^t			 
	//Construct the film tensor by adding the azimuthal rotations 
	MatrixOP/O  RotTensor0X   =  rotMat0   x tiltedTensorX x rotMat0^t
	MatrixOP/O  RotTensor90X  =  rotMat90  x tiltedTensorX x rotMat90^t
	MatrixOP/O  RotTensor180X =  rotMat180 x tiltedTensorX x rotMat180^t
	MatrixOP/O  RotTensor270X =  rotMat270 x tiltedTensorX x rotMat270^t
	MatrixOP/O  filmTensor = RotTensor0X + RotTensor90X + RotTensor180X + RotTensor270X
		
	Wave filmTensorTrunc = truncateMatrix(targetDP,filmTensor)
	//print sqrt(filmTensorTrunc[0][0]^2+filmTensorTrunc[1][1]^2+filmTensorTrunc[2][2]^2)	
	//Apply electric field to current tensor to obtain the mass absorbance
	for(i=0;i<n;i+=1)
		Variable thetaVal = thw[i]*(Pi/180)
		Make/O/D/N=(3) eField = {sin(thetaVal)*cos(phi),sin(thetaVal)*sin(phi),cos(thetaVal)}	
		MatrixOP/O tempMA = (eField^t x filmTensorTrunc x eField)//Check if this is correct!
		tensorInt[i] = tempMA[0]
		
		String pkName = "t_pk"+ num2str(thw[i])
		Make/O/N=2000 $pkName
		Wave z = $pkName
		SetScale/i x,280,360,z
		z = 0
		z = (tensorInt[i]) * gauss(x,pos,wid)
	endfor
	
	KillWaves/Z rotMat0,rotMat90,rotMat180,rotMat270,RotTensor0,RotTensor90,RotTensor180,RotTensor270,TEMPma
	
	DoWindow plotTensorInt
	if(!V_Flag)
		Display/N=plotTensorInt/K=1 tensorInt vs thw
		ModifyGraph mode=4,marker=19,lsize=1.5,grid=2,mirror=1,minor=1
		Label left "Tensor Intensity[a.u.]"
		Label bottom "θ[°]"
		TextBox/C/N=text1/A=RB "\\JCα="+num2str(npw[0])+"° \rModAmp = "+num2str(npw[5])
	else
		TextBox/C/N=text1/A=RB/W=plotTensorInt "\\JCα="+num2str(npw[0])+"° \rModAmp = "+num2str(npw[5])
	endif
End

Function/WAVE truncateMatrix(targetDP,mat)
	
	Variable targetDP
	Wave mat
	
	Variable k,m
	for(k=0;k<=2;k+=1)
		for(m=0;m<=2;m+=1)
			Variable currComp = mat[k][m]
			targetDP = round(targetDP)
			currComp = round(currComp * (10^targetDP)) / (10^targetDP)
			if(currComp == -0)
				currComp = 0
			endif
			mat[k][m] = currComp
		endfor
	endfor
	
	return mat
End

Function calcVectorRatio(v)

	Wave v
	Variable ratio
	//Ratio is defined in terms of z component(v[1]) to xy-component(v[0])
	//What to do when z or xy- component is 0??
	if((v[0] == 0) || (v[1] == 0))
		ratio = 1
	else
		ratio = v[1]/v[0]
	endif
	
	return ratio

End

Function/WAVE truncateVector(targetDP,vec)
	
	Variable targetDP
	Wave vec
	
	Variable n = numpnts(vec)
	Variable k
	for(k=0;k<n;k+=1)
		Variable currComp = vec[k]
		targetDP = round(targetDP)
		currComp = round(currComp * (10^targetDP)) / (10^targetDP)
		if(currComp == -0)
			currComp = 0
		endif
		vec[k] = currComp
	endfor
	return vec
End

Function calcVecNorm(v)

	Wave v
	
	Variable n = numpnts(v),i,mag=0
	for(i=0;i<n;i+=1)
		mag += v[i]^2
	endfor
	mag = sqrt(mag)
	//print mag
	return mag
End

Function/WAVE normVec(v)

	Wave v
	Duplicate/O v,nv
	
	Variable mag = calcVecNorm(v)
	nv = v/mag
	Variable mag2 = calcVecNorm(nv)
	return nv
End

Function/WAVE rotVector(v,deg,axis)

	Wave v
	Variable deg
	String axis
	
	Variable n = DimSize(v,0),ndeg
	Duplicate/O v,rv
	ndeg = deg * (Pi/180)
	if(n == 3)
		if(StringMatch(axis,"Z"))
			Make/O/N=(3,3) rotMat = {{cos(ndeg),sin(ndeg),0},{-sin(ndeg),cos(ndeg),0},{0,0,1}}
		elseif(StringMatch(axis,"X"))
			Make/O/N=(3,3) rotMat = {{1,0,0},{0,cos(ndeg),sin(ndeg)},{0,-sin(ndeg),cos(ndeg)}}
		elseif(StringMatch(axis,"Y"))
			Make/O/N=(3,3) rotMat = {{cos(ndeg),0,-sin(ndeg)},{0,1,0},{sin(ndeg),0,cos(ndeg)}}
		endif
	else
		Make/O/N=(2,2) rotMat = {{cos(ndeg),sin(ndeg)},{-sin(ndeg),cos(ndeg)}}
	endif
	MatrixOP/O rv = rotMat x v 
	Wave rv2 = truncateVector(10,rv)
	Variable mag = calcVecNorm(rv)
	return rv2
End

Function modelTensorvsVector(tensorpw,theta,tensor,alpha,modAmp)

	Wave tensorpw,theta,tensor
	Variable alpha,modAmp
	
	tensorpw[0] = alpha
	tensorpw[5] = modAmp
	modelTensorInt(tensorpw,theta,tensor)
	orbitalWrapper(theta,alpha) 
End

Function/WAVE vecOpsWrap(v,deg,axis)
	Wave v
	Variable deg
	String axis
	
	//Determine magnitude of initial vector and normalize it (creates nv)
	Variable mag = calcVecNorm(v)
	Wave nv = normVec(v)
	Variable mag2 = calcVecNorm(nv)
	//Rotate normalized 2D vector by some degree 
	Wave rv = rotVector(nv,deg,"x")//axis)
	Variable mag3 = calcVecNorm(rv)
	//Calculate the ratio between the rotated vector components --> ModAmp(amp2)
	//print mag,mag2,mag3
	return rv
End

Function calcDeterminant(tensor,rotAngle)

	Wave tensor
	Variable rotAngle
	
	Variable detIni,detFin,As,Bs,Cs,Ds,Es,Fs,Gs,Hs,Is,Ae,Be,Ce,De,Ee,Fe,Ge,He,Ie 
	rotAngle = rotAngle * (Pi/180)
	Make/O/D/N=(3,3) rotMat = {{cos(rotAngle)  ,sin(rotAngle),0}  ,{-sin(rotAngle)  ,cos(rotAngle),0}  ,{0,0,1}}
	MatrixOP/O  newTensor = rotMat x tensor x rotMat^t	
	As = tensor[0][0];Bs = tensor[0][1];Cs = tensor[0][2]
	Ds = tensor[1][0];Es = tensor[1][1];Fs = tensor[1][2]
	Gs = tensor[2][0];Hs = tensor[2][1];Is = tensor[2][2]
	detIni = As * (Es*Is - Fs*Hs) - Bs * (Ds*Is - Fs*Gs) + Cs * (Ds*Hs - Es*Gs)
	Ae = newtensor[0][0];Be = newtensor[0][1];Ce = newtensor[0][2]
	De = newtensor[1][0];Ee = newtensor[1][1];Fe = newtensor[1][2]
	Ge = newtensor[2][0];He = newtensor[2][1];Ie = newtensor[2][2]
	detFin = Ae * (Ee*Ie - Fe*He) - Be * (De*Ie - Fe*Ge) + Ce * (De*He - Ee*Ge)
	print detIni
	print detFin
End

Function calcNormTensor(tensor)
	Wave tensor
	
	Variable i,j
	for(i=0;i<3;i+=1)
		Variable normT = (tensor[i][j])^2 
		if(j==2)
			break
		elseif(i<2)
			j+=1
		else
			i=0
		endif
	endfor
	print normT
End

Function/WAVE truncate2D(w,targetDP)

	Wave w
	Variable targetDP
	
	WaveStats/Q w
	Variable k,m
	//If the value of a component is less than tol times the maximum tensor component then set it to 0 
	for(k=0;k<=2;k+=1)
		for(m=0;m<=2;m+=1)
			Variable currComp = w[k][m]
			targetDP = round(targetDP)
			currComp = round(currComp * (10^targetDP)) / (10^targetDP)
			if(currComp == -0)
				currComp = 0
			endif
			w[k][m] = currComp
		endfor
	endfor	
	return w
End

Function calcRedChiSq(ORI,FIT)
	Wave ORI,FIT
	
	Variable i,n=numpnts(ori),redChiSq=0,num,den
	for(i=0;i<n;i+=1)
		num = (ori[i] - fit[i])^2 
		den = ori[i]^2
		redChiSq += (num/den)
	endfor
	return redChiSq
End

Function/WAVE scaleExptoBA(expSpec,expEn,anchorStep1,anchorStep2,mol)
	
	Wave expSpec,expEn
	Variable anchorStep1,anchorStep2
	String mol
		
	//Find the bare atom waves
	Wave mu_energy = findBAenergy()
	Wave mu = findBAMA(mol)
	
	//Find energies in mu 
	Variable lo = BinarySearch(mu_energy,anchorStep1)
	Variable hi = BinarySearch(mu_energy,anchorStep2)
	WaveStats/R=[lo,hi]/Q mu
	Variable mins = V_min
	Duplicate/O expSpec,test
	
	//Find the experimental energies used for scaling   
	Variable exp280 = findWaveValAtEergy(expEn,expSpec,anchorStep1)
	Variable exp360 = findWaveValAtEergy(expEn,expSpec,anchorStep2)
	
	//Find the bare atom energies used for scaling   
	Variable bareStepLo = findWaveValAtEergy(mu_energy,mu,anchorStep1)
	Variable bareStepHi = findWaveValAtEergy(mu_energy,mu,anchorStep2)
			
	Variable expScale = (exp360-exp280) / (bareStepHi-bareStepLo)
	test = test/ expScale + V_min
	return test
End

Function MAtoF2F1(maw,ew,rho,mol)
	Wave maw,ew	//Mass absorbance and energy waves
	Variable rho //Density of molecule
	String mol//Name of molecule
	
	Variable Na=6.0221415e23 //Avogadro's number [at/d.mol]
	Variable re=2.81794e-13 //classical electron radius [cm]
	Variable hc=1.97263e-05//eV*cm //1.23984e-4//Product of speed of light and Planck's constant
	Wave mu = findBAMA(mol)
	Variable Mw=NumberByKey("Mw",Note(mu))//Molecular weight of molecule
	Variable zstar=NumberByKey("Zstar",Note(mu))//Zstar for molecule
	String betaWaveStr=removeending(nameofwave(maw),"_ma")+"_beta"
	String f2WaveStr=removeending(betaWaveStr,"_beta")+"_f2"
	Duplicate/O maw,$betaWaveStr,$f2WaveStr
	Wave betas = $betaWaveStr
	Wave F2 = $f2WaveStr
	//Convert mass absorbance to beta and then to f2
	betas = maw*hc/(4*pi*ew)
	F2=(Mw*maw*ew)/(2*re*hc*Na)
	//Calculate f1 from f2 wave 
	//Wave f1 = NXA_KK_FFT2(ew, F2, zstar)
	Wave f1 = BAC_KK_FFT(ew, F2, zstar,19)//
	String deltaWaveStr=removeending(nameofwave(f1),"_f1")+"_delta"
	Duplicate/O/D f1,$deltaWaveStr
	Wave delta = $deltaWaveStr
	Variable n=numpnts(ew),i
	for(i=0;i<n;i+=1)
	Variable lambda = hc/ew[i]
		delta[i] = (F1[i]*rho*Na*re*lambda^2)/(2*Pi*Mw)
	endfor
End

Function/WAVE BAC_KK_FFT(EnergyWave, f2Wave, Z_star,pad)
	wave EnergyWave, f2Wave
	variable Z_star//sum of effective number of protons of repeating unit
	//Check Z* calculation
	variable pad
	string f2WaveStr	
	redimension/D EnergyWave, f2Wave//, EnergyWave_ori
	f2WaveStr=removeending(nameofwave(f2Wave),"_f2")+"_F1"
	make/D/O/N=(2^pad) tempF2
	setscale/I x, 0, 30000, "eV", tempf2//30000eV is max energy
	Interpolate2/T=1/I=3/Y=tempf2 EnergyWave, f2Wave
	redimension/C/D tempf2
	Insertpoints 2^pad, 2^pad, tempf2
	Duplicate/d/c/o tempf2, isign
	isign[0,2^pad]=cmplx(0,-1)
	isign[2^pad,2^(pad+1)-1]=cmplx(0,1)
	Variable t0=startMStimer
	If( 0)
		MatrixOp/O outWave=real(IFFT(cmplx(real(FFT(tempf2,0)*isign),0),0))
		SetScale/i x, 0, 60000, outwave
		Interpolate2/T=1/I=3/Y=tempf1/X=EnergyWave outwave
	else
		FFT tempf2
		MatrixOp/O tempf2=tempf2*isign
		redimension/R tempf2
		redimension/C tempf2
		IFFT tempf2
		Redimension/R tempf2
		Interpolate2/T=1/I=3/Y=tempf1/X=EnergyWave tempf2
	endif
	print stopMStimer(t0)/1e6
	tempf1=tempf1*2+Z_star
	duplicate/O tempf1 $(f2WaveStr)
	Wave w = $(f2WaveStr)
	return w
End

Function/WAVE NXA_KK_FFT2(EnergyWave_ori, f2Wave_ori, Z_star)
	wave EnergyWave_ori, f2Wave_ori
	variable Z_star//sum of effective number of protons of repeating unit

	string f2WaveStr
	duplicate/O EnergyWave_ori EnergyWave
	duplicate/O f2Wave_ori f2Wave
	redimension/D EnergyWave, f2Wave, EnergyWave_ori
	f2WaveStr=removeending(nameofwave(f2Wave_ori),"_f2")+"_F1"
	make/D/O/N=(2^21) sampleEn
//1
	sampleEn=x*30000/(2^21-1)
	InsertPoints 0,2, EnergyWave
	EnergyWave[0]=0
	EnergyWave[1]=10-1e-10
	InsertPoints 0,2, f2Wave
	//f2wave[0,1]=f2wave[3]*5
	Interpolate2/T=1/N=2000/I=3/Y=tempf2intp/X=sampleEn EnergyWave, f2Wave
//2
	killwaves sampleEn
	redimension/C/D tempf2intp
	setscale/I x, 0, 30000, "eV", tempf2intp
	FFT/PAD={(2^22)}/DEST=tempf2intpFFT/OUT=1 tempf2intp
//3
	killwaves tempf2intp
	duplicate/O tempf2intpFFT tempf2intpFFTimag
	redimension/R tempf2intpFFTimag
	tempf2intpFFTimag=imag(tempf2intpFFT)
	killwaves tempf2intpFFT
	tempf2intpFFTimag[(2^21), (2^22)]*=-1
	redimension/C tempf2intpFFTimag
	IFFT/DEST=tempf2intpIFFT tempf2intpFFTimag
//4
	killwaves tempf2intpFFTimag
	Redimension/R tempf2intpIFFT
	DeletePoints (2^21),(2^22), tempf2intpIFFT
	make/D/O/N=(2^21) sampleEn
//5
	sampleEn=x*30000/(2^21-1)
	Interpolate2/T=1/N=2000/I=3/Y=tempf1/X=EnergyWave_ori sampleEn, tempf2intpIFFT
//6
	tempf1=tempf1*2+Z_star
	duplicate/O tempf1 $(f2WaveStr)
	Wave w = $(f2WaveStr)
	killwaves tempf2intpIFFT, sampleEn, EnergyWave, f2wave//,tempf1
	return w
End

Function MAtoF2(maw,ew,rho,mol)
	Wave maw,ew	//Mass absorbance and energy waves
	Variable rho //Density of molecule
	String mol//Name of molecule
	
	Variable Na=6.0221415e23 //Avogadro's number [at/d.mol]
	Variable re=2.81794e-13 //classical electron radius [cm]
	Variable hc=1.97263e-05//eV*cm //1.23984e-4//Product of speed of light and Planck's constant
	Wave mu = findBAMA(mol)
	Variable Mw=NumberByKey("Mw",Note(mu))//Molecular weight of molecule
	Variable zstar=NumberByKey("Zstar",Note(mu))//Zstar for molecule
	String betaWaveStr=removeending(nameofwave(maw),"_ma")+"_beta"
	String f2WaveStr=removeending(betaWaveStr,"_beta")+"_f2"
	Duplicate/O maw,$betaWaveStr,$f2WaveStr
	Wave betas = $betaWaveStr
	Wave F2 = $f2WaveStr
	//Convert mass absorbance to beta and then to f2
	betas = maw*hc/(4*pi*ew)
	F2=(Mw*maw*ew)/(2*re*hc*Na)
End

Function makeTensorSpecsOCs(pw,yw,ew,thw,sw)

	Wave pw	//1D Parameter wave. Has peak positions,widths, and amplitudes. has initial guess for alpha and IPs required to build step edge
	Wave yw	//Wave containing the different spectra to be fit
	Wave ew	//Wave containing the energies of the various spectra to be fit (i.e. the x wave)
	Wave thw	//Wave containing the various values of theta
	Wave sw	//Wave containing the concatenated steps from the DFT IPs
		
	Variable i,j=0,k,m
	Variable targetDP = 10
	//Make the absorption tensor for each transition, normalize it, and place each transition into a layer of a 3D wave
	makeTensor1D(pw)
	normalizeTensor()	
	Wave norm3D   = make3DnormWave()
	Wave tensor3D = make3DResonance()
	
	//Add the step edge
	yw = sw
	Duplicate/O yw,output
	//Make waves that will contain tensor element NEXAFS
	Duplicate/O yw,xxSpec,yySpec,zzSpec,xySpec,xzSpec,yzSpec
	//Reference hold Wave. To be used to determine how the peak amplitudes will be defined
	Wave holdWave
	
	Variable nTransitions = DimSize(norm3D,2)	
	Variable pos,wid,amp,alpha,i0,phi,localalpha,amp2	
	Variable p0,pf,nSpec=howManySpec(thw),cSpec=1,thetaVal	
	
	Make/O/N=(3,3,nTransitions) filmTensor3D,r0x3D,r90x3D,r180x3D,r270x3D,tilt3D,correctedMolTensor3D
	
	Variable nAtoms = pw[0]
		
	Make/O/D/N=(3,3) currentTensor = 0 
	
	Variable nPnts = numpnts(yw)/nSpec
	
	Make/O/N=(nTransitions) MA = 0,vecAngle = 0,vecAngle2 = 0
	 
	//Wave that will populate the XX/YY and ZZ tensor elements into 2D vector
	Make/O/N=2 vector
	//Make rotation matrices for rotation of 90,180,and 270 degrees around z-axis
	Make/O/D/N=(3,3) rotMat0   = {{cos(0*(Pi/180))  ,sin(0*(Pi/180)),0}  ,{-sin(0*(Pi/180))  ,cos(0*(Pi/180)),0}  ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat90  = {{cos(90*(Pi/180)) ,sin(90*(Pi/180)),0} ,{-sin(90*(Pi/180)) ,cos(90*(Pi/180)),0} ,{0,0,1}}
	Make/O/D/N=(3,3) rotMat180 = {{cos(180*(Pi/180)),sin(180*(Pi/180)),0},{-sin(180*(Pi/180)),cos(180*(Pi/180)),0},{0,0,1}}
	Make/O/D/N=(3,3) rotMat270 = {{cos(270*(Pi/180)),sin(270*(Pi/180)),0},{-sin(270*(Pi/180)),cos(270*(Pi/180)),0},{0,0,1}}
	
	//Duplicated parameter wave that will contain the results of rotating the molecular tensor
	Duplicate/O pw,pwMolAdj,pwFilm
	Wave pwMolAdj,pwFilm

	i0         = pw[1] 
	alpha      = pw[2] * (Pi/180)
	phi        = pw[3] * (Pi/180)
	Variable count = 0
	for(i=nAtoms + 4 ;i<11*nTransitions + nAtoms;i+=11)
		pos        = pw[i + 0]
		wid        = ABS(pw[i + 1])
		amp        = ABS(pw[i + 2])
		localalpha = pw[i + 9]
		amp2       = pw[i + 10]
		
		currentTensor = tensor3D[p][q][j]// norm3D[p][q][j]
		vector[0] = currentTensor[0][0]
		vector[1] = currentTensor[2][2]
		Wave rv = vecOpsWrap(vector,amp2,"Z")
		currentTensor[0][0] =  abs(rv[0])
		currentTensor[1][1] =  abs(rv[0])
		currentTensor[2][2] =  abs(rv[1])
		
		if(WaveExists(pwMolAdj))
			pwMolAdj[i+3] = currentTensor[0][0]
			pwMolAdj[i+4] = currentTensor[0][1]
			pwMolAdj[i+5] = currentTensor[0][2]
			pwMolAdj[i+6] = currentTensor[1][1]
			pwMolAdj[i+7] = currentTensor[1][2]
			pwMolAdj[i+8] = currentTensor[2][2]
		endif
		
		Variable ang2 = atan(0.5*(tensor3D[2][2][j]/tensor3D[0][0][j])) * (180/pi)//atan(0.5*(rv[1]/rv[0])) * (180/pi)
		Variable ang3 = atan((0.5*rv[1])/rv[0]) * (180/pi)
		vecAngle[j] = ang2
		vecAngle2[j] = ang3
		currentTensor *= i0 * amp// * currentTensor
		Make/O/D/N=(3,3) rotMatAlignX = {{1,0,0},{0,cos(alpha),sin(alpha)},{0,-sin(alpha),cos(alpha)}}
			
		//Tilt the tensor by angle alpha
		MatrixOP/O  tiltedTensorX = rotMatAlignX x currentTensor x rotMatAlignX^t			 
		//Construct the film tensor by adding the azimuthal rotations 
		MatrixOP/O  RotTensor0X   =  rotMat0   x tiltedTensorX x rotMat0^t
		MatrixOP/O  RotTensor90X  =  rotMat90  x tiltedTensorX x rotMat90^t
		MatrixOP/O  RotTensor180X =  rotMat180 x tiltedTensorX x rotMat180^t
		MatrixOP/O  RotTensor270X =  rotMat270 x tiltedTensorX x rotMat270^t
		MatrixOP/O filmTensor = RotTensor0X + RotTensor90X + RotTensor180X + RotTensor270X
		
		
		Wave filmTensorT = truncate2D(filmTensor,targetDP)
		
		if(WaveExists(pwFilm))
			pwFilm[i+3] = filmTensorT[0][0]
			pwFilm[i+4] = filmTensorT[0][1]
			pwFilm[i+5] = filmTensorT[0][2]
			pwFilm[i+6] = filmTensorT[1][1]
			pwFilm[i+7] = filmTensorT[1][2]
			pwFilm[i+8] = filmTensorT[2][2]
		endif
		
		correctedMolTensor3D[][][j] = currentTensor[p][q]	//Molecular tensor after rotating ZZ/XX elements by ModTheta
		filmTensor3D[][][j]         = filmTensorT[p][q]		//Film Tensor after tilting molecular tensor by alpha and 4-fold addition
		r0x3D[][][j]                = RotTensor0X[p][q]
		r90x3D[][][j]               = RotTensor90X[p][q]
		r180x3D[][][j]              = RotTensor180X[p][q]
		r270x3D[][][j]              = RotTensor270X[p][q]
		tilt3D[][][j]               = tiltedTensorX[p][q]
		
		if(j<nTransitions)
			if(cSpec == 1)
				p0 = 0
				pf = nPnts-1
				thetaVal = thw[p0] * (Pi/180)
			elseif(cSpec == nSpec)
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0] * (Pi/180)
			else
				p0 = (cSpec-1)*nPnts
				pf = cSpec*nPnts-1
				thetaVal = thw[p0] * (Pi/180)
			endif
		endif
		
		//Apply electric field to current tensor to obtain the mass absorbance
		Make/O/D/N=(3) eField = {sin(thetaVal)*cos(phi),sin(thetaVal)*sin(phi),cos(thetaVal)}			
		MatrixOP/O tempMA = eField^t x filmTensor x eField
		
		
		if(numtype(tempMA[0]) !=0)
			print "Transition " + num2str(j) + " with a theta value of " + num2str(thetaVal) + " has a problem"
			//MA[j] = 0
		else
			MA[j] = tempMA[0] //abs(tempMA[0])
		endif
		//Make the BB peaks for display
		String pkName = "pk" + num2str(j) + "_spec" + num2str(cSpec)
		Make/O/N=2000 $pkName
		Wave z = $pkName
		SetScale/i x,280,360,z
		z = 0
		z = MA[j] * gauss(x,pos,wid)

		//yw[p0,pf] += sqrt(2*Pi) * wid * MA[j] * gauss(ew,pos,wid)
		output[p0,pf] += MA[j] * gauss(ew,pos,wid)
		xxSpec[p0,pf] += currentTensor[0][0] * gauss(ew,pos,wid)
		yySpec[p0,pf] += currentTensor[1][1] * gauss(ew,pos,wid)
		zzSpec[p0,pf] += currentTensor[2][2] * gauss(ew,pos,wid)
		xySpec[p0,pf] += currentTensor[0][1] * gauss(ew,pos,wid)
		xzSpec[p0,pf] += currentTensor[0][2] * gauss(ew,pos,wid)
		yzSpec[p0,pf] += currentTensor[1][2] * gauss(ew,pos,wid)
		j+=1
		
		if(j<=nTransitions)
			if(cSpec == 1)
				if(j==(nTransitions) && cSpec < nSpec)
					i=nAtoms + 4 - 11
					cSpec+=1
					j=0
				endif 			
			else
				if(j==(nTransitions) && cSpec < nSpec)
					i=nAtoms + 4 - 11
					cSpec+=1
					j=0
				endif 
			endif
		endif
	endfor
	
	KillWaves/Z rotMat0,rotMat90,rotMat180,rotMat270,RotTensor0,RotTensor90,RotTensor180,RotTensor270,TEMPma
	splitSpec(ew,4,"en")
	splitSpec(xxSpec,4,"NXFS_XX")
	splitSpec(yySpec,4,"NXFS_YY")
	splitSpec(zzSpec,4,"NXFS_ZZ")
	splitSpec(xySpec,4,"NXFS_XY")
	splitSpec(xzSpec,4,"NXFS_XZ")
	splitSpec(yzSpec,4,"NXFS_YZ")
	Wave NXFS_XX1,NXFS_XX2,NXFS_XX3,NXFS_XX4,NXFS_YY1,NXFS_YY2,NXFS_YY3,NXFS_YY4,NXFS_ZZ1,NXFS_ZZ2,NXFS_ZZ3,NXFS_ZZ4,en1,en2,en3,en4
	MAtoF2(NXFS_XX1,en2,1.6,"ZnPc")
	MAtoF2(NXFS_XX2,en2,1.6,"ZnPc")
	MAtoF2(NXFS_XX3,en3,1.6,"ZnPc")
	MAtoF2(NXFS_XX4,en4,1.6,"ZnPc")
	MAtoF2(NXFS_YY1,en2,1.6,"ZnPc")
	MAtoF2(NXFS_YY2,en2,1.6,"ZnPc")
	MAtoF2(NXFS_YY3,en3,1.6,"ZnPc")
	MAtoF2(NXFS_YY4,en4,1.6,"ZnPc")
	MAtoF2(NXFS_ZZ1,en2,1.6,"ZnPc")
	MAtoF2(NXFS_ZZ2,en2,1.6,"ZnPc")
	MAtoF2(NXFS_ZZ3,en3,1.6,"ZnPc")
	MAtoF2(NXFS_ZZ4,en4,1.6,"ZnPc")
//	MAtoF2(NXFS_XY1,en2,1.6,"ZnPc")
//	MAtoF2(NXFS_XY2,en2,1.6,"ZnPc")
//	MAtoF2(NXFS_XY3,en3,1.6,"ZnPc")
//	MAtoF2(NXFS_XY4,en4,1.6,"ZnPc")
//	MAtoF2(NXFS_XZ1,en2,1.6,"ZnPc")
//	MAtoF2(NXFS_XZ2,en2,1.6,"ZnPc")
//	MAtoF2(NXFS_XZ3,en3,1.6,"ZnPc")
//	MAtoF2(NXFS_XZ4,en4,1.6,"ZnPc")
//	MAtoF2(NXFS_YZ1,en2,1.6,"ZnPc")
//	MAtoF2(NXFS_YZ2,en2,1.6,"ZnPc")
//	MAtoF2(NXFS_YZ3,en3,1.6,"ZnPc")
//	MAtoF2(NXFS_YZ4,en4,1.6,"ZnPc")
End