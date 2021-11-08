#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function seqThresholds(atomName,fnum,osWave,ovpWave,maxRot,LUMOwave,[broadShift,erange,broad1,broad2,Eini,Efin,d,gRes,buildSym,corrWave,dClusters,realign,thx,thy,thz,rotOrder,modelAlpha,IPwave,expSpecName,expEnergyName,expFolderPath,mol,alpha,i0,phi,fit,rigidShift,stepShift,thetaList,NEXAFStype,startPre,startPost,endPre,endPost,justModel,justFit,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,holdPos,holdWidths,holdAmps,refinement,maskEnergy1,maskEnergy2,pkToRefine,refitPWave,holdTensorElems])
	
	String atomName
	Variable fnum
	Wave osWave,ovpWave
	Variable maxRot	//What's the symmetry of the molecule? Currently supported (1 = no symmetry, 2 = 2-fold rotation [180 degrees], 3 = 3 fold rotation [120 degrees], 4 = 4-fold rotation [90 degrees])
	Wave LUMOwave
	//Optional Parameters for Clustering/Filtering Algorithm
	Variable broadShift //Allows for the adjustment of the energies at which to apply the broadening scheme
	Variable erange //Maximum energy range to consider transitions in
	Variable broad1 //Peak width at lower energy (ewid1)
	Variable broad2 //Peak width at higher energy (ewid2)
	Variable Eini	//Starting energy at which to make transition peaks
	Variable Efin	//Ending energy at which to make transition peaks
	
	Variable d //Display the data?
	Variable gRes	//Energy Resolution of the Gaussian peaks 
	String buildSym	//Calculate TDM tensor and determine theta,symmetry type and x,y,z components of OS? Yes or no
	Wave corrWave	////Energy correction wave. Must be provided if every atom has different energy corrections. If all atoms have the same energy corrections then, rigidShift is enough.
	Variable dClusters	//Display cluster plots from symmetry/overlap routines?
	Variable realign
	Variable thx
	Variable thy
	Variable thz
	String rotOrder
	//Fitting of DFT optical model to experimental NEXAFS
	String modelAlpha
	Wave IPwave
	String expSpecName
	String expEnergyName	
	String expFolderPath
	String mol
	
	Variable startPre 		//At what energy should the preedge start to be defined?
	Variable endPre 			//At what energy should the preedge stop to be defined?
	Variable startPost 		//At what energy should the postedge start to be defined?
	Variable endPost			//At what energy should the postedge stop to be defined?
	//Optional Parameters for fitting of DFT to experiment
	Variable alpha
	Variable i0
	Variable phi
	String fit
	Variable rigidShift
	Variable stepShift
	String thetaList
	String NEXAFStype
	
	Variable justModel
	Variable justFit
	Variable anchorStep1
	Variable anchorStep2
	Variable anchorExp1
	Variable anchorExp2
	Variable stepWid1,stepWid2,stepE1,stepE2
	
	Variable holdPos,holdWidths,holdAmps,holdTensorElems
	Variable refinement
	Variable maskEnergy1
	Variable maskEnergy2
	Variable pkToRefine
	Wave refitPWave
	
	Variable seqTimer = StartMSTimer
	
	if(ParamIsDefault(d))
		d = 0	//Don't make cluster plots
	endif
	
	if(ParamIsDefault(buildSym))
		buildSym = "no"	//Take into account molecule symmetry?
	endif
	
	String iniSequenceFolder = GetDataFolder(1)
	
	Variable nOS = numpnts(osWave)
	Variable nOVP = numpnts(ovpWave)
		
	Variable i,j
	Variable currentOS, currentOVP
	
	Make/O/N=(nOS,nOVP) pDiff,rChiSq,nPeaks
	for(i=0;i<=nOS-1;i+=1)
		currentOS = osWave[i]
		for(j=0;j<=nOVP-1;j+=1)
			currentOVP = ovpWave[j] 
			filterDFT(atomName,fnum,currentOS,currentOVP,maxRot,LUMOwave,broadShift=broadShift,erange=erange,broad1=broad1,broad2=broad2,Eini=Eini,Efin=Efin,d=d,gRes=gRes,buildSym=buildSym,corrWave=corrWave,dClusters=dClusters,realign=realign,thx=thx,thy=thy,thz=thz,rotOrder=rotOrder,modelAlpha=modelAlpha,IPwave=IPwave,expSpecName=expSpecName,expEnergyName=expEnergyName,expFolderPath=expFolderPath,mol=mol,alpha=alpha,i0=i0,phi=phi,fit=fit,rigidShift=rigidShift,stepShift=stepShift,thetaList=thetaList,NEXAFStype=NEXAFStype,startPre=startPre,startPost=startPost,endPre=endPre,endPost=endPost,justModel=justModel,justFit=justFit,anchorStep1=anchorStep1,anchorStep2=anchorStep2,anchorExp1=anchorExp1,anchorExp2=anchorExp2,stepWid1=stepWid1,stepWid2=stepWid2,stepE1=stepE1,stepE2=stepE2,holdPos=holdPos,holdWidths=holdWidths,holdAmps=holdAmps,refinement=refinement,maskEnergy1=maskEnergy1,maskEnergy2=maskEnergy2,pkToRefine=pkToRefine,refitPWave=refitPWave,holdTensorElems=holdTensorElems)
			NVAR perDiff = root:perDiff
			NVAR numpeaks = root:numpeaks
			NVAR reducedChiSq = root:redchiSq
			pDiff[i][j]    = perDiff
			nPeaks[i][j]   = numpeaks
			rChiSq[i][j] = reducedChiSq		
		endfor
	endfor
	
	SetDataFolder iniSequenceFolder
	//Calculate the Bayesian Inference Coefficient
	String pathToExpWave = expFolderPath + "allExpSpec"
	Wave expSpecWave = $pathToExpWave
	Variable n = numpnts(expSpecWave)
	Wave BICwave = calculateBIC(rChiSq,nPeaks,3,n)
	plotParamSpaceResults(osWave,ovpWave,BICwave,pDiff,nPeaks,rChiSq)
	Variable endTimer = StopMSTimer(seqTimer)/1000000
	
	if(endTimer <= 60)
		print "Process has ended after " + num2str(endTimer) + " seconds."
	elseif(endTimer > 60)
		print "Process has ended after " + num2str(endTimer/60) + " minutes."
	endif
End

Function/WAVE calculateBIC(chSqWave,nPeaksWave,nParams,n)

	Wave chSqWave //2d wave containing chiSq values from p-space exploration
	Wave nPeaksWave //2d wave containing the number of BB peaks resulting from each OS/OVP combination
	Variable nParams //How many parameters were opened during the fit?
	Variable N //How many points does the experimental data have?
	
	Variable x = DimSize(chSqWave,0),y=DimSize(chSqWave,1),i,j
	Make/O/N=(x,y) BICWave
	for(i=0;i<x;i+=1)
		for(j=0;j<y;j+=1)
			Variable GF =  chSqWave[i][j]
			Variable k = nPeaksWave[i][j] * nParams
			Variable BIC = (n - k) *(GF/n) + k * ln(n)
			BICWave[i][j] = BIC
		endfor
	endfor
	
	return BICWave
End

Function plotParamSpaceResults(osw,ovpw,BICwave,pDiff,nPeaks,rChiSq)
	Wave osw,ovpw,BICwave,pDiff,nPeaks,rChiSq
	
	Variable nx = numpnts(osw),ny=numpnts(ovpw),i=0
	//Make axis waves for plotting image
	Make/O/N=(nx)/T ost
	Make/O/N=(ny)/T ovpt
	
	Make/O/N=(nx) osScale
	Make/O/N=(ny) ovpScale

	for(i=0;i<nx;i+=1)
		ost[i] = num2str(osw[i])
		osScale[i] = i
	endfor
	
	for(i=0;i<ny;i+=1)
		ovpt[i] = num2str(ovpw[i])
		ovpScale[i] = i
	endfor
	
	//Make plots
	DoWindow BICPlot
	if(!V_Flag)
		Display/N=BICPlot/K=1/W=(0,0,300,500)  
		AppendImage/W=BICPlot BICWave
		ModifyGraph userticks(left)={ovpScale,ovpt},userticks(bottom)={osScale,ost}
		SetAxis left ny-0.5,-0.5
		SetAxis bottom -0.5,nx-0.5
		ModifyGraph mirror=1,fStyle=1,margin(right)=86
		Label left "OVP[%]"
		Label bottom "OS[%]"
		ModifyImage BICWave ctab= {*,*,YellowHot256,1}
		ColorScale/C/N=text0/A=RC/X=1.00/Y=5.00/E image=BICWave,fstyle=1
		ColorScale/C/N=text0 "BIC"
	endif
	DoWindow nPeaksPlot
	if(!V_Flag)
		Display/N=nPeaksPlot/K=1/W=(300,0,600,500) 
		AppendImage/W=nPeaksPlot nPeaks
		ModifyGraph userticks(left)={ovpScale,ovpt},userticks(bottom)={osScale,ost}
		SetAxis left ny-0.5,-0.5
		SetAxis bottom -0.5,nx-0.5
		ModifyGraph mirror=1,fStyle=1,margin(right)=86
		Label left "OVP[%]"
		Label bottom "OS[%]"
		ModifyImage nPeaks ctab= {*,*,YellowHot256,1}
		ColorScale/C/N=text1/A=RC/X=1.00/Y=5.00/E image=nPeaks,fstyle=1
		ColorScale/C/N=text1 "nPeaks"
	endif
	
	DoWindow pDiffPlot
	if(!V_Flag)
		Display/N=pDiffPlot/K=1/W=(600,0,900,500)  
		AppendImage/W=pDiffPlot pDiff
		ModifyGraph userticks(left)={ovpScale,ovpt},userticks(bottom)={osScale,ost}
		SetAxis left ny-0.5,-0.5
		SetAxis bottom -0.5,nx-0.5
		ModifyGraph mirror=1,fStyle=1,margin(right)=86
		Label left "OVP[%]"
		Label bottom "OS[%]"
		ModifyImage pDiff ctab= {0,10,YellowHot256,1}
		ColorScale/C/N=text2/A=RC/X=1.00/Y=5.00/E image=pDiff,fstyle=1
		ColorScale/C/N=text2 "%Diff"
	endif
	
	DoWindow nChiSqPlot
	if(!V_Flag)
		Display/N=nChiSqPlot/K=1/W=(900,0,1200,500)  
		AppendImage/W=nChiSqPlot rChiSq
		ModifyGraph userticks(left)={ovpScale,ovpt},userticks(bottom)={osScale,ost}
		SetAxis left ny-0.5,-0.5
		SetAxis bottom -0.5,nx-0.5
		ModifyGraph mirror=1,fStyle=1,margin(right)=86
		Label left "OVP[%]"
		Label bottom "OS[%]"
		ModifyImage rChiSq ctab= {*,*,YellowHot256,1}
		ColorScale/C/N=text3/A=RC/X=1.00/Y=5.00/E image=rChiSq,fstyle=1
		ColorScale/C/N=text3 "rChiSq"
	endif
End