#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function symmetryOps(TDMx,TDMy,TDMz,E,tol,muWave0,normTDMs0,symmetryParams,dotProdNew,thetaNew,newVecNorm,sumRotOSW,oriOSW,oriOSadjW,targetDP,zVec,maxRot)
	Wave TDMx,TDMy,TDMz,E,muWave0,normTDMs0,symmetryParams,dotProdNew,thetaNew,newVecNorm,sumRotOSW,oriOSW,oriOSadjW,zVec
	Variable tol,targetDP,maxRot
	
	Duplicate/O 	TDMx,originalSymmetries
	
	//Make rotation matrices around z-axis and rotated TDMs
	makeRotMats(maxRot,muWave0,normTDMs0)
	
	Variable i,j=0,k,m,x=DimSize(TDMx,0)
	//Make TDM tensor for molecular frame for each transition at each equivalent symmetry center (0,90,180,270 degrees TDM) and convert tensor back into vector via diagonal components
	for(i=0;i<=x-1;i+=1)
		
		makeTDMTensor(j,i,maxRot)
		String tdmTensorName = "tdmTensor0_" + num2str(i)
		String totTensorName = "totTensor_"    + num2str(i)
		Wave wTot = $totTensorName , w0 = $tdmTensorName 
		
		//Truncate tensor elements to prevent floating/rounding problems
		truncateSym(wTot,targetDP,tol)
		Wave wTot = $totTensorName
		Duplicate/O wTot, wTot2	
		
		//Check transition symmetry: Isotropic = 0 , Uniaxial = 1 , Biaxial = 2 , Triaxial = 3
		//Determine whether tensor is isotropic, uniaxial, biaxial, or triaxial
		determineTensorSym(wTot2,symmetryParams,originalSymmetries,i,0)
		
		//Truncate tensor elements to prevent floating/rounding problems
		truncateSym(w0,targetDP,tol)
		Wave w0 = $tdmTensorName 
		
		//Determine whether original tensor is isotropic, uniaxial, biaxial, or triaxial
		determineTensorSym(w0,symmetryParams,originalSymmetries,i,1)
		
		//Populate the parameter wave containing the symmetry information
		populateSymPWave(symmetryParams,wTot2,E,i,TDMx,TDMy,TDMz,sumRotOSW,oriOSW,oriOSadjW,wTot,dotProdNew,thetaNew,newVecNorm,zVec)
		
		killRotationWaves(maxRot,i)
	endfor
	make3Dtensor()	
End

Function make3Dtensor()
	String tensorList  = WaveList("totTensor_*",";","")
	String tensorList2 = WaveList("tdmTensor0_*",";","")
	Variable nTensors = ItemsInList(tensorList)
	
	Wave allTensors,allIniTensors
	if(WaveExists(allTensors))
		KillWaves/Z allTensors
	endif
	
	if(WaveExists(allIniTensors))
		KillWaves/Z allIniTensors
	endif
	
	Make/O/N=(3,3,nTensors) allTensors,allIniTensors
	Concatenate/O/NP=2 tensorList,allTensors
	Concatenate/O/NP=2 tensorList2,allIniTensors
	
	Variable i
	for(i=0;i<=nTensors-1;i+=1)
		String  tenName = StringFromList(i,tensorList)
		String  tenName2 = StringFromList(i,tensorList2)
		Wave w = $tenName
		Wave x = $tenName2
		KillWaves/Z w,x
	endfor
End

//Make rotation matrices for rotation around z-axis and normalized rotated tensors
Function makeRotMats(maxRotation,muWave0,normTDMs0)
	
	Variable maxRotation
	Wave muWave0,normTDMs0
	
	Variable rotDegrees = 360/maxRotation
	Variable i
	for(i=1;i<maxRotation;i+=1)
		Variable cRot =  rotDegrees*i
		String rotMatName = "rotMat" + num2str(cRot)
		Wave rotMat = $rotMatName
		Make/O/N=(3,3) rotMat   = {{cos(cRot*(Pi/180)),sin(cRot*(Pi/180)),0},{-sin(cRot*(Pi/180)),cos(cRot*(Pi/180)),0},{0,0,1}}
		String normTDMName = "normTDMs" + num2str(cRot)
		Duplicate/O muWave0,$normTDMName
		Wave normTDM = $normTDMName
		MatrixOP/O normTDM = normTDMs0 x rotMat 
	endfor
	
End

Function makeTDMTensor(j,i,maxRot)
	
	Variable j,i,maxRot
	
	Variable rotDegrees = 360/maxRot
	
	String totTensor     = "totTensor_"    + num2str(i)
	Make/O/D/N=(3,3) $totTensor
	Wave wTot = $totTensor
	wTot = 0
	Variable k=0
	//Add up the TDM tensors for each symmetry equivalent center
	for(k=0;k<maxRot;k+=1)
		Variable cRot =  rotDegrees*k
		String tdmTensorName = "tdmTensor" + num2str(cRot) + "_" + num2str(i)
		String normTDMName   = "normTDMs" + num2str(cRot)
		Make/O/D/N=(3,3) $tdmTensorName
		Wave w = $tdmTensorName
		Wave normTDM = $normTDMName

		w[j][0]   = normTDM[i][0]^2             ;w[j][1]   = normTDM[i][0]*normTDM[i][1] ;w[j][2]     = normTDM[i][0]*normTDM[i][2]
		w[j+1][0] = normTDM[i][1]*normTDM[i][0] ;w[j+1][1] = normTDM[i][1]^2             ;w[j+1][2]   = normTDM[i][1]*normTDM[i][2]
		w[j+2][0] = normTDM[i][2]*normTDM[i][0] ;w[j+2][1] = normTDM[i][2]*normTDM[i][1] ;w[j+2][2]   = normTDM[i][2]^2 
		wTot += w
	endfor	
End

Function truncateSym(w,targetDP,tol)

	Wave w
	Variable targetDP,tol
	
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
		//	elseif(currComp < V_max*tol)
		//		currComp = 0 //If current tensor element is less than product of  tolerance and  max tensor element then set it to 0
			endif
			w[k][m] = currComp
		endfor
	endfor	
End

Function determineTensorSym(w,symmetryParams,originalSymmetries,i,iniOrFin)
	
	Wave w,symmetryParams,originalSymmetries
	Variable i,iniOrFin
	
	Variable xx,xy,xz,yx,yy,yz,zx,zy,zz
	xx = w[0][0] ; xy = w[0][1] ; xz = w[0][2] 	
	yx = w[1][0] ; yy = w[1][1] ; yz = w[1][2]  	
	zx = w[2][0] ; zy = w[2][1] ; zz = w[2][2]
		
	if(iniOrFin == 0)
		if(xy == 0 || xz == 0 || yz == 0 || yz == 0 || zx == 0 || zy == 0)
			if(xx==yy && xx==zz)     									    									//Isotropic Case
				symmetryParams[i][4] = 0
			elseif(((xx==yy) && (xx!=zz)) || ((xx==zz) && (xx!=yy)) || ((yy==zz) && (yy!=xx)))	//Uniaxial case
				symmetryParams[i][4] = 1 
			elseif((xx!=yy) && (yy!=zz) && (xx!=zz))														//Biaxial Case
				symmetryParams[i][4] = 2
			endif
		else  						//If the transition has nonzero offdiagonal elements then its triaxial		    								
			symmetryParams[i][4] = 3														
		endif
	elseif(iniOrFin == 1)
		if(xy == 0 || xz == 0 || yz == 0 || yz == 0 || zx == 0 || zy == 0)
			if(xx==yy && xx==zz)     									    									//Isotropic Case
				originalSymmetries[i] = 0
			elseif(((xx==yy) && (xx!=zz)) || ((xx==zz) && (xx!=yy)) || ((yy==zz) && (yy!=xx)))	//Uniaxial case
				originalSymmetries[i] = 1 
			elseif((xx!=yy) && (yy!=zz) && (xx!=zz))														//Biaxial Case
				originalSymmetries[i] = 2
			endif
		else  						//If the transition has nonzero offdiagonal elements then its triaxial		    								
			originalSymmetries[i] = 3														
		endif
	endif
End

Function populateSymPWave(symmetryParams,w,E,i,TDMx,TDMy,TDMz,sumRotOSW,oriOSW,oriOSadjW,wTot,dotProdNew,thetaNew,newVecNorm,zVec)
	
	Wave symmetryParams,w,E,TDMx,TDMy,TDMz,sumRotOSW,oriOSW,oriOSadjW,wTot,dotProdNew,thetaNew,newVecNorm,zVec
	Variable i
	
	//Place diagonal components of summed tensor into 2D wave
	symmetryParams[i][0] = wTot[0][0]	//TDMx after rotation and summation
	symmetryParams[i][1] = wTot[1][1]   //TDMy after rotation and summation
	symmetryParams[i][2] = wTot[2][2]   //TDMz after rotation and summation
		
	//Normalize the summed up TDM
	newVecNorm[i]  = sqrt(symmetryParams[i][0]^2 + symmetryParams[i][1]^2 + symmetryParams[i][2]^2)
	symmetryParams[i][0] = symmetryParams[i][0]/newVecNorm[i]
	symmetryParams[i][1] = symmetryParams[i][1]/newVecNorm[i]
	symmetryParams[i][2] = symmetryParams[i][2]/newVecNorm[i]
		
	dotProdNew[i] = zVec[0]*symmetryParams[i][0] + zVec[1]*symmetryParams[i][1] + zVec[2]*symmetryParams[i][2]
	thetaNew[i]   = acos(dotProdNew[i])*(180/Pi)
		
	symmetryParams[i][3] = thetaNew[i] //theta
		
	Variable xx,xy,xz,yx,yy,yz,zx,zy,zz
	xx = w[0][0] ; xy = w[0][1] ; xz = w[0][2] 	
	yx = w[1][0] ; yy = w[1][1] ; yz = w[1][2]  	
	zx = w[2][0] ; zy = w[2][1] ; zz = w[2][2]
	
	//Compute the oscillator strength coming from each component AFTER adding TDM together
	Variable osXX,osYY,osZZ,osXY,osXZ,osYZ, en ,oriOS, sumRotOS
	
	//Need to convert energy from eV to Hartree. 1 Hartree = 27.2114 eV
	en   = E[i]/27.2114
	osXX = (2/3)*en*xx^2
	osYY = (2/3)*en*yy^2
	osZZ = (2/3)*en*zz^2
	
	//Need to check if this is right...
	osXY = (2/3)*en*xy^2
	osXZ = (2/3)*en*xz^2
	osYZ = (2/3)*en*yz^2

	sumRotOS = osXX + osYY + osZZ
	oriOS = (2/3)*en*(TDMx[i]^2 + TDMy[i]^2 + TDMz[i]^2)
	Variable ratio = sumRotOS/oriOS
			
	symmetryParams[i][5] = osXX/ratio
	symmetryParams[i][6] = osYY/ratio
	symmetryParams[i][7] = osZZ/ratio
		
	symmetryParams[i][8]  = osXY/ratio
	symmetryParams[i][9]  = osXZ/ratio
	symmetryParams[i][10] = osYZ/ratio
		
	symmetryParams[i][11] = ratio
		
	sumRotOSW[i] = sumRotOS
	oriOSW[i]    = oriOS
	oriOSadjW[i] = symmetryParams[i][5] + symmetryParams[i][6] + symmetryParams[i][7]
End

Function killRotationWaves(maxRot,j)
	
	Variable maxRot,j
	
	Variable rotDegrees = 360/maxRot,i
	for(i=1;i<maxRot;i+=1)
		Variable cRot =  rotDegrees*i
		String tdmTensorName = "tdmTensor" + num2str(cRot) + "_" + num2str(j)
		//String normTDMName   = "normTDMs" + num2str(cRot)
		KillWaves/Z $tdmTensorName
	endfor
End

Function tdmSymmetry2(TDMx,TDMy,TDMz,E,maxRot,[realign,thx,thy,thz,rotOrder])
	
	Wave TDMx,TDMy,TDMz,E
	Variable maxRot
	Variable realign,thx,thy,thz
	String rotOrder				//The code is setup to calculate polar angles w.r.t. z-axis {0 0 1}. If molecular orientation is not coincident with
										// this then the results will be misinterpreted. Use these parameters to specify whether the TDM orientation should be realigned.
										//thx, thy, thz represent the angle in degrees to rotate about each axis, rotOrder is a list that tells the code the order of rotations
	if(ParamIsDefault(realign))
		realign = 0
	endif
	Variable tol = 1/100	//What is the minimum threshold of values to consider for the TDM tensor? If component is less than the V_max*tol then make it 0
	Variable targetDP = 10
	//Make 2D wave containing x,y, and z components of TDM
	Duplicate/O TDMx,TDMx_Ori
	Duplicate/O TDMy,TDMy_Ori
	Duplicate/O TDMz,TDMz_Ori
	Concatenate/O {TDMx,TDMy,TDMz}, muWave0
	Concatenate/O {TDMx,TDMy,TDMz}, muWave0_Ori
	
	//Define principal symmetry axis
	Make/O/N=3 zVec
	zVec = {0,0,1}
	
	//Align P3HT TDMs onto desired coordinate frame (i.e. thiophene plane orthogonal to z-axis)
	if(realign==1)
		Make/O/N=(3,3) rotMatAlignZ = {{cos(thz*(Pi/180)),sin(thz*(Pi/180)),0},{-sin(thz*(Pi/180)),cos(thz*(Pi/180)),0},{0,0,1}}
		Make/O/N=(3,3) rotMatAlignY = {{cos(thy*(Pi/180)),0,-sin(thy*(Pi/180))},{0,1,0},{sin(thy*(Pi/180)),0,cos(thy*(Pi/180))}}
		Make/O/N=(3,3) rotMatAlignX = {{1,0,0},{0,cos(thx*(Pi/180)),sin(thx*(Pi/180))},{0,-sin(thx*(Pi/180)),cos(thx*(Pi/180))}}
		if(StringMatch(rotOrder,"x,y,z"))
			MatrixOP/O muWave0 = muWave0 x rotMatAlignX x rotMatAlignY x rotMatAlignZ 
		elseif(StringMatch(rotOrder,"x,z,y"))
			MatrixOP/O muWave0 = muWave0 x rotMatAlignX x rotMatAlignZ x rotMatAlignY 
		elseif(StringMatch(rotOrder,"y,x,z"))
			MatrixOP/O muWave0 = muWave0 x rotMatAlignY x rotMatAlignX x rotMatAlignZ 
		elseif(StringMatch(rotOrder,"y,z,x"))
			MatrixOP/O muWave0 = muWave0 x rotMatAlignY x rotMatAlignZ x rotMatAlignX 
		elseif(StringMatch(rotOrder,"z,x,y"))
			MatrixOP/O muWave0 = muWave0 x rotMatAlignZ x rotMatAlignX x rotMatAlignY 
		elseif(StringMatch(rotOrder,"z,y,x"))
			MatrixOP/O muWave0 = muWave0 x rotMatAlignZ x rotMatAlignY x rotMatAlignX 
		else
			print "Incorrect input for rotation order of rotation"
			return 0 
		endif
	endif
	 
	//Make TDMs into unit vectors
	Duplicate/O TDMx, vecNorm, newVecNorm
	Duplicate/O TDMx,sumRotOSW,oriOSW,oriOSadjW
	Duplicate/O muWave0,normTDMs0
	Variable i,j=0,k,m,x=DimSize(TDMx,0)
	Make/O/N=(x,12) symmetryParams
	Make/O/N=(x) dotProdNew,thetaNew,unitVecNorm
	
	//Calculate norm of each TDM and convert each TDM to a unit vector
	for(i=0;i<=x-1;i+=1)
		vecNorm[i] = sqrt(TDMx[i]^2 + TDMy[i]^2 + TDMz[i]^2)
		normTDMs0[i][0] = muWave0[i][0]/vecNorm[i]
		normTDMs0[i][1] = muWave0[i][1]/vecNorm[i]
		normTDMs0[i][2] = muWave0[i][2]/vecNorm[i]
		unitVecNorm[i] = sqrt(normTDMs0[i][0]^2 + normTDMs0[i][1]^2 + normTDMs0[i][2]^2)
	endfor
	
	//Determine transition symmetries and TDM thetas by taking into account the symmetry of the molecule.	
	symmetryOps(TDMx,TDMy,TDMz,E,tol,muWave0,normTDMs0,symmetryParams,dotProdNew,thetaNew,newVecNorm,sumRotOSW,oriOSW,oriOSadjW,targetDP,zVec,maxRot)
End

Function prepThetaSym2(fnum,iniFolder,atomName,maxRot,[realign,thx,thy,thz,rotOrder])
	Variable fnum,maxRot
	String iniFolder,atomName
	Variable realign,thx,thy,thz
	String rotOrder
	
	if(ParamIsDefault(realign))
		realign = 0
	endif
	
	Variable i,totSymTime=0
	for(i=1;i<=fnum;i+=1)
		Variable symTimer = StartMSTimer
		String cAtom = atomName + num2str(i)
		String cFolder = iniFolder + cAtom
		SetDataFolder $cFolder
		Wave TDMx_,TDMy_,TDMz_,eV_
		if(realign)
			print "Realigning transition dipole moment frame of reference"
			tdmSymmetry2(TDMx_,TDMy_,TDMz_,eV_,maxRot,realign=realign,thx=thx,thy=thy,thz=thz,rotOrder=rotOrder)
		elseif(!realign)
			//print "Calculating symmetry"
			tdmSymmetry2(TDMx_,TDMy_,TDMz_,eV_,maxRot)
		endif
		
		//Make histogram comparing transition symmetries before and after azimuthal rotation
		Wave originalSymmetries,symmetryParams
		String finSymName = cFolder + ":finalSymmetries"
		Duplicate/O originalSymmetries,$finSymName
		Wave finalSymmetries = $finSymName
		finalSymmetries = symmetryParams[x][4]
		Make/N=8/O originalSymmetries_Hist,finalSymmetries_Hist
		Histogram/N/B=2 originalSymmetries,originalSymmetries_Hist
		Histogram/N/B=2 finalSymmetries,finalSymmetries_Hist
		String histName = "SymmetryHistogram_" + cAtom
		DoWindow $histName
		if(!V_Flag)
			Display/N=$histName originalSymmetries_Hist,finalSymmetries_Hist
			SetAxis bottom -0.1,4.1
			ModifyGraph mode=5,lsize=2,rgb(originalSymmetries_Hist)=(0,0,0),tick=2,mirror=1,fSize=14
			Label left "Number of Transitions";DelayUpdate
			Label bottom "Transition Symmetry";DelayUpdate
			Legend/C/N=text0/A=LT ""+cAtom+"\r\\s(finalSymmetries_Hist) Final\r\\s(originalSymmetries_Hist) Initial"
		endif	
		Variable symTimerEnd = stopMSTimer(symTimer)/1000000
		totSymTime += symTimerEnd
		print "Symmetry operations for " + cAtom + " have been completed after " + num2str(symTimerEnd) + " seconds."
	endfor
	
	symmetryLayout(fnum)
	
	print "All symmetry operations completed after " + num2str(totSymTime) + " seconds."
	SetDataFolder iniFolder
End

function truncate(inValue,targetDP)
    Wave inValue
    Variable targetDP
    targetDP = round(targetDP)
    inValue = round(inValue * (10^targetDP)) / (10^targetDP)
end

Function symmetryLayout(nAtoms)

	Variable nAtoms
	
	DoWindow Symmetry_Summary
	if(V_Flag)
		KillWindow/Z Symmetry_Summary
	endif
	
	NewLayout/N=Symmetry_Summary/W=(0,0,410,460)
	Variable nPages = round(nAtoms/8)
	Variable i,j=1
	for(i=1;i<nPages;i+=1)
		LayoutPageAction appendPage
	endfor
	
	String compareSymplots = SortList(WinList("SymmetryHistogram_*",",",""),",",16)
	String arrange = "Tile/O=1"
	for(i=0;i<=nAtoms-1;i+=1)	
		String cSymPlot = StringFromList(i,compareSymplots,",")
		
		if(i<=7)
			AppendLayoutObject/F=1/PAGE=1 graph $cSymPlot
			LayoutPageAction PAGE=1
			Execute/Q/Z arrange
		elseif(i>7)
			if(mod(i,8)==0)
				j+=1
			endif
			LayoutPageAction PAGE=(j)
			AppendLayoutObject/F=1/PAGE=(j) graph $cSymPlot
			Execute/Q/Z arrange
		endif
		DoWindow/HIDE=1 $cSymPlot
	endfor	
End
