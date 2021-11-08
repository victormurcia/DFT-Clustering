#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <KBColorizeTraces>

Function clusteringTransitions(pwave,ovpwave,os,ovp,Eini,Efin,sym,cMethod,mol)
	
	Wave pwave,ovpWave
	Variable os
	Variable ovp
	Variable Eini,Efin
	String sym,mol
	Variable cMethod
		
	Variable nTrans=DimSize(pWave,0),r

	String tpwName = "TotalPWave_" + sym
	Make/O/N=(nTrans,16) $tpwName
	Wave TotalPWave = $tpwName
	TotalPWave = pWave
	MDsort2(pWave,0)	//Sort TotalpWave by increasing Gauss position
	Wave oriNEXAFS = getOriginalDFTNEXAFS(os,ovp,sym,mol)
	Wave filNEXAFS = getFilteredDFTNEXAFS(os,ovp,sym,mol)
	//Add up all the peaks "gTotal"
	Wave gTotal = makeTotalWave(Eini,Efin,TotalPWave,sym)
	String gTotName = "gTotal_" + sym
	Wave gTtotal = $gTotName 
	//Calculate overlap between peaks and merge overlapping peaks
	Variable ncl = 1	//Counter for clustering loops
	peakOverlap(pWave,sym,ncl = ncl)
	String ovpwName = "ovpWave"  + sym  + num2str(ncl)
	String acpwName = "combClusterPW" + sym + "_" + num2str(ncl)
	Wave ovpWave = $ovpwName
	
	if(cMethod == 1)
		simpleCluster2(ovpWave,pWave,ovp,sym,ncl,Eini,Efin)
	elseif(cMethod == 2)
		simpleCluster3(ovpWave,pWave,ovp,sym,ncl,Eini,Efin)
	endif
	Wave combClusterPW = $acpwName
	//Check overlap on clustered peaks
	peakOverlap(combClusterPW,sym,ncl = ncl)
	Wave ovpWave = $ovpwName
	 
	//Continue merging transition clusters until the overlap matrix has no values above OT
	Variable i=0,done = determineMaxOVP(ovpWave,ovp)

	print "Clustering transitions with " + sym + " symmetry"
	if(done == 0)
		do
			acpwName = "combClusterPW" + sym + "_" + num2str(ncl)
			Wave combClusterPW = $acpwName
			ncl+=1
			
			if(cMethod == 1)
				simpleCluster2(ovpWave,combClusterPW,ovp,sym,ncl,Eini,Efin)
			elseif(cMethod == 2)
				simpleCluster3(ovpWave,combClusterPW,ovp,sym,ncl,Eini,Efin)
			endif
			acpwName = "combClusterPW" + sym + "_" + num2str(ncl)
			Wave combClusterPW = $acpwName
	//5/18/2021 ******* Try not fitting amplitudes at each iteration
	//		Variable chiSq = fitAmplitudes3(combClusterPW,oriNEXAFS,280,320,sym,os,ovp)
	//		Wave combClusterPW = $acpwName
	//*********
			Wave ovpWave = peakOverlap(combClusterPW,sym,ncl = ncl)
			done = determineMaxOVP(ovpWave,ovp)
			if(ncl >= 20)//Break clustering loop if 10 iterations are exceeded
				break
			endif
		while(done == 0)
	endif
	
	print "Clustering completed for transitions with " + sym + " symmetry after " + num2str(ncl) + " iterations."
	
	//Make peaks from various trials
	makePeaks(combClusterPW,Eini,Efin,sym)
	Variable nPeaksRand = DimSize(combClusterPW,0)

	//Determine which transitions to cluster together
	removeOldPW(sym,ncl)
	
	//Organize Peaks into Folders
	organizeClusters(sym)
	organizeClusterWaves(sym,os,ovp,ncl)
	
	//Sort the transitions in the cluster waves in ascending order
	sortClusterWaves()
	
	//Identify the DFT transitions with the final transition clusters
	consolidateTransitionClusters(pwave)
End

//This function clusters peaks that are directly adjacent to each other by checking up to two rows
Function simpleCluster2(ow,pw,ot,sym,ncl,Eini,Efin)
	
	Wave ow,pw
	Variable ot,ncl,Eini,Efin
	String sym
	
	Variable i=0,j=0,k=0,n,nTrans = DimSize(pw,0)
	String cpkwName = "clusteredPks_" + sym + num2str(ncl)
	String upkwName = "usedPks_" + num2str(ncl)
	Make/O/T/N=(nTrans) $cpkwName
	Wave/T clusteredPks = $cpkwName
	Make/O/N=(nTrans) $upkwName
	Wave usedPks = $upkwName
	usedPks = 0
	clusteredPks = ""
	String pwName
	
	usedPks = 0
	clusteredPks = ""
	do
		Variable ovp = ow[i][j]
		if(ovp >= ot)
			if(usedPks[j] == 0)
				usedPks[j] = 1
				clusteredPks[k] += num2str(j) + ";"
				j+=1
			else
				j+=1
			endif
		else
			i=j-1
			j=i
			do
				ovp = ow[i][j]
				if(ovp >= ot)
					if(usedPks[j] == 0)
						usedPks[j] = 1
						clusteredPks[k] += num2str(j) + ";"
						j+=1
					else
						j+=1
					endif
				else
					i=j
					k+=1
					break
				endif
				
				if(i >= nTrans || j >= nTrans)
					break
				endif
			while(i < nTrans || j < nTrans)
		endif
		
		if(i >= nTrans || j >= nTrans)
			break
		endif
	while(i < nTrans || j < nTrans)
	
	//Remove empty points
	for(i=nTrans-1;i>=0;i-=1)
		String val = clusteredPks[i]
		if(StringMatch(val,""))
			DeletePoints i,1,clusteredPks
		endif
	endfor
	
	//Combine 1d parameter waves into a single 2d wave
	makePWCluster(clusteredPks,pw,sym)
	makeMergedPW(pw,sym,ncl,Eini,Efin)
	removeOldPW(sym,ncl)
End

//This function clusters peaks that may be separated by at most one matrix element by checking up to two rows
Function simpleCluster3(ow,pw,ot,sym,ncl,Eini,Efin)
	
	Wave ow,pw
	Variable ot,ncl,Eini,Efin
	String sym
	
	Variable i=0,j=0,k=0,n,nTrans = DimSize(pw,0)
	
	String cpkwName = "clusteredPks_" + sym + num2str(ncl)
	String upkwName = "usedPks_" + num2str(ncl)
	//****6/2/2021****
	//1. Populate this wave with the oscillator strengths correponding to the jth transition in a cluster. 
	//Step 1 should be done until we are ready to jump into a new row.
	//2. Find the position of the transition with the highest oscillator strength in this wave(V_maxloc)
	//3. This position should correspond to a list element of the kth point in the wave clusteredPks
	//4. The value of this list element should be the new "i" 
	//****************
	Make/O/T/N=(nTrans) $cpkwName
	Wave/T clusteredPks = $cpkwName
	Make/O/N=(nTrans) $upkwName
	Wave usedPks = $upkwName
	
	usedPks = 0
	clusteredPks = ""
	String pwName
	
	do
		Variable ovp = ow[i][j]
		if(ovp >= ot)	
			if(usedPks[j] == 0) //If peak is not clustered yet
				usedPks[j] = 1	//Mark as clustered
				clusteredPks[k] += num2str(j) + ";"//Record peak in current cluster list
				j+=1	
			else
				j+=1
			endif
		else	//Check the next peak over. Degree of adjacency is 2
			j+=1
			if(i >= nTrans || j >= nTrans)
				//Check that all peaks have been used. If any is missing then set iterator to that value.
				i = checkAllPeaksAreInUse(usedPks,nTrans)
				j = i
				k+=1
				
				if(i >= nTrans || j >= nTrans)
					break							
				endif
			endif
			ovp = ow[i][j]
			if(ovp >= ot)
				if(usedPks[j] == 0)
					usedPks[j] = 1
					clusteredPks[k] += num2str(j) + ";"
					j+=1
				else
					j+=1
				endif
			else	//Switch to a different row. Only check first degree of adjacency here.
				//i=findLastPkUsed(usedPks,nTrans)	//Set this number to equal the last pk that had an overlap greater than ot
				//j=i	
				//**6/3/2021
				i=findMaxNumInList(pw,clusteredPks[k])//Set the row value to the transition in the peaks clustered so far that has highest OS
				j=findLastPkUsed(usedPks,nTrans)//Set this number to equal the last pk that had an overlap greater than ot
				//**********
				do
					ovp = ow[i][j]
					if(ovp >= ot)
						if(usedPks[j] == 0)
							usedPks[j] = 1
							clusteredPks[k] += num2str(j) + ";"
							j+=1
						else
							j+=1
						endif
					else
						i=findFirstUnusedPk(usedPks,nTrans)//Set this value to be the first unused peak
						j=i
						k+=1
						break
					endif
										
					if(i >= nTrans || j >= nTrans)
						break							
					endif
					
				while(i < nTrans || j < nTrans)
			endif
		endif
		
		if(i >= nTrans || j >= nTrans)
		//Check that all peaks have been used. If any is missing then set iterator to that value.
			i = checkAllPeaksAreInUse(usedPks,nTrans)
			j = i
			k+=1
			if(i >= nTrans || j >= nTrans)
				break							
			endif
		endif
	while(i < nTrans || j < nTrans)
	
	//Remove empty points
	for(i=nTrans-1;i>=0;i-=1)
		String val = clusteredPks[i]
		if(StringMatch(val,""))
			DeletePoints i,1,clusteredPks
		endif
	endfor
	
	//Combine 1d parameter waves into a single 2d wave
	makePWCluster(clusteredPks,pw,sym)
	makeMergedPW(pw,sym,ncl,Eini,Efin)
	removeOldPW(sym,ncl)
End

Function makePWCluster(clusteredPks,pw,sym)
	
	Wave/T clusteredPks
	Wave pw
	String sym
	
	//Determine number of clusters
	Variable nClusters = numpnts(clusteredPks),i,j
	for(i=0;i<nClusters;i+=1)
		String pwName = "clusterPW_" + num2str(i) + "_" + sym	
		Variable pksInCluster = ItemsInList(clusteredPks[i])
		Make/O/N=(pksInCluster,16) $pwName
		Wave w = $pwName
		for(j=0;j<pksInCluster;j+=1)
			Variable cpk = str2num(StringFromList(j,clusteredPks[i]))
			w[j][0]  = pw[cpk][0]	//Peak position
			w[j][1]  = pw[cpk][1]
			w[j][2]  = pw[cpk][2]
			w[j][3]  = pw[cpk][3]
			w[j][4]  = pw[cpk][4]
			w[j][5]  = pw[cpk][5]
			w[j][6]  = pw[cpk][6]
			w[j][7]  = pw[cpk][7]
			w[j][8]  = pw[cpk][8]
			w[j][9]  = pw[cpk][9]
			w[j][10] = pw[cpk][10]
			w[j][11] = pw[cpk][11]
			w[j][12] = pw[cpk][12]
			w[j][13] = pw[cpk][13]
			w[j][14] = pw[cpk][14]
			w[j][15] = pw[cpk][15] 
		endfor
	endfor
End

Function/WAVE makeMergedPW(pw,sym,ncl,Eini,Efin)
	
	Wave pw
	String sym
	Variable ncl,Eini,Efin
		
	String pwList = WaveList("clusterPW*",";","")
	Variable n = ItemsInList(pwList),i,j
	String maxVName = "MaxValues" + sym
	String clusterPDName = "percentDifference_" + num2str(ncl)
	Make/O/N=(n) $maxVName,numpks,maxVAdjusted,$clusterPDName
	Wave mvw = $maxVName,clPD = $clusterPDName
	for(i=0;i<n;i+=1)
		String cpw = StringFromList(i,pwList)
		Wave w = $cpw
		//Identify which peak in the cluster has the strongest OS.
		//This peak will be the one used to populate the rest of the information in the final parameter wave
		Variable nPks = DimSize(w,0)
		WaveStats/RMD=[][1,1]/Q w
		numpks[i] = nPks
		if(i==0)
			mvw[i] = V_maxloc - nPks
			maxVAdjusted[i] = V_maxloc - nPks
		else
			mvw[i] = V_maxloc - nPks + sum(numpks,0,i-1)
			maxVAdjusted[i] = V_maxloc - nPks
		endif
		
		Variable amp,wid,pos
		Variable xx=0,yy=0,zz=0,xy=0,xz=0,yz=0,mux=0,muy=0,muz=0,theta=0

		//Make summed peak
		Make/O/N=2000 SPK=0,CPK=0
		SetScale/I x,Eini,Efin,SPK,CPK
		for(j=0;j<nPks;j+=1)
			amp = w[j][1]
			wid = w[j][2]
			pos = w[j][0]
			SPK += amp * gauss(x,pos,wid)//Fix amplitude/normalization problem
		endfor
		//Determine summed peak amplitude, position and width
		WaveStats/Q SPK
		Variable maxAmp = V_max, maxAmploc = V_maxloc,halfMaxAmp = maxAmp/2
		Findlevel/Q/R=(Eini,maxAmploc) SPK,halfMaxAmp
		Variable en1 = V_levelX//,dif = abs(maxAmploc-en1),en2 = maxAmploc+dif
		Findlevel/Q/R=(maxAmploc,Efin) SPK,halfMaxAmp
		Variable en2 = V_levelX,posAdj=(en1+en2)/2
		Variable sigma = (en2-en1)/fwhmconversion//fwhm divided by conversion factor
		
		for(j=0;j<nPks;j+=1)//Sum the TDM components, and tensor elements present in this cluster
			mux += w[j][3]
			muy += w[j][4]
			muz += w[j][5]
			xx  += w[j][8]
			yy  += w[j][9]
			zz  += w[j][10]
			xy  += w[j][11]
			xz  += w[j][12]
			yz  += w[j][13]
		endfor
		
		String pwName = "mergedPW_" + num2str(i)
		Make/O/N=12 $pwName	//This wave contains the position,width,amplitude,xx,yy,zz,xy,xz,yz values of the clustered peaks
		Wave mpw = $pwName
		mpw = 0
		mpw[0] = maxAmploc//posNum/posDen					//Position
		mpw[1] = sqrt(2*pi) * sigma * maxAmp//newAmp	//Amplitude
		mpw[2] = sigma//(widNum/widDen) - mpw[0]	//Width
		mpw[3] = mux		//mux
		mpw[4] = muy		//muy
		mpw[5] = muz		//muz
		mpw[6] = xx		//xx
		mpw[7] = yy		//yy
		mpw[8] = zz		//zz
		mpw[9] = xy		//xy
		mpw[10] = xz		//xz
		mpw[11] = yz		//yz
	//	CPK = mpw[1] * gauss(x,mpw[0],mpw[2])
	//	Variable  pd = calcPercentDiff(SPK,CPK)
	//	if(pd>=0)
	//		Variable amp1,wid1,pos1,posNum=0,widNum=0,Den=0,newAmp=0
	//		for(j=0;j<nPks;j+=1)
	//			amp1 = w[j][1]
	//			wid1 = w[j][2]
	//			pos1 = w[j][0]
		
	//			posNum += amp1*pos1
	//			Den += amp1
	//			widNum += amp1*(wid1 + pos1)
	//		endfor
	//		mpw[0] = posNum/Den					//Position
	//		mpw[2] = (widNum/Den) - mpw[0]  //Width
	//		mpw[1] = sqrt(2*pi) * mpw[2] * maxAmp//ampSum//newAmp				//Amplitude
	//		CPK = mpw[1] * gauss(x,mpw[0],mpw[2])
	//		pd = calcPercentDiff(SPK,CPK)
	//	endif
	//	clPD[i] = pd
	endfor
	
	String pwMergedList = WaveList("mergedPW_*",";","")
	String ccpwName = "combClusterPW" + sym + "_" + num2str(ncl)
	Make/O/N=(n,16) $ccpwName
	Wave combClusterPW = $ccpwName
	for(i=0;i<n;i+=1)
		String cmpw = StringFromList(i,pwMergedList)
		String ccpw = StringFromList(i,pwList)
		Wave w = $cmpw,y = $ccpw
		combClusterPW[i][0]  = w[0]						 //Position		
		combClusterPW[i][1]  = w[1]						 //Amplitude		
		combClusterPW[i][2]  = w[2] 						 //Width
		combClusterPW[i][3]  = w[3]//y[maxVAdjusted[i]][3] //TDMx
		combClusterPW[i][4]  = w[4]//y[maxVAdjusted[i]][4] //TDMy
		combClusterPW[i][5]  = w[5]//y[maxVAdjusted[i]][5] //TDMz
		combClusterPW[i][6]  = y[maxVAdjusted[i]][6] //Theta
		combClusterPW[i][7]  = y[maxVAdjusted[i]][7] //Symmetry		
		combClusterPW[i][8]  = w[6] //xx
		combClusterPW[i][9]  = w[7] //yy
		combClusterPW[i][10] = w[8]//zz
		combClusterPW[i][11] = w[9]//xy
		combClusterPW[i][12] = w[10]//xz
		combClusterPW[i][13] = w[11]//yz
		combClusterPW[i][14] = y[maxVAdjusted[i]][14]//Atom ID
		combClusterPW[i][15] = y[maxVAdjusted[i]][15]//MO
	endfor	
	Wave combClusterPW2 = calcTDMTheta(combClusterPW)//Updates the TDM polar angle for the current set of clusters
	return combClusterPW2
End

Function/WAVE calcTDMTheta(pw)
	Wave pw
	
	Variable mag,th, n=DimSize(pw,0),i
	Make/O/N=3 zVec = {0,0,1}
	for(i=0;i<n;i+=1)
		mag = sqrt(pw[i][8]^2 + pw[i][9]^2 +pw[i][10]^2)
		th = acos((zVec[0]*pw[i][8]+zVec[1]*pw[i][9]+zVec[2]*pw[i][10])/mag)*(180/pi)
		pw[i][6] = th
	endfor
	KillWaves zVec
	return pw
End

//This function returns the last peak that was used after clustering.
Function findLastPkUsed(usedPks,nTrans)

	Wave usedPks
	Variable nTrans
	
	Variable i
	
	for(i=0;i<nTrans;i+=1)
		Variable pk = usedPks[i]
		if(pk == 0)
			break
		endif
	endfor
	
	return i-1
End

//This function returns the first peak that was unused after clustering.It tells the clustering function what row/column to start the next cluster on.
Function findFirstUnusedPk(usedPks,nTrans)

	Wave usedPks
	Variable nTrans
	
	Variable i
	
	for(i=0;i<nTrans;i+=1)
		Variable pk = usedPks[i]
		if(pk == 0)
		//	print i
			break
		endif
	endfor
	
	return i
End

//This function checks whether all the peaks have been used. It tells the clustering function to stop clustering if at the end of the loop all peaks are accounted for
Function checkAllPeaksAreInUse(usedPks,nTrans)
	
	Wave usedPks
	Variable nTrans
	
	Variable i
	
	for(i=0;i<nTrans;i+=1)
		Variable pk = usedPks[i]
		if(pk == 0)
		//	print i
			break
		endif
	endfor
	
	return i
End

Function removeEmptyPointsFromPW()
	
	String pwList = WaveList("clusterPW*",";","")
	Variable n = ItemsInList(pwList),i,j
	
	for(i=0;i<n;i+=1)
		String cpw = StringFromList(i,pwList)
		Wave w = $cpw
		WaveStats/Q w
		Variable np = DimSize(w,0)
		for(j=np-1;j>=0;j-=1)
			Variable val = w[j][0]
			if(V_max == 0)
				KillWaves/Z w
				break
			elseif(val == 0)
				DeletePoints j,1,w
			endif
		endfor
	endfor
End

Function removeOldPW(sym,ncl)
	
	String sym
	Variable ncl
	
	String iniFolder = GetDataFolder(1)
	String cpwfName = iniFolder + "CPW_" + sym + num2str(ncl)
	String mpwfName = iniFolder + "MPW_" + sym + num2str(ncl) 
	NewDataFolder/O $cpwfName
	NewDataFolder/O $mpwfName
	String cpwList = WaveList("clusterPW*",";","")
	String mpwList = WaveList("mergedPW_*",";","")
	Variable ncp = ItemsInList(cpwList),nmp = ItemsInList(mpwList),i
	for(i=0;i<ncp;i+=1)
		String cpw = StringFromList(i,cpwList)
		Wave w = $cpw
		String cpwfName2 = iniFolder + "CPW_" + sym + num2str(ncl) + ":" + cpw
		Duplicate/O w,$cpwfName2
		KillWaves w
	endfor
	
	for(i=0;i<nmp;i+=1)
		String mpw = StringFromList(i,mpwList)
		Wave w = $mpw
		String mpwfName2 = iniFolder + "MPW_" + sym + num2str(ncl) + ":" + mpw
		Duplicate/O w,$mpwfName2
		KillWaves w
	endfor
End

//This function checks the overlap matrix for values above the overlap threshold
//If there are any values above the threshold ir returns 0 and the clustering loop continues
//If there are no values above the threshold it returns 1 and the clustering loop ends
Function determineMaxOVP(ow,ot)

	Wave ow	//Overlap Matrix
	Variable ot	//Overlap Threshold
	
	Duplicate/O ow, tempOW,maskOW
	maskOW = 0

	Variable n = DimSize(ow,0),i=0,j=0
	for(i=j;i<n;i+=1)
		for(j=i;j<n;j+=1)
			if(i==j)
				maskOW[i][j] = 100
			endif
		endfor
	endfor
	
	MatrixOP/O tempOW = ow - maskOW
	WaveStats/Q tempOW
	Variable done
	if(V_max >= ot)
		done = 0
		print "OVP threshold exceeded!"
	else
		done = 1
		print "OVP threshold not exceeded."	
	endif
	KillWaves maskOW,tempOW
	
	return done
End

Function MDsort2(w,keycol)	//Sorts a 2D wave by keycol
    Wave w
    variable keycol
 
    variable ii
  
    make/o/n=(dimsize(w,0)) key
    make/o/n=(dimsize(w,0)) valindex
 
    key[] = w[p][keycol]
    valindex=p
 
    sort key,key,valindex
 
    duplicate/o w, M_newtoInsert
 
    for(ii=0;ii<dimsize(w,0);ii+=1)
        M_newtoInsert[ii][] = w[valindex[ii]][q]
    endfor
 
    duplicate/o M_newtoInsert,w
    killwaves/z key,valindex,M_newtoInsert
End

Function/WAVE makeTotalWave(Eini,Efin,totalPWave,sym)

	Variable Eini,Efin
	Wave totalPwave
	String sym
	
	Variable i

	String gTotName = "gTotal_" + sym 
	Make/O/N=(2000) $gTotName 
	Wave gTotal = $gTotName 
	gTotal = 0
	SetScale/i x,Eini,Efin,gTotal

	Variable npks = DimSize(totalPWave,0),pos,wid,amp
	for(i=0;i<=npks-1;i+=1)
		pos = totalPWave[i][0]
		wid = totalPWave[i][2]
		amp = totalPWave[i][1]
		gTotal += sqrt(2*Pi)*amp*wid*gauss(x,pos,wid)
	endfor
	
	return gTotal
End

Function makePeaks(pWave2D,Eini,Efin,sym)

	Wave pWave2D
	Variable Eini,Efin
	String sym
	
	Variable totPeaks = DimSize(pWave2D,0),i
	Variable pos, wid,amp	
	for(i=0;i<=totPeaks-1;i+=1)
		String pkName = "pk" + num2str(i) + "_Trial" + num2str(i) + "_" + sym
		Make/O/N=2000 $pkName
		Wave w = $pkName
		pos = abs(pWave2D[i][0])
		wid = abs(pWave2D[i][2])
		amp = abs(pWave2D[i][1])
		SetScale/i x,Eini,Efin,w
		w = sqrt(2*Pi)*wid*amp * Gauss(x,pos,wid)
	endfor
End

Function organizeClusters(sym)

	String sym
	
	String fName = "ClusterPeaks_" + sym
	String iniFolder = GetDataFolder(1)
	
	String pkList = WaveList("pk*"+sym,";","")
	Variable npks = ItemsInList(pkList),i
	
	NewDataFolder/O/S $fName
	String pkFolder = GetDataFolder(1)
	
	SetDataFolder $iniFolder
	for(i=0;i<nPks;i+=1)
		String cPk = StringFromList(i,pkList)
		Wave w = $cPk
		String nwName = pkFolder + cPk
		Duplicate/O w, $nwName 
		KillWaves w
	endfor
End

Function/WAVE getOriginalDFTNEXAFS(os,ovp,sym,mol)
	
	Variable os,ovp
	String sym,mol
	
	String iniFolder = GetDataFolder(1)
	String rootFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":"
	String algoFolder = rootFolder + "TransitionFiltering_" + replaceString(".",num2str(os),"p") + "OS_" + replaceString(".",num2str(ovp),"p") + "OVP:"
	String filterFolder = algoFolder + "originalPeaks:"
	String specFolder = filterFolder + "Spectra:tdmCompSpecs:"
	SetDataFolder $specFolder
	Wave Total_Specf0
	SetDataFolder $iniFolder
	
	return Total_Specf0
End

Function/WAVE getFilteredDFTNEXAFS(os,ovp,sym,mol)
	
	Variable os,ovp
	String sym,mol
	
	String iniFolder = GetDataFolder(1)
	String rootFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":"
	String algoFolder = rootFolder + "TransitionFiltering_" + replaceString(".",num2str(os),"p") + "OS_" + replaceString(".",num2str(ovp),"p") + "OVP:"
	String filterFolder = algoFolder + "firstFilter:"
	String specFolder = filterFolder + "Spectra:tdmCompSpecs:"
	SetDataFolder $specFolder
	Wave Total_Specf1
	SetDataFolder $iniFolder
	
	return Total_Specf1
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

Function organizeClusterWaves(sym,tval,ovpMax,ncl)

	String sym
	Variable tval,ovpMax,ncl
	String iniFolder = GetDataFolder(1)
	String fName = "cluster_" + replacestring(".",num2str(tval),"p") + "_OS_" + replacestring(".",num2str(ovpMax),"p") + "_OVP_" + sym
	NewDataFolder/O $fName
	
	String pwName = "combClusterPW" + sym + "_" + num2str(ncl)
	String gtName = "gTotal_" + sym
	String mvName = "MaxValues" + sym
	String owName = "ovpWave" + sym
	
	Wave a = $pwName
	Wave d = $gtName
	Wave f = $mvName
	Wave g = $owName
	
	String pwName2 = iniFolder + fName + ":" + "combClusterPW" + sym// + "_" + num2str(ncl)
	String gtName2 = iniFolder + fName + ":" +  "gTotal_" + sym// + num2str(ncl)
	String mvName2 = iniFolder + fName + ":" +  "MaxValues" + sym// + num2str(ncl)
	String owName2 = iniFolder + fName + ":" +  "ovpWave" + sym// + num2str(ncl)
	
	//Move cluster waves into cluster folder
	Duplicate/O a,$pwName2
	Duplicate/O d,$gtName2
	Duplicate/O f,$mvName2
	Duplicate/O g,$owName2
	
	//Move folder containing cluster peaks into merged folder
	String pkFolderName = "ClusterPeaks_" + sym
	MoveDataFolder/O=1 $pkFolderName, $fName
	
	//Remove unnecessary parameter waves
	String tpwName = "TotalPWave_" + sym
	Wave b = $tpwName
	
	//Remove cluster waves from original folder
	KillWaves a,b,d,f,g
End

Function sortClusterWaves()
	String cPkList = WaveList("clusteredPks_ALL*",";","")//List of clustering iteration waves
	Variable nw = ItemsInList(cPkList),i,j
	for(i=0;i<nw;i+=1)
		String cwName = StringFromList(i,cPkList)
		Wave/T w = $cwName
		Variable n = numpnts(w)
		for(j=0;j<n;j+=1)
			w[j] = SortList(w[j],";",2)
		endfor
	endfor
End

Function consolidateTransitionClusters(pw)

	Wave pw	//Parameter wave after filtering
	
	String cPkList = WaveList("clusteredPks_ALL*",";","")//List of clustering iteration waves
	Variable nw = ItemsInList(cPkList)
	Variable nTransDFT = DimSize(pw,0)//How many DFT transitions are there
	String iniClsName = StringFromList(0,cPkList)//Initial clusters
	Wave/T wI = $iniClsName
	Variable nIniCls = numpnts(wI)
	String finClsName = StringFromList(nw-1,cPkList)//Final clusters
	Wave/T wF = $finClsName
	Variable nFinCls = numpnts(wF),i,j,k //nFinCls defines final number of clusters
	
	for(i=0;i<nTransDFT;i+=1)	//Iterate through DFT Transitions
		String currTrans =  num2str(i)
		for(j=0;j<nIniCls;j+=1)//Iterate through initial clusters
			if(StringMatch(wI[j],"*"+currTrans+"*"))
				pw[i][16] = j
				break
			endif
		endfor
		
		for(j=1;j<nw;j+=1)	//This iterates through subsequent cluster text waves
			String currClusterWave = StringFromList(j,cPkList)
			Wave/T cw = $currClusterWave
			Variable nw2 = numpnts(cw)
			for(k=0;k<nw2;k+=1)//Iterate through current cluster
				if(StringMatch(cw[k],"*"+num2str(pw[i][16])+"*"))
					pw[i][16] = k	//Update cluster number based on previous loop results
					break
				endif
			endfor
		endfor
	endfor
End

Function findMaxNumInList(w,list)
	Wave w	//This should be the 2D parameter wave containing the transition/cluster info
	String list//This should be the list of transitions that make up the cluster
	
	Variable n = ItemsInList(list),i
	Make/O/N=(n) test = 0
	for(i=0;i<n;i+=1)
		//1. Iterate through every transition in the cluster
		//2. Evaluate the 1st column of the 2D parameter wave [OS] at each transition in cluster
		//3. Place OS value from each transition in cluster in wave "test"
		test[i] = w[str2num(StringFromList(i,list))][1]
	endfor
	WaveStats/Q test
	//4. Determine the point in wave "test" that has the maximum value
	//5. This should be the list element that has the maximum OS. This is the output.
	Variable maxpnt = str2num(StringFromList(V_maxloc,list))
	return maxpnt
End