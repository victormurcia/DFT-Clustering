#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//This function calculates the percent overlap between two unnormalized Gaussian peaks.
//It does this by determining the intersection point between the two Gaussians and then evaluating the
//error function for each peak at the intersection point. 
Function/WAVE peakOverlap(pWave,ovpWaveID,[ncl])
	
	Wave pWave
	String ovpWaveID	//What amplitude to use?		
	Variable ncl //Clustering loop stage
	//print "Starting calculations of overlap between peaks..."
	
	Variable nTrans = DimSize(pWave,0)
	String ovpWaveName
	if(ParamIsDefault(ncl))//If the overlap matrix is being made outside of the clustering stage name the overlap wave as so
		ovpWaveName = "ovpWave" + ovpWaveID
	else//Otherwise, include the clustering iteration that we are on as part of the name of the overlap matrix
		ovpWaveName = "ovpWave" + ovpWaveID + num2str(ncl)
	endif
	Make/O/N=(nTrans,nTrans) $ovpWaveName
	Wave ovpWave    = $ovpWaveName
	
	Variable i,j,ovpArea,mu1,mu2,std1,std2,amp1,amp2,pOVP,wid=0.4
	
	for(i=0;i<=nTrans-1;i+=1)
		mu1  = pWave[i][0]
		if(ParamIsDefault(ncl))//If first calculation of overlap, use a narrower and constant width for the peaks, else use standard width values from broad scheme
			std1 = wid/fwhmConversion
		else
			std1 = pWave[i][2]
		endif
		amp1 = pWave[i][1]
	
		for(j=0;j<=nTrans-1;j+=1)
			mu2 = pWave[j][0]
			if(ParamIsDefault(ncl))
				std2 = wid/fwhmConversion
			else
				std2 = pWave[j][2]
			endif
			amp2 = pWave[j][1]
  		
			if(i==j)//A peak always has 100% overlap with itself
				ovpWave[i][j] = 100
			elseif((mu1 == mu2) && (std1 == std2) && (amp1 == amp2))//If two peaks have the same position, width and amplitude then their overlap is 100%
				ovpWave[i][j] = 100
			elseif(amp1==0 || amp2==0)	//If the amplitude for either transition is 0, then there is no overlap
				ovpWave[i][j] = 0
			elseif(mu1 != mu2)	//If none of the above conditions are met AND the peak positions are different then calculate the overlap using the error function method [numerically]
				pOVP = compGaussErf3(mu1,std1,amp1,mu2,std2,amp2,2000,1)
				if(pOVP < 0.01)
					pOVP = 0
				elseif(pOVP > 100)
					pOVP = 100
				endif
				ovpWave[i][j] = pOVP
			else	//If none of the above conditions are met AND the positions are the same then calculate the overlap analytically
				Variable area1 = amp1 * std1 * sqrt(2*Pi)
				Variable area2 = amp2 * std2 * sqrt(2*Pi)
				Variable totArea = min(area1 , area2)//totArea = (area1 + area2)/2
				pOVP = (area1/totArea) * 100
				if(pOVP < 0.01)
					pOVP = 0
				elseif(pOVP > 100)
					pOVP = 100	
				endif
				ovpWave[i][j] = pOVP
			endif
		endfor
	endfor
	
	return ovpWave
End

//This function calculates the percent overlap between two unnormalized Gaussian peaks.
//It does this by determining the intersection point between the two Gaussians and then evaluating the
//error function for each peak at the intersection point. 
Function compGaussErf3(mu1,sd1,amp1,mu2,sd2,amp2,gRes,normal)

	Variable mu1,sd1,amp1,mu2,sd2,amp2,gRes,normal

	Variable c,a1,a2,totAreaA,areaA1,areaA2,ovpArea
	//Different cases for solving for c		
	if(normal)
		if((mu1 == mu2) && (sd1 == sd2) && (amp1 == amp2))
			ovpArea = 100
		elseif(sd1 == sd2)
			c = (mu1+mu2)/2
			if(mu1 < mu2)
				a1 = (1/2)*(1 + erf((mu1 - c)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((c - mu2)/(sd2*sqrt(2))))
			elseif(mu1 > mu2)
				a1 = (1/2)*(1 + erf((c-mu1)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((mu2-c)/(sd2*sqrt(2))))
			endif
		else
			if(mu1 < mu2)
				c = (mu2*sd1^2 - sd2*(mu1*sd2 + sd1*sqrt((mu1-mu2)^2 + 2*(sd1^2 - sd2^2)*ln(sd1/sd2))))/(sd1^2 - sd2^2)
				a1 = (1/2)*(1 + erf((mu1 - c)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((c - mu2)/(sd2*sqrt(2))))
			elseif(mu1 > mu2)
				c = (mu1*sd2^2 - sd1*(mu2*sd1 + sd2*sqrt((mu2-mu1)^2 + 2*(sd2^2 - sd1^2)*ln(sd2/sd1))))/(sd2^2 - sd1^2)
				a1 = (1/2)*(1 + erf((c-mu1)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((mu2-c)/(sd2*sqrt(2))))
			endif
		endif
	//	Variable a1 = (1/2) + (1/2)*erf((c - mu1)/(sd1*sqrt(2)))	//Area of Gauss1
	//	Variable a2 = (1/2) + (1/2)*erf((c - mu2)/(sd2*sqrt(2)))	//Area of Gauss2
		areaA1 = sd1*sqrt(2*Pi)
		areaA2 = sd2*sqrt(2*Pi)
		totAreaA = ((1/(sd1*sqrt(2*Pi)))*sd1*sqrt(2*Pi) + (1/(sd2*sqrt(2*Pi)))*sd2*sqrt(2*Pi))/2
		ovpArea  = ((a1 + a2)/totAreaA)*100
	elseif(!normal)
	//*/*/*/*
		if(sd1 == sd2)
			c = (-2*ln(amp1/amp2)*sd1^2 + mu1^2 - mu2^2)/(2*mu1 - 2*mu2)
		else
			Variable radical
			//Decide which root to take based on transition energies of peaks
  			if(mu1 < mu2)
  				radical = -sqrt(2*ln(amp1/amp2)*(sd2^2 - sd1^2) + mu1^2 -2*mu1*mu2 + mu2^2)
  		 	elseif(mu1 > mu2)  
  		 		radical =  sqrt(2*ln(amp1/amp2)*(sd2^2 - sd1^2) + mu1^2 -2*mu1*mu2 + mu2^2)
  		 	endif	
  			Variable numerator = (mu2*sd1/sd2) - (mu1*sd2/sd1) + radical
  			Variable denominator = (sd1/sd2)- (sd2/sd1)
 			c = numerator/denominator		
		endif
		
		//Check if "c" is +/-Inf or a NaN
		//If c is +/-Inf or a NaN, then there is no intersection point between Gaussians and the overlap is 100%
		if(numtype(c) != 0)
			a1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((mu1)/(sd1*sqrt(2))))
			a2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((mu2)/(sd2*sqrt(2))))
			
			areaA1 = amp1*sd1*sqrt(2*Pi)
			areaA2 = amp2*sd2*sqrt(2*Pi)
			totAreaA = min(areaA1 , areaA2)
						
			if(amp1 > amp2)
				ovpArea= ((a2)/totAreaA)*100
			elseif(amp1 < amp2)
				ovpArea= ((a1)/totAreaA)*100
			endif
		
		else
			if(mu1 < mu2)
				a1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((mu1 - c)/(sd1*sqrt(2))))
				a2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((c - mu2)/(sd2*sqrt(2))))
			elseif(mu1 > mu2)
				a1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((c-mu1)/(sd1*sqrt(2))))
				a2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((mu2-c)/(sd2*sqrt(2))))
			endif
		
			areaA1 = amp1*sd1*sqrt(2*Pi)
			areaA2 = amp2*sd2*sqrt(2*Pi)
			totAreaA = min(areaA1 , areaA2)
			//totAreaA = (amp1*sd1*sqrt(2*Pi) + amp2*sd2*sqrt(2*Pi))/2
			ovpArea  = ((a1 + a2)/totAreaA)*100
		endif
	endif
	//ovpArea = round(ovpArea)
	return ovpArea
End

//This function also calculates the percent overlap between two Gaussians. However, it also creates
//the error function and gaussian waves for visualization purposes
Function compGaussErf2(mu1,sd1,amp1,mu2,sd2,amp2,gRes,normal,[d,pk1,pk2])

	Variable mu1,sd1,amp1,mu2,sd2,amp2,gRes,normal,d,pk1,pk2
	Make/O/N=(gRes) xwG,gWave1,erfWave1,gWave2,erfWave2,gWave1Int,gWave2Int,gWaveDif
	Variable eIni = 280,eFin = 360
	SetScale/i x,eIni,eFin,gWave1
	SetScale/i x,eIni,eFin,erfWave1
	SetScale/i x,eIni,eFin,gWave2
	SetScale/i x,eIni,eFin,erfWave2
	SetScale/i x,eIni,eFin,gWaveDif
	SetScale/i x,eIni,eFin,gWave1Int
	SetScale/i x,eIni,eFin,gWave2Int
	
	Variable res =(eFin-eIni)/gRes
	Variable i
	for(i=0;i<=gRES-1;i+=1)
		xwG[i] = 280 + i*res
	endfor
	
	//Is it a standard normal distribution or is it an amplitude modified Gaussian?
	if(normal)
		if(mu1 <= mu2)
			gWave1   = gauss(xwG,mu1,sd1)
			erfWave1 = (1/2)+(1/2)*erf((mu1-xwG)/(sd1*sqrt(2)))
			gWave2   = gauss(xwG,mu2,sd2)
			erfWave2 = (1/2)+(1/2)*erf((xwG-mu2)/(sd2*sqrt(2)))
		elseif(mu1 >= mu2)
			gWave1   = gauss(xwG,mu1,sd1)
			erfWave1 = (1/2)+(1/2)*erf((xwG-mu1)/(sd1*sqrt(2)))
			gWave2   = gauss(xwG,mu2,sd2)
			erfWave2 = (1/2)+(1/2)*erf((mu2-xwG)/(sd2*sqrt(2)))
		endif
	elseif(!normal)
		if(mu1 <= mu2)
			//Make Gaussian Peak and error function using mu1,sd1,amp1
			gWave1   = amp1*sd1*sqrt(2*Pi)*gauss(x,mu1,sd1)
			erfWave1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((mu1-x)/(sd1*sqrt(2))))
			//Make Gaussian Peak and error function using mu2,sd2,amp2
			gWave2   = amp2*sd2*sqrt(2*Pi)*gauss(x,mu2,sd2)
			erfWave2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((x-mu2)/(sd2*sqrt(2))))
		elseif(mu1 >= mu2)
			//Make Gaussian Peak and error function using mu1,sd1,amp1
			gWave1   = amp1*sd1*sqrt(2*Pi)*gauss(x,mu1,sd1)
			erfWave1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((x-mu1)/(sd1*sqrt(2))))
			//Make Gaussian Peak and error function using mu2,sd2,amp2
			gWave2   = amp2*sd2*sqrt(2*Pi)*gauss(x,mu2,sd2)
			erfWave2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((mu2-x)/(sd2*sqrt(2))))
		endif
	endif
	
	//Calculate the integral of gWave1. The integral of a gaussian is the error function
	//Therefore, gWaveInt should be identical to errWave
	Integrate gWave1/D=gWave1Int
	Integrate gWave2/D=gWave2Int
	
	if(mu1 <= mu2)
		gWave1Int -= WaveMax(gWave1Int) 
		gWave1Int*=-1
	elseif(mu1 >= mu2)
		gWave2Int -= WaveMax(gWave2Int) 
		gWave2Int*=-1
	endif
	
	gWaveDif = min(gWave1,gWave2)
	Variable difArea = area(gWaveDif)
	Variable g1Area  = area(gWave1)
	Variable g2Area  = area(gWave2)
	Variable totArea = (g1Area+g2Area)/2
	print "The average area between the two Gaussians calculated by integration is:",totArea
	Variable pOVP = (difArea/totArea)*100,a1,a2
	print "The percent overlap calculated by integration is:",pOVP,"%"
	if(normal)
		if(sd1 == sd2)
			Variable c = (mu1+mu2)/2
			if(mu1 < mu2)
				a1 = (1/2)*(1 + erf((mu1 - c)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((c - mu2)/(sd2*sqrt(2))))
			elseif(mu1 > mu2)
				a1 = (1/2)*(1 + erf((c-mu1)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((mu2-c)/(sd2*sqrt(2))))
			endif
		else
			if(mu1 <= mu2)
				c = (mu2*sd1^2 - sd2*(mu1*sd2 + sd1*sqrt((mu1-mu2)^2 + 2*(sd1^2 - sd2^2)*ln(sd1/sd2))))/(sd1^2 - sd2^2)
				a1 = (1/2)*(1 + erf((mu1 - c)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((c - mu2)/(sd2*sqrt(2))))
			elseif(mu1 > mu2)
				c = (mu1*sd2^2 - sd1*(mu2*sd1 + sd2*sqrt((mu2-mu1)^2 + 2*(sd2^2 - sd1^2)*ln(sd2/sd1))))/(sd2^2 - sd1^2)
				a1 = (1/2)*(1 + erf((c-mu1)/(sd1*sqrt(2))))
				a2 = (1/2)*(1 + erf((mu2-c)/(sd2*sqrt(2))))
			endif
		endif
	//	Variable a1 = (1/2) + (1/2)*erf((c - mu1)/(sd1*sqrt(2)))	//Area of Gauss1
	//	Variable a2 = (1/2) + (1/2)*erf((c - mu2)/(sd2*sqrt(2)))	//Area of Gauss2
		Variable areaA1 = sd1*sqrt(2*Pi)
		Variable areaA2 = sd2*sqrt(2*Pi)
		Variable totAreaA = ((1/(sd1*sqrt(2*Pi)))*sd1*sqrt(2*Pi) + (1/(sd2*sqrt(2*Pi)))*sd2*sqrt(2*Pi))/2
		Variable ovpArea  = ((a1 + a2)/totAreaA)*100
		print "The percent overlap calculated via the CDF is:   ",ovpArea//(totAreaA - a1 + a2)*100,"%"
	elseif(!normal)
	//Different cases for solving for c		
		if(sd1 == sd2)
			c = (-2*ln(amp1/amp2)*sd1^2 + mu1^2 - mu2^2)/(2*mu1 - 2*mu2)
		else
			Variable radical
			//Decide which root to take based on transition energies of peaks
  	 		if(mu1 < mu2)
  	 			radical = -sqrt(2*ln(amp1/amp2)*(sd2^2 - sd1^2) + mu1^2 -2*mu1*mu2 + mu2^2)
  		 	elseif(mu1 > mu2)  
  		 		radical =  sqrt(2*ln(amp1/amp2)*(sd2^2 - sd1^2) + mu1^2 -2*mu1*mu2 + mu2^2)
  		 	endif	
  			Variable numerator = (mu2*sd1/sd2) - (mu1*sd2/sd1) + radical
  			Variable denominator = (sd1/sd2)- (sd2/sd1)
 			c = numerator/denominator		
		endif
		
		//Check if "c" is +/-Inf or a NaN
		//If c is +/-Inf or a NaN, then there is no intersection point between Gaussians and the overlap is 100%
		if((mu1 == mu2) && (sd1 == sd2) && (amp1 == amp2))
			ovpArea = 100
			print "Peaks are identical. The percent overlap is: 100%"
		elseif(numtype(c) != 0)
			print "There is no intersection point between the two Gaussians"
			a1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((mu1)/(sd1*sqrt(2))))
			a2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((mu2)/(sd2*sqrt(2))))
			
			areaA1 = amp1*sd1*sqrt(2*Pi)
			areaA2 = amp2*sd2*sqrt(2*Pi)
			totAreaA = min(areaA1 , areaA2)
			//Variable totAreaA2 = (areaA1 + areaA2)/2
			
			if(amp1 > amp2)
				ovpArea= ((a2)/totAreaA)*100
			elseif(amp1 < amp2)
				ovpArea= ((a1)/totAreaA)*100
			endif
			
			print "The percent overlap calculated via the CDF is:   ",ovpArea,"%"
		else
			print "The intersection point between the two Gaussians is: ",c
			if(mu1 < mu2)
				a1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((mu1 - c)/(sd1*sqrt(2))))
				a2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((c - mu2)/(sd2*sqrt(2))))
			elseif(mu1 > mu2)
				a1 = amp1*sd1*sqrt(2*Pi)*(1/2)*(1 + erf((c-mu1)/(sd1*sqrt(2))))
				a2 = amp2*sd2*sqrt(2*Pi)*(1/2)*(1 + erf((mu2-c)/(sd2*sqrt(2))))
			endif
			
			areaA1 = amp1*sd1*sqrt(2*Pi)
			areaA2 = amp2*sd2*sqrt(2*Pi)
			totAreaA = min(areaA1 , areaA2) //
			//totAreaA = (amp1*sd1*sqrt(2*Pi) + amp2*sd2*sqrt(2*Pi))/2
			ovpArea  = ((a1 + a2)/totAreaA)*100
			print "The percent overlap calculated via the CDF is:   ",ovpArea,"%"
		endif
	endif
	//Need to account for cases when difference in peak energy[mu] is very small
	//Perhaps use OS_Threshold? Evaluate pk intensity at c? 
	if(d)
		DoWindow ovpVisualizerPlot
		if(!V_Flag)
			Display/N=ovpVisualizerPlot/K=1 gWaveDif,gWave1,gWave2
		//	ModifyGraph/W=ovpVisualizerPlot margin(right)=108,lstyle(erfWave1)=3,lsize(erfWave1)=2,lstyle(erfWave2)=3,lsize(erfWave2)=2,rgb(erfWave2)=(0,0,0)
			ModifyGraph/W=ovpVisualizerPlot mode(gWaveDif)=7,hbFill(gWaveDif)=4,rgb(gWaveDif)=(1,39321,39321,32768),lsize(gWave1)=2,lsize(gWave2)=2,rgb(gWave2)=(0,0,0)
			ModifyGraph/W=ovpVisualizerPlot mirror=1,minor=1,fStyle=1,fSize=12
			Label/W=ovpVisualizerPlot left "Transition Intensity[a.u.]\\U"
			Label/W=ovpVisualizerPlot bottom "Transition Energy[eV]"
			if(mu1 > mu2)
				SetAxis/W=ovpVisualizerPlot bottom mu2 - sd2*3,mu1+sd1*3
			elseif(mu1 < mu2)
				SetAxis/W=ovpVisualizerPlot bottom mu1 - sd1*3,mu2+sd2*3
			elseif(mu1 == mu2)
				SetAxis/W=ovpVisualizerPlot bottom mu1 - sd1*3,mu1+sd1*3
			endif
			Legend/W=ovpVisualizerPlot/C/N=text0/J/X=0.75/Y=5.76/A=RC/E "\\s(gWave1) Pk"+num2str(pk1)+" \\s(gWave2) Pk"+num2str(pk2)+"\r\\s(gWaveDif) OVP = "+num2str(round(ovpArea))+"%"
		else
			if(mu1 > mu2)
				SetAxis/W=ovpVisualizerPlot bottom mu2 - sd2*3,mu1+sd1*3
			elseif(mu1 < mu2)
				SetAxis/W=ovpVisualizerPlot bottom mu1 - sd1*3,mu2+sd2*3
			elseif(mu1 == mu2)
				SetAxis/W=ovpVisualizerPlot bottom mu1 - sd1*3,mu1+sd1*3
			endif
			Legend/W=ovpVisualizerPlot/C/N=text0/J/X=0.75/Y=5.76/A=RC/E "\\s(gWave1) Pk"+num2str(pk1)+" \\s(gWave2) Pk"+num2str(pk2)+"\r\\s(gWaveDif) OVP = "+num2str(round(ovpArea))+"%"
		endif
	endif
	return ovpArea
End
