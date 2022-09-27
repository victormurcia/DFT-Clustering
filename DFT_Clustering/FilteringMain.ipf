#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include "SymmetryOperations"
#include "FilteringWrapper"
#include "peakOverlap"
#include "ClusteringProcs"
#include "fittingAmplitudes"
#include "filteringWrapper"
#include "sequentialThresholding"
#include "DFT Plotting Utilities"
#include "makeBareAtomAbsorption"
#include "alphaModeling"

//Developed by:
//Victor Manuel Murcia Ruiz, PhD
//Washington State University 
//2018
//**********************************************************************************************************************
///Routines used to carry out the filtering and subsequent clustering of NEXAFS transitions calculated via Density Functional Theory calculations
//using the StoBe software package as well as fitting these transitions to experimental NEXAFS.
//Filtering is done by defining the following conditions:
//
//1. Energy range -- Within what energy range do we want to consider DFT transitions?
//2. Oscillator Strength -- What's the minimum intensity a transition must have to be considered relevant?
//3. Peak Overlap -- What's the maximum overlap that a pair of peaks can have before being considered as indistinguishable?
//4. Transition Symmetry -- What kind of symmetry defines a transition as determined from the absorption tensor for each resonance?
//
//The transitions that remain after the filtering are then clustered together into a single representative peak via an algorithm that pools the peak position, width and amplitude
//of the transitions that make that cluster.
//
//The amplitudes of the clustered transitions are then fitted iteratively so that they match the intensity of the original (i.e. unfiltered) DFT spectrum.
//
//The NEXAFS generated at this point is for the molecular orientation that was defined by the molecular configuration used to run the DFT calculations. The next step of this algorithm
//involves generating the absorption tensor for ALL molecular orientations which allows for the generation of the FILM absorption tensor. This is done by applying a tilt angle (alpha) 
//with respect to the z-axis to the absorption tensor. This tilted tensor is then rotated azimuthally AROUND the z-axis. Each of the azimuthally rotated tensors is added to generate the
//TOTAL film absorption tensor.   
//
//The film tensor can then be impinged upon with an electric field (E) coming from an orientation specified by a polar angle (theta) and azimuthal angle (phi) which allows for the simulation
//of angle resolved NEXAFS. The equation for the electric field is  E = ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) )
//
//The amplitude fitted DFT spectrum then has a step edge superimposed that is made from the bare atom absorption profile made for the molecule being studied 

Function filterDFT(atomName,fnum,tval,ovpmax,maxRot,LUMOwave,[broadShift,erange,broad1,broad2,Eini,Efin,d,gRes,buildSym,corrWave,dClusters,realign,thx,thy,thz,rotOrder,modelAlpha,IPwave,expSpecName,expEnergyName,expFolderPath,mol,alpha,i0,phi,fit,rigidShift,stepShift,thetaList,NEXAFStype,startPre,startPost,endPre,endPost,justModel,justFit,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,holdPos,holdWidths,holdAmps,refinement,maskEnergy1,maskEnergy2,pkToRefine,refitPWave,holdTensorElems])
	
	String atomName //What kind of atom are we processing?
	Variable fnum //How many atoms to process?
	Variable tval //Set a value between 0 and 100%. 
	Variable ovpMax //Peak overlap threshold, number between 0 and 100
	Variable maxRot	//What's the symmetry of the molecule? Currently supported (1 = no symmetry, 2 = 2-fold rotation [180 degrees], 3 = 3 fold rotation [120 degrees], 4 = 4-fold rotation [90 degrees])
	Wave LUMOwave
	//Optional Parameters for Clustering/Filtering Algorithm
	Variable erange //Maximum energy range to consider transitions in
	Variable broad1 //Peak width at lower energy (ewid1)
	Variable broad2 //Peak width at higher energy (ewid2)
	Variable Eini	//Starting energy at which to make transition peaks
	Variable Efin	//Ending energy at which to make transition peaks
	Variable broadShift //Allows for the adjustment of the energies at which to apply the broadening scheme
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
	
	print "---------------------------------------------------------------------------------------------------------------------------------------------"
	print "--------------------------------Optical Tensor Modeling of NEXAFS via Filtering and Clustering of DFT calculations---------------------------"
	print "Victor Manuel Murcia Ruiz ---------Department of Materials Science and Engineering, Washington State University---------victor.murcia@wsu.edu"
	print "Brian A. Collins -----------------------Department of Physics and Astronomy, Washington State University----------------brian.collins@wsu.edu"
	print "---------------------------------------------------------------------------------------------------------------------------------------------" 
	
	print "Initializing algorithm..."
	
	if(!ParamIsDefault(modelAlpha))
		if(StringMatch(modelAlpha,"yes"))
			print "Molecular tilt angle will be modeled."
			print "Checking that necessary waves are present..."
			if(ParamIsDefault(IPwave))
				print "Please provide wave containing Ionization Potentials before proceeding."
				Abort
			endif
		
			if(ParamIsDefault(expSpecName))
				print "Please provide the base name for the experimental angle resolved NEXAFS waves before proceeding."
				Abort
			endif
		
			if(ParamIsDefault(expEnergyName))
				print "Please provide the base name for the experimental energy waves before proceeding."
				Abort
			endif
		
			if(ParamIsDefault(expFolderPath))
				print "Please provide the IGOR path to the experimental NEXAFS waves and experimental energy waves before proceeding."
				Abort
			endif
		
			if(ParamIsDefault(mol))
				print "Please specify the name of the molecule to be processed before proceeding. The name of the molecule should be the same name used when building the bare atom absorption."
				Abort
			endif
		
			String baWave = "root:Packages:NXA:" + mol + "_mu"
			if(!WaveExists($baWave))
				print "The bare atom absorption for " + mol + " was not in the folder root:Packages:NXA: . Please make sure that the bare atom absorption wave was constructed prior to starting algorithm."
			endif
		
			if(ParamIsDefault(thetaList))
				print "Please enter a list of semi-colon separated values corresponding to the NEXAFS alignments for the experimental data before proceeding."
				Abort
			endif
		
			print "All necessary waves are present. Starting filtering algorithm"
		else
			print "Molecular tilt angle will not be modeled."
			print "Starting filtering algorithm."
		endif
	endif
	
	
	//Make Global Variables to store results to be used in paramspace exploration
	Variable/G root:perDiff
	Variable/G root:numpeaks
	Variable/G root:redchiSqBBCl
	Variable/G root:redchiSqBBCluref
	Variable/G root:redchiSqExpCl
	Variable/G root:chiSqBBCl
	Variable/G root:chiSqBBCluref
	Variable/G root:chiSqExpCl
	Variable/G root:GFBBCl
	Variable/G root:GFBBCluref
	Variable/G root:GFExpCl
	Variable/G root:compTime
	
	if(ParamIsDefault(modelAlpha))
		modelAlpha = "no"
	endif
	
	//holdSecAmp=holdTensorElems
	String startFolder = "root:Packages:DFTClustering:PolarAngles_" + mol +":"
	String baseFolderName = startFolder + "TransitionFiltering_" + replacestring(".", num2str(tval),"p") + "OS_" + replacestring(".", num2str(ovpMax),"p") + "OVP"
	String fitFolder = baseFolderName + ":AmplitudeFitting:"
	print fitFolder
	//Just model the NEXAFS with the BB model using the results of a previous clustering run
	if(justModel)
		SetDataFolder $fitFolder
		print "Modeling the NEXAFS using results obtained by clustering with a " + num2str(tval) + "% OS threshold and a " + num2str(ovpMax) + "% OVP threshold."
		simDFT(tval,ovpmax,IPwave,mol,alpha,i0,phi,expSpecName,expEnergyName,expFolderPath,"no",rigidShift,thetaList,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,holdAmps=holdAmps,holdWidths=holdWidths,holdPos=holdPos,d=1,NEXAFStype=NEXAFStype,holdModTheta=holdTensorElems)
	//Just fit the NEXAFS with the BB model using the results of a previous clustering run
	elseif(justFit)
		print "Fitting the NEXAFS using results obtained by clustering with a " + num2str(tval) + "% OS threshold and a " + num2str(ovpMax) + "% OVP threshold."
		//NEED TO FIX THIS SO THAT IT'S JUST FOR SINGLE PEAK
		if(refinement)
			String refitFolder = fitFolder + "Alpha_" + num2str(alpha) + ":TensorMisc:"
			SetDataFolder $refitFolder
			simDFT(tval,ovpmax,IPwave,mol,alpha,i0,phi,expSpecName,expEnergyName,expFolderPath,"yes",rigidShift,thetaList,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,holdAmps=holdAmps,holdWidths=holdWidths,holdPos=holdPos,d=1,NEXAFStype=NEXAFStype,refine=refinement,maskEnergy1=maskEnergy1,maskEnergy2=maskEnergy2,pkToRefine=pkToRefine,refpw=refitPWave,holdModTheta=holdTensorElems)
		else
			SetDataFolder $fitFolder
			String dfName = "Alpha_" + num2str(alpha)
			if(DataFolderExists(dfName))
				String summaryPlot = "Final_Comparison" + replaceString(".",num2str(alpha),"p") + "_OS" + replaceString(".",num2str(tval),"p") + "_OVP" + replaceString(".",num2str(ovpMax),"p") + "_" + NEXAFStype
				String gName2 = "FitAmplitudes_Alpha_" + replaceString(".",num2str(alpha),"p") + "OS_"  + replaceString(".",num2str(tval),"p") + "OVP_" + replaceString(".",num2str(ovpMax),"p")
				String corrGName = "CORRELATION_" + replaceString(".",num2str(alpha),"p")
				DoWindow/K $summaryPlot
				DoWindow/K $gName2
				DoWindow/K $corrGName
				KillDataFolder $dfName
			endif	
			
			simDFT(tval,ovpmax,IPwave,mol,alpha,i0,phi,expSpecName,expEnergyName,expFolderPath,"yes",rigidShift,thetaList,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,holdAmps=holdAmps,holdWidths=holdWidths,holdPos=holdPos,d=1,NEXAFStype=NEXAFStype,maskEnergy1=maskEnergy1,maskEnergy2=maskEnergy2,pkToRefine=pkToRefine,holdModTheta=holdTensorElems)
		endif 
	else//Run the clustering algorithm from the start and potentially fit/model the NEXAFS using the DFT generated BB model.
		NewDataFolder/O $baseFolderName//Make data folder to contain results from algorithm
	
		//Get wave with relevant symmetry info ready (TDM components, theta, symmetry group and OS derived from each TDM component)
		if(StringMatch(buildSym,"yes"))
			print "Determining symmetry elements"
			if(!realign)	
				prepThetaSym2(fnum,startFolder,atomName,maxRot)
			else
				prepThetaSym2(fnum,startFolder,atomName,maxRot,realign=1,thx=thx,thy=thy,thz=thz,rotOrder=rotOrder)
			endif
		endif
	
		//Extract the energies, oscillator strengths and TDM thetas from DFT importing stage
		extractDFT(fnum,startFolder,baseFolderName,atomName)
	
		//This wrapper function takes care of all the operations involving the filtering of transitions via the thresholding of the Oscillator Strength, the transition overlap, and the transition symmetry
		filteringRoutines(fnum,startFolder,baseFolderName,atomName,erange,tval,broad1,broad2,Eini,Efin,gRes,ovpMax,corrWave,dClusters,LUMOwave,IPwave,mol)
		
		Wave allParamsF1Sorted
		WaveStats/Q/RMD=[][7,7] allParamsF1Sorted
		Variable symmetry = V_max
		if(symmetry == 0)
			print "System has isotropic symmetry"
		elseif(symmetry == 1)
			print "System has uniaxial symmetry"
		elseif(symmetry == 2)
			print "System has biaxial symmetry"
		elseif(symmetry == 3)
			print "System has triaxial symmetry"
		else
			print "Symmetry couldn't be determined. Check parameter wave."	
		endif
		
		SetDataFolder $baseFolderName
		
		if(d)
			print "Making plots..."
			dftPlots(baseFolderName,tval,ovpMax)
		endif
	
		//Wrapper function that fits clustered BB peaks to original DFT NEXAFS
		Variable pDiff1 = ampFittingWrapper(baseFolderName,tval,ovpMax,Eini,erange,gRes)
		//Variable pDiff1 = ampFittingWrapper(baseFolderName,tval,ovpMax,Eini,Efin,gRes)
		
		//Wrapper function that fits DFT BB Model to Experimental NEXAFS
		if(StringMatch(modelAlpha,"yes"))				
			//****Change this into a structure --> More user accessible****
			simDFT(tval,ovpmax,IPwave,mol,alpha,i0,phi,expSpecName,expEnergyName,expFolderPath,fit,rigidShift,thetaList,anchorStep1,anchorStep2,anchorExp1,anchorExp2,stepWid1,stepWid2,stepE1,stepE2,holdAmps=holdAmps,holdWidths=holdWidths,holdPos=holdPos,d=1,NEXAFStype=NEXAFStype,maskEnergy1=maskEnergy1,maskEnergy2=maskEnergy2,pkToRefine=pkToRefine,holdModTheta=holdTensorElems)
		endif
	endif
	if(justFit)
		plotRawDFTvsClusterDFT(tval,ovpMax,mol)
		plotParamChanges(tval,ovpMax,alpha,mol)
		plotDFTBBvsEXP(thetaList,tval,ovpMax,alpha,mol)
	elseif(StringMatch(modelAlpha,"yes"))	
		plotRawDFTvsClusterDFT(tval,ovpMax,mol)
		plotParamChanges(tval,ovpMax,alpha,mol)
		plotDFTBBvsEXP(thetaList,tval,ovpMax,alpha,mol)	
	else	
		plotRawDFTvsClusterDFT(tval,ovpMax,mol)
	endif
	
	return pDiff1
	
	SetDataFolder root:
End

Function extractDFT(fnum,iniFolder,baseFolderName,atomName)
	
	Variable fnum
	String iniFolder,baseFolderName,atomName

	Variable i
	for(i=1;i<=fnum;i+=1)
		String cAtom = atomName + num2str(i)
		String cFolder = iniFolder + cAtom
		SetDataFolder $cFolder
		String cEnName = "eV_" + cAtom
		String osWaveName = "OS_" //+ cAtom
		String symmWaveName = "symmetryParams" + cAtom
		Wave eV_, symmetryParams,currOS = $osWaveName
		String enOrig = baseFolderName  + ":" + "eV_" + cAtom
		String osOrig = baseFolderName  + ":" + "OS_" + cAtom
		String symOrig = baseFolderName + ":" + "symmetryParams" +cAtom
		if(WaveExists($osWaveName) ==1)
			Duplicate/O currOS, $osOrig
			Duplicate/O eV_, $enOrig
			Duplicate/O symmetryParams, $symOrig
		elseif(WaveExists($osWaveName) == 0)
			print "Couldn't find Oscillator Strength wave"
		elseif(WaveExists(eV_) == 0)
			print "Couldn't find Energy wave"
		elseif(WaveExists(symmetryParams) == 0)
			print "Couldn't find symmetry wave"
		endif
	endfor
	
	SetDataFolder baseFolderName
End

Function ApplyColorTableToTopGraph(ctabname)
    String ctabname

    String graphName = WinName(0, 1)
    if (strlen(graphName) == 0)
        return -1
    endif

    Variable numTraces = ItemsInList(TraceNameList(graphName,";",3))

    if (numTraces <= 0)
        return -1
    endif
   
    Variable denominator= numTraces-1
    if( denominator < 1 )
        denominator= 1    // avoid divide by zero, use just the first color for 1 trace
    endif

    ColorTab2Wave $ctabname // creates M_colors
    Wave M_colors
    Variable numRows= DimSize(M_colors,0)
    Variable red, green, blue
    Variable i, index
    for(i=0; i<numTraces; i+=1)
        index = round(i/denominator * (numRows-1))  // spread entire color range over all traces.
        ModifyGraph/W=$graphName rgb[i]=(M_colors[index][0], M_colors[index][1], M_colors[index][2])
    endfor
    return 0
End
