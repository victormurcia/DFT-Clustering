#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <WaveSelectorWidget>
#include <PopupWaveSelector>
#include "Chem3dhooks1_1"
#include "Chem3dprocs1_3a"
#include "FilteringMain"
#include "StoBeImportProcs"

Menu "Macros"
	"Clustering Algorithm", clusterPanel()
end

Function clusterPanelSimple()

	DoWindow ClusteringAlgorithmBasic
	if(V_Flag)
		KillWindow/Z ClusteringAlgorithmBasic
	endif
	
	NewPanel/K=1/N=ClusteringAlgorithmBasic/W=(0,0,500,500)
	
	//Section background color definitions
	Variable bc1=21588//25700
	Variable bc2=2827//0
	Variable bc3=3598//8995
	
	//Font color definition
	Variable fc1=65535
	Variable fc2=65535
	Variable fc3=65535
	
	//Load Button color definitions
	Variable lb1 = 53713
	Variable lb2 = 42148
	Variable lb3 = 13878
	
	//Set "SetVariable" valueback color
	Variable sv1 = 57825
	Variable sv2 = 53970
	Variable sv3 = 46774
	
	String atomName	 = "C"	//What kind of atom are we processing?
	Variable nAtoms	 = 1		//How many atoms to process?
	Variable osThre	 = 10		//Set a value between 0 and 100%. 
	Variable ovpThre	 = 50	//Peak overlap threshold, number between 0 and 100
	Variable rotSym	 = 1		//What's the symmetry of the molecule? Currently supported (1 = no symmetry, 2 = 2-fold rotation [180 degrees], 3 = 3 fold rotation [120 degrees], 4 = 4-fold rotation [90 degrees])
	Variable nPasses = 1
	//Optional Parameters for Clustering/Filtering Algorithm
	Variable erange = 320					//Maximum energy range to consider transitions in
	Variable broad1 = 0.5					//Peak width at lower energy (ewid1)
	Variable broad2 = 12					//Peak width at higher energy (ewid2)
	Variable ewid1  = 288						//Lower energy at which to apply broad1
	Variable ewid2	= 320 				//Higher energy at which to apply broad2
	Variable Eini 	= 280						//Starting energy at which to make transition peaks
	Variable Efin 	= 320						//Ending energy at which to make transition peaks
	Variable fwhmConversion = 2.355 		//Conversion factor between FWHM and standard deviation
	Variable broadShift = 0				//Allows for the adjustment of the energies at which to apply the broadening scheme
	Variable d = 0							//Display the data?
	Variable gRes = 2000					//Energy Resolution of the Gaussian peaks 
	Variable buildSym	= 1				//Calculate TDM tensor and determine theta,symmetry type and x,y,z components of OS? Yes or no
	Variable ovpOri		= 0 					//Determine overlap matrix for original transitions? ***Can be a bit of a lengthy process***
	Variable dClusters	= 0					//Display cluster plots from symmetry/overlap routines?
	Variable realign		= 0
	Variable thx = 0
	Variable thy = 0
	Variable thz = 0
	String rotOrder = "x,y,z"
	//Fitting of DFT optical model to experimental NEXAFS
	Variable modelAlpha 	= 0
	String expSpecName 		= "ScanID"
	String expEnergyName 	= "E"
	String expFolderPath 	= "root:Experiment"
	String mol				= "myMolecule"
	Variable holdAmps     = 0
	Variable holdWidths   = 1
	Variable holdPos      = 1
	Variable holdAlpha    = 1
	//Optional Parameters for fitting of DFT to experiment
	Variable alpha	= 54.7
	Variable i0		= 1.0E5
	Variable phi		= 0 
	Variable fitARNEXAFS	= 0
	Variable rigidShift		= 0
	Variable stepShift		= 0
	String thetaList			= "20;40;54.7;70;90"
	String NEXAFStype
	
	Variable modelOnly		= 0
	Variable fitOnly			= 0
	Variable anchorStep1	= 280
	Variable anchorStep2	= 360
	Variable anchorExp1		= 280
	Variable anchorExp2		= 340
	Variable transitionSym = 1	//What's the transition symmetry? Isotropic = 0, Uniaxial = 1 , Biaxial = 2, Triaxial = 3
	Variable clusterID = 0
	Variable pass = 1
	Variable stage = 0
	String clusterSym = "XX"
	String transSym = "XX"
	Variable nFiles = 1
	Variable stepWid1 = 0.5
	Variable stepWid2 = 12
	Variable stepE1 = 288
	Variable stepE2 = 320
	Variable refine = 0
	
	//DFT Data Loading section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00
	DrawRect 5,3,259,113
	TitleBox dftLoadingTitle title="DFT Data Loading"                   ,pos={79,8}             ,fstyle=1                         ,fColor=(fc1,fc2,fc3)
	SetVariable nFiles       title="Number of Atoms"      ,size={145,20},pos={11,37} ,value=_NUM:nFiles,limits={1,100,1},valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"How many atoms are being processed?"}
	SetVariable atomName     title="Element"              ,size={75,20} ,pos={176,36},value=_STR:atomName               ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What element are we processing?"}
	Button loadDFT           title="Load DFT Data"        ,size={200,30},pos={33,81} ,proc=loadDFTButton                ,fColor=(lb1,lb2,lb3),fSize=16,fstyle=3, appearance={os9}
	SetVariable mol          title="Molecule Name"        ,size={175,20},pos={44,61} ,value=_STR:mol                    ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3)
	SetDrawEnv fillfgc= (49151,53155,65535),linethick= 2.00 
	DrawRect 226,6,251,33
	SetDrawEnv fstyle= 1,fsize= 20
	DrawText 234,32,"1"	
End

Function clusterPanel()
	
	if(DataFolderExists("root:Packages:"))
		NewDataFolder/O/S root:Packages:DFTClustering
	else
		NewDataFolder/O/S root:Packages
		NewDataFolder/O/S root:Packages:DFTClustering
	endif
	
	DoWindow ClusteringAlgorithm
	if(V_Flag)
		KillWindow/Z ClusteringAlgorithm
	endif	
	
	//W=(left,top,right,bottom) (0,0,410,460)
	NewPanel/K=1/N=ClusteringAlgorithm/W=(680,0,1393,753)
	//SetWindow ClusteringAlgorithm graphicsTech=3 
	//MoveWindow/W=ClusteringAlgorithm 360,0,1110,420
	String atomName	 = "C"	//What kind of atom are we processing?
	Variable nAtoms	 = 1		//How many atoms to process?
	Variable osThre	 = 10		//Set a value between 0 and 100%. 
	Variable ovpThre	 = 50	//Peak overlap threshold, number between 0 and 100
	Variable rotSym	 = 1		//What's the symmetry of the molecule? Currently supported (1 = no symmetry, 2 = 2-fold rotation [180 degrees], 3 = 3 fold rotation [120 degrees], 4 = 4-fold rotation [90 degrees])
	Variable nPasses = 1
	//Optional Parameters for Clustering/Filtering Algorithm
	Variable erange = 320					//Maximum energy range to consider transitions in
	Variable broad1 = 0.5					//Peak width at lower energy (ewid1)
	Variable broad2 = 12					//Peak width at higher energy (ewid2)
	Variable ewid1  = 288						//Lower energy at which to apply broad1
	Variable ewid2	= 320 				//Higher energy at which to apply broad2
	Variable Eini 	= 280						//Starting energy at which to make transition peaks
	Variable Efin 	= 320						//Ending energy at which to make transition peaks
	Variable fwhmConversion = 2.355 		//Conversion factor between FWHM and standard deviation
	Variable broadShift = 0				//Allows for the adjustment of the energies at which to apply the broadening scheme
	Variable d = 0							//Display the data?
	Variable gRes = 2000					//Energy Resolution of the Gaussian peaks 
	Variable buildSym	= 1				//Calculate TDM tensor and determine theta,symmetry type and x,y,z components of OS? Yes or no
	Variable ovpOri		= 0 					//Determine overlap matrix for original transitions? ***Can be a bit of a lengthy process***
	Variable dClusters	= 0					//Display cluster plots from symmetry/overlap routines?
	Variable realign		= 0
	Variable thx = 0
	Variable thy = 0
	Variable thz = 0
	String rotOrder = "x,y,z"
	//Fitting of DFT optical model to experimental NEXAFS
	Variable modelAlpha 	= 0
	String expSpecName 		= "ScanID"
	String expEnergyName 	= "E"
	String expFolderPath 	= "root:Experiment"
	String mol				= "myMolecule"
	Variable holdAmps     = 0
	Variable holdWidths   = 1
	Variable holdPos      = 1
	Variable holdTensorElems = 1
	//Optional Parameters for fitting of DFT to experiment
	Variable alpha	= 54.7
	Variable i0		= 1.0E5
	Variable phi		= 0 
	Variable fitARNEXAFS	= 0
	Variable rigidShift		= 0
	Variable stepShift		= 0
	String thetaList			= "20;40;54.7;70;90"
	String NEXAFStype
	
	Variable modelOnly		= 0
	Variable fitOnly			= 0
	Variable anchorStep1	= 280
	Variable anchorStep2	= 360
	Variable anchorExp1		= 280
	Variable anchorExp2		= 340
	Variable transitionSym = 1	//What's the transition symmetry? Isotropic = 0, Uniaxial = 1 , Biaxial = 2, Triaxial = 3
	Variable clusterID = 0
	Variable transitionID = 0
	Variable pass = 1
	Variable stage = 0
	String clusterSym = "XX"
	String transSym = "XX"
	Variable nFiles = 1
	Variable stepWid1 = 0.5
	Variable stepWid2 = 12
	Variable stepE1 = 288
	Variable stepE2 = 320
	Variable refine = 0
	String clusterPltName = "TEST"
	Variable alphaFitVal = 0	
	//Fit refinement stuff
	Variable maskEnergyIni=270
	Variable maskEnergyFin=360
	Variable pkToRefine=0
		
	//Color palette used: https://coolors.co/540b0e-9e2a2b-335c67-fff3b0-e09f3e
	
	//Section background color definitions
	Variable bc1=21588//25700
	Variable bc2=2827//0
	Variable bc3=3598//8995
	
	//Font color definition
	Variable fc1=65535
	Variable fc2=65535
	Variable fc3=65535
	
	//Load Button color definitions
	Variable lb1 = 53713
	Variable lb2 = 42148
	Variable lb3 = 13878
	
	//Set "SetVariable" valueback color
	Variable sv1 = 57825
	Variable sv2 = 53970
	Variable sv3 = 46774 
	
	//Display functions button color
	Variable df1 = 13107
	Variable df2 = 23644
	Variable df3 = 26471
	
	//Run Algorithm button color
	Variable rb1 = 0
	Variable rb2 = 43176
	Variable rb3 = 30840
	
	//DFT Data Loading section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00
	DrawRect 5,3,259,113
	TitleBox dftLoadingTitle title="DFT Data Loading"                   ,pos={79,8}             ,fstyle=1                         ,fColor=(fc1,fc2,fc3)
	SetVariable nFiles       title="Number of Atoms"      ,size={145,20},pos={11,37} ,value=_NUM:nFiles,limits={1,100,1},valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"How many atoms are being processed?"}
	SetVariable atomName     title="Element"              ,size={75,20} ,pos={176,36},value=_STR:atomName               ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What element are we processing?"}
	Button loadDFT           title="Load DFT Data"        ,size={200,30},pos={33,81} ,proc=loadDFTButton                ,fColor=(lb1,lb2,lb3),fSize=16,fstyle=3, appearance={os9}
	SetVariable mol          title="Molecule Name"        ,size={175,20},pos={44,61} ,value=_STR:mol                    ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3)
	SetDrawEnv fillfgc= (49151,53155,65535),linethick= 2.00 
	DrawRect 226,6,251,33
	SetDrawEnv fstyle= 1,fsize= 20
	DrawText 234,32,"1"	
	
	//Reorient Molecular Coordinate System Section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	 
	DrawRect 264,3,555,113
	TitleBox reorienting  title="Reorient Molecular Coordinate System"                            ,pos={292,8}              ,fstyle=1                                      ,fColor=(fc1,fc2,fc3)		
	CheckBox realign      title="Reorient the Coordinate Frame of the Molecule?"                  ,pos={270,80}                            ,value = realign      ,fColor=(fc1,fc2,fc3)
	SetVariable thx       title="thx"                                            ,size={60,20}    ,pos={295,57} ,limits={0,360,1}          ,value=_NUM:thx       ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable thy       title="thy"                                            ,size={60,20}    ,pos={371,57} ,limits={0,360,1}          ,value=_NUM:thy       ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable thz       title="thz"                                            ,size={60,20}    ,pos={443,57} ,limits={0,360,1}          ,value=_NUM:thz       ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable rotOrder  title="Order of Axis Rotations"                        ,size={200,20}   ,pos={294,35}                            ,value=_STR:rotOrder  ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	
	//3D structure loader, Molecule JPEG loader, and Support Website link
	Button load3D       title="Load 3D Model\rfor Molecule from .xyz",pos={561,5}  ,size={150,50},proc=load3Dstruct      ,fColor=(lb1,lb2,lb3),fSize=13,fstyle=2		, appearance={os9}
	Button LoadStruct   title="Load Molecule Structure"              ,pos={561,59} ,size={150,20},proc=loadMolStructure	,fColor=(lb1,lb2,lb3)	, appearance={os9}
	Button DownloadLink title="Support Website"                      ,pos={561,85} ,size={150,30},proc=supportButton, fsize=12,Help={"Link to website to download analysis panel"},fcolor=(52428,1,1), appearance={os9}
	
	//Experimental Data Inputs Section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00		
	DrawRect 4,120,555,225	
	TitleBox experimentParams title="Experimental Data Inputs"     ,pos={270,122} ,fstyle=1                                                       ,fColor=(fc1,fc2,fc3)		
	PopupMenu spectraType     title="NEXAFS Type"                  ,pos={137,151}              ,value="BL;TEY;PEY"                     ,fColor=(fc1,fc2,fc3)
	SetVariable rigidShift    title="DFT-EXP Shift"                ,pos={13,151} ,size={120,20},value=_NUM:rigidShift,limits={-5,5,0.1},fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable expSpecName   title="Experiment NEXAFS Base Name"  ,pos={14,174},size={250,20},value=_STR:expSpecName                 ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable expEnergyName title="Experiment Energy Base Name"  ,pos={272,174},size={250,20},value=_STR:expEnergyName               ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable expFolderPath title="Experiment Folder Path"       ,pos={13,201},size={510,20},value=_STR:expFolderPath               ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable anchorExp1    title="Low Exp Anchor"               ,pos={249,151},size={130,20},value=_NUM:anchorExp1                  ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable anchorExp2    title="High Exp Anchor"              ,pos={385,151},size={135,20},value=_NUM:anchorExp2                  ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)

	//Main Clustering Parameters Section 
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00
	DrawRect 4,231,386,355
	TitleBox clusteringTitle title="Main Clustering Parameters"         ,pos={114,236}         ,fstyle=1                                             ,fColor=(fc1,fc2,fc3)
	SetVariable rotSym       title="Molecule Symmetry"                  ,pos={14,286}  ,size={145,20},limits={1,6,1}   ,value=_NUM:rotSym  ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What is the maximum axis of rotation of the molecule?"}
	SetVariable osThre       title="OS %"                               ,pos={17,263} ,size={75,20} ,limits={0,100,1} ,value=_NUM:osThre  ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What is the oscillator strength threshold value?"}
	SetVariable ovpThre      title="Overlap %"                          ,pos={93,263} ,size={100,20},limits={0,100,1} ,value=_NUM:ovpThre ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What is the peak overlap threshold value?"}
	SetVariable erange       title="Max Transition Energy"              ,pos={200,263} ,size={175,20}                  ,value=_NUM:erange   ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What is the maximum DFT energy to consider transitions in?"}
	//CheckBox ovpOri          title="Construct Original Overlap Matrix?" ,pos={513,313}                                 ,value = ovpOri     ,fColor=(fc1,fc2,fc3),help={"Calculate the Overlap Matrix for the initial transitions? Can take several minutes."}
	CheckBox buildSym        title="Symmetrize Transitions?"            ,pos={165,288}                                  ,value = buildSym   ,fColor=(fc1,fc2,fc3),help={"Enforce molecule symmetry?"}
	//Select Energy correction wave
	
	Button selECorrWave                                                 ,pos={204,307} ,size={150,20}
	MakeButtonIntoWSPopupButton("ClusteringAlgorithm", "selECorrWave", "notifyEcorrWaveSelection", options=PopupWS_OptionFloat)
	TitleBox WSPopupTitle1                                              ,pos={15,307} ,size={115,12},title="Select DFT Energy Correction Wave",frame=0,fColor=(fc1,fc2,fc3)
	String/G pathtoEcorr
	//Select LUMO wave
	Button selLUMOwave                                                  ,pos={204,327} ,size={150,20}
	MakeButtonIntoWSPopupButton("ClusteringAlgorithm", "selLUMOwave", "notifyLUMOWaveSelection", options=PopupWS_OptionFloat)
	TitleBox WSPopupTitle3                                              ,pos={98,327} ,size={115,12},title="Select LUMO Wave",frame=0,fColor=(fc1,fc2,fc3)
	String/G pathtoLUMO
	
	//Peak Broadening Parameters Section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00
	DrawRect 390,231,538,355			
	TitleBox BroadeningParams title="Broadening Parameters"              ,pos={399,236}           ,fstyle=1                              ,fColor=(fc1,fc2,fc3)		
	SetVariable broad1        title="Width 1"                              ,size={100,20},pos={414,261},limits={0,50,0.1}  ,value=_NUM:broad1 ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What is the peak broadening before Width Energy 1?"}
	SetVariable broad2        title="Width 2"                              ,size={100,20},pos={414,282},limits={0,50,1}    ,value=_NUM:broad2 ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"What is the peak broadening before Width Energy 2?"}
	
	//Peak Making Parameters Section	
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00
	DrawRect 544,231,708,355
	TitleBox PeakParams         title="Wave Scaling"                 ,pos={593,236}       ,fstyle=1                                            ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)	
	SetVariable Eini            title="Min Energy"   ,size={120,20},pos={561,263},limits={270,12000,10}   ,value=_NUM:Eini           ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable Efin            title="Max Energy"   ,size={120,20},pos={561,286},limits={270,12000,10}   ,value=_NUM:Efin           ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable gRes            title="Points"   ,size={90,20},pos={591,308},limits={100,10000,100},value=_NUM:gRes           ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)	
		
	//Modeling Tilt Angle Section	
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	
	DrawRect  4,361,428,571
	TitleBox modelTiltName   title="Modeling Molecular Tilt Angle"   ,pos={87,365}   ,fstyle=1                                                               ,fColor=(fc1,fc2,fc3)		
	SetVariable thetaList    title="NEXAFS θ List"                   ,pos={10,395} ,size={175,20}                  ,value=_STR:thetaList                      ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	CheckBox fitARNEXAFS     title="Fit AR NEXAFS to DFT-BB Model?"  ,pos={80,463}                                 ,value = fitARNEXAFS 							 ,fColor=(fc1,fc2,fc3),proc=CheckProc1
	CheckBox modelAlpha      title="Model alpha for DFT-BB Model?"   ,pos={80,444}                                 ,value = modelAlpha                        ,fColor=(fc1,fc2,fc3),proc=CheckProc4
	SetVariable alpha        title="α[°]"                   			,pos={269,395} ,size={75,20},limits={0,90,5}   ,value=_NUM:alpha                          ,fColor=(fc1,fc2,fc3),help={"Molecular tilt angle in degrees"},valueBackColor=(sv1,sv2,sv3)
	SetVariable i0           title="I0"                              ,pos={348,395},size={75,20}                   ,value=_NUM:i0                             ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable phi          title="φ[°]"                     			,pos={189,395},size={75,20} ,limits={0,180,45},value=_NUM:phi                            ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	CheckBox justModel       title="Model Only"                      ,pos={276,442},size={200,20}                  ,value=modelOnly                           ,fColor=(fc1,fc2,fc3),proc=CheckProc3
	CheckBox justFit         title="Fit Only"                        ,pos={276,463},size={200,20}                  ,value=fitOnly                             ,fColor=(fc1,fc2,fc3),proc=CheckProc2
	CheckBox holdAmps        title="Hold amps?"                      ,pos={206,487},size={200,20}                  ,value=holdAmps                            ,fColor=(fc1,fc2,fc3)
	CheckBox holdWidths      title="Hold widths?"                    ,pos={114,487},size={200,20}                  ,value=holdWidths                          ,fColor=(fc1,fc2,fc3)
	CheckBox holdPos         title="Hold positions?"                 ,pos={9,487},size={200,20}                  ,value=holdPos                             ,fColor=(fc1,fc2,fc3)
	CheckBox refine          title="Refine Fit?"                     ,pos={54,510},size={200,20}                  ,value=refine                              ,fColor=(fc1,fc2,fc3)
	SetVariable MaskE1       title="Mask Energy1"                   	,pos={69,418},size={140,20}                  ,value=_NUM:maskEnergyIni                  ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable MaskE2       title="Mask Energy2"                    ,pos={218,418},size={140,20}                  ,value=_NUM:maskEnergyFin                  ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable PeakToRefine title="Peak To Refine"            	   ,pos={137,510},size={120,20}                   ,value=_NUM:pkToRefine                     ,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	CheckBox holdTensorElems title="Hold modθ?"                      ,pos={288,487},size={200,20}                  ,value=holdTensorElems                      ,fColor=(fc1,fc2,fc3)
	
	//Select Refinement Paramter wave
	Button selRefPWave                                                 ,pos={275,547} ,size={150,20}
	MakeButtonIntoWSPopupButton("ClusteringAlgorithm", "selRefPWave", "notifyselRefPWaveSelection", options=PopupWS_OptionFloat)
	TitleBox WSPopupTitle4                                              ,pos={16,547} ,size={115,12},title="Select Refinement Parameter Wave",frame=0,fColor=(fc1,fc2,fc3)
	String/G pathtoRefPWave 
	
	//DFT Step Edge Section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	
	DrawRect 433,361,708,509
	TitleBox stepParamsTitle title="DFT Step Edge Parameters"                        ,pos={495,365},fColor=(fc1,fc2,fc3),fstyle=1 
	SetVariable stepWid1     title="σ\\BS1"                           ,size={75,20}  ,pos={437,392},limits={0,100,1}    ,value=_NUM:stepWid1    ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"How many atoms are being processed?"}
	SetVariable stepWid2     title="σ\\BS2"                           ,size={75,20}  ,pos={437,416},limits={0,100,1}    ,value=_NUM:stepWid2    ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"How many atoms are being processed?"}
	SetVariable stepE1       title="E\\BS1"                           ,size={75,20}  ,pos={521,392},limits={200,900,1}  ,value=_NUM:stepE1      ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"How many atoms are being processed?"}
	SetVariable stepE2       title="E\\BS2"                           ,size={75,20}  ,pos={521,416},limits={200,900,1}  ,value=_NUM:stepE2      ,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3),help={"How many atoms are being processed?"}
	Button stepBuild         title="Construct Molecule Step Edge"     ,size={250,25} ,pos={446,473}     ,proc=makeBareAtomAbs,fSize=16,fstyle=3                              ,fColor=(lb1,lb2,lb3),appearance={os9}
	Button selIPwave,pos={529,438},size={150,20}
	MakeButtonIntoWSPopupButton("ClusteringAlgorithm", "selIPwave", "notifyIPwaveSelection", options=PopupWS_OptionFloat)
	TitleBox WSPopupTitle2,pos={442,438},size={115,12},title="Select IP Wave",frame=0,fColor=(fc1,fc2,fc3)
	String/G pathtoIP
	SetVariable anchorStep1  title="Anchor\\BLO"  ,size={100,20},pos={605,392},value=_NUM:anchorStep1,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	SetVariable anchorStep2  title="Anchor\\BHI"  ,size={100,20},pos={605,416},value=_NUM:anchorStep2,fColor=(fc1,fc2,fc3),valueBackColor=(sv1,sv2,sv3)
	
	//Display fit results section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	
	DrawRect 435,513,707,571
	TitleBox fitResults title="Show Fit Results",pos={526,515},fColor=(fc1,fc2,fc3),fstyle=1 
	SetVariable bbPeak title = "BB Peak",pos={443,540},size={80,20},limits={0,inf,1},value=_NUM:0,valueBackColor=(sv1,sv2,sv3)    ,fColor=(fc1,fc2,fc3),proc=SetVarProc4
	SetVariable resAlpha title = "Alpha",pos={529,540},size={80,20},limits={0,90,1},value=_NUM:0,valueBackColor=(sv1,sv2,sv3),fColor=(fc1,fc2,fc3)
	Button fitPlot        title="Plot Fits"       ,pos={618,540} ,size={80,20}                                               ,valueColor=(65535,65535,65535)  ,fColor=(df1,df2,df3)      ,proc=displayfitResults,fSize=12,fstyle=3 , appearance={os9}

	//This button displays information about a particular cluster
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	
	DrawRect 4,575,223,746
	TitleBox clusterInfoTitle   title="Cluster Info Panel" ,pos={54,582}  ,fstyle=1                                                                                              ,fColor=(fc1,fc2,fc3)
	Button clusterInfo          title="Cluster Info"       ,pos={73,724} ,size={100,20}                                               ,valueColor=(65535,65535,65535)  ,fColor=(df1,df2,df3)      ,proc=displayClusterInfo,fSize=12,fstyle=3 , appearance={os9}
	SetVariable clusterID       title="Cluster ID"         ,pos={62,610} ,size={100,20},limits={0,150,1},value=_NUM:clusterID         ,valueBackColor=(sv1,sv2,sv3)    ,fColor=(fc1,fc1,fc1)      ,proc=ChangeClusterView
	SetVariable transitionID    title="Transition ID"      ,pos={53,635} ,size={150,20},limits={0,150,1},value=_NUM:transitionID         ,valueBackColor=(sv1,sv2,sv3)    ,fColor=(fc1,fc1,fc1)   ,proc=ChangeTransitionView
	SetVariable alphaFitVal     title="Alpha From Fit"     ,pos={53,658} ,size={130,20} ,value=_NUM:alphaFitVal
	SetVariable clusterPlotName title="Plot Name"          ,pos={23,682} ,size={160,20} ,value=_STR:clusterPltName
	
	//This button displays the overlap matrix
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	
	DrawRect 238,645,383,746
	TitleBox ovpMatTitle  title="Display Overlap Matrix"  ,pos={248,649}    ,fstyle=1                                                                                         ,fColor=(fc1,fc2,fc3)
	Button ovpMat         title="Display Matrix"          ,pos={268,724},size={100,20}                                               ,valueColor=(65535,65535,65535),fColor=(df1,df2,df3),proc=displayOVPMatButton,fSize=12,fstyle=3, appearance={os9}
	SetVariable stage     title="Algorithm Stage"         ,pos={252,679},size={120,20},limits={0,inf,1},value=_NUM:pass              ,valueBackColor=(sv1,sv2,sv3)  ,fColor=(fc1,fc2,fc3),proc=SetVarProc2
	
	//This button plots the OS for a given symmetry at the different stages of the algorithm
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00	
	DrawRect 392,645,526,746        
	TitleBox otherPlots   title="Plot OS"      ,pos={436,649}  ,fstyle=1                                                                                           ,fColor=(fc1,fc2,fc3)
	Button plotStuff      title="Display OS"   ,pos={410,724},size={100,20}                                               ,valueColor=(65535,65535,65535) ,fColor=(df1,df2,df3),proc=displayOSPlot ,fSize=12,fstyle=3  , appearance={os9} 
	SetVariable stage2    title="Merge Stage"  ,pos={400,679},size={120,20},limits={1,inf,1},value=_NUM:pass              ,valueBackColor=(sv1,sv2,sv3)   ,fColor=(fc1,fc2,fc3),proc=SetVarProc3
	
	//Explore Parameter Space Section
	SetDrawEnv fillfgc= (bc1,bc2,bc3),fillpat= 1,linejoin= 2,linethick= 3.00
	DrawRect 239,579,705,640   	
	TitleBox paramSpace   title="Explore Parameter Space"      ,pos={348,582}  ,fstyle=1  ,fColor=(fc1,fc2,fc3)
	Button selOSTWave,pos={305,616},size={150,20}
	MakeButtonIntoWSPopupButton("ClusteringAlgorithm", "selOSTWave", "notifyOSTwaveSelection", options=PopupWS_OptionFloat)
	TitleBox WSPopupTitle5,pos={245,616},size={115,12},title="OST Wave",frame=0,fColor=(fc1,fc2,fc3)
	String/G pathtoOST
	Button selOVPTWave,pos={548,616},size={150,20}
	MakeButtonIntoWSPopupButton("ClusteringAlgorithm", "selOVPTwave", "notifyOVPTwaveSelection", options=PopupWS_OptionFloat)
	TitleBox WSPopupTitle6,pos={482,616},size={115,12},title="OVPT Wave",frame=0,fColor=(fc1,fc2,fc3)
	String/G pathtoOVPT
	CheckBox expParamSpace title="Explore Parameter Space?"      ,pos={550,588},size={200,20}                  ,value=0    ,fColor=(fc1,fc2,fc3)

	//This is the Run Algorithm button
	Button run title="Run \rAlgorithm",pos={532,643},size={179,104},proc=RunAlgorithmButton,fColor=(rb1,rb2,rb3),fSize=32,fstyle=3, appearance={os9}
	
	//Display WSU logo in bottom left corner of panel
	DisplayWSULogo()	
End

Function runAlgorithmButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			String startFolder = GetdataFolder(1)
			ControlInfo/W=ClusteringAlgorithm atomName
			String atomName = S_Value
			
			ControlInfo/W=ClusteringAlgorithm nFiles
			Variable nAtoms = V_Value
			
			ControlInfo/W=ClusteringAlgorithm rotSym
			Variable rotSym = V_Value
			
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			
			ControlInfo/W=ClusteringAlgorithm erange
			Variable erange = V_Value			
			
			ControlInfo/W=ClusteringAlgorithm buildSym
			Variable buildSym = V_Value
			String bSym
			if(buildSym)
				bSym = "yes"
			else
				bSym = "no"
			endif
			
			ControlInfo/W=ClusteringAlgorithm broad1
			Variable broad1 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm broad2
			Variable broad2 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm Eini
			Variable Eini = V_Value
			
			ControlInfo/W=ClusteringAlgorithm Efin
			Variable Efin = V_Value
			
			ControlInfo/W=ClusteringAlgorithm gRes
			Variable gRes = V_Value
			
			Wave Ecorr = $GetEcorrPath()
			
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			
			ControlInfo/W=ClusteringAlgorithm spectraType
			String spectraType = S_Value
			
			ControlInfo/W=ClusteringAlgorithm rigidShift
			Variable rigidShift = V_Value
			
			ControlInfo/W=ClusteringAlgorithm expSpecName
			String expSpecName = S_Value
			
			ControlInfo/W=ClusteringAlgorithm expEnergyName
			String expEnergyName = S_Value
			
			ControlInfo/W=ClusteringAlgorithm expFolderPath
			String expFolderPath = S_Value
			
			ControlInfo/W=ClusteringAlgorithm fitARNEXAFS
			Variable fitARNEXAFS = V_Value
			String fitExpNEXAFS
			if(fitARNEXAFS)
				fitExpNEXAFS = "yes"
			else
				fitExpNEXAFS = "no"
			endif
			
			ControlInfo/W=ClusteringAlgorithm alpha
			Variable alpha = V_Value
			
			ControlInfo/W=ClusteringAlgorithm i0
			Variable i0 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm phi
			Variable phi = V_Value
			
			ControlInfo/W=ClusteringAlgorithm thetaList
			String thetaList = S_Value	
			
			Wave IP    = $GetIPPath()
			
			ControlInfo/W=ClusteringAlgorithm realign
			Variable realign = V_Value
			
			ControlInfo/W=ClusteringAlgorithm thx
			Variable thx = V_Value
			
			ControlInfo/W=ClusteringAlgorithm thy
			Variable thy = V_Value
			
			ControlInfo/W=ClusteringAlgorithm thz
			Variable thz = V_Value
			
			ControlInfo/W=ClusteringAlgorithm rotOrder
			String rotOrder = S_Value	
	
			ControlInfo/W=ClusteringAlgorithm modelAlpha
			Variable modelAlpha = V_Value
			String modAlpha
			if(modelAlpha)
				modAlpha = "yes"
			else
				modAlpha = "no"
			endif
			
			ControlInfo/W=ClusteringAlgorithm justModel
			Variable justModel = V_Value
			
			ControlInfo/W=ClusteringAlgorithm justFit
			Variable justFit = V_Value
					
			Wave LUMOwave = $GetLUMOPath()
			
			ControlInfo/W=ClusteringAlgorithm anchorStep1
			Variable anchorStep1 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm anchorStep2
			Variable anchorStep2 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm anchorExp1
			Variable anchorExp1 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm anchorExp2
			Variable anchorExp2 = V_Value
				
			ControlInfo/W=ClusteringAlgorithm stepWid1
			Variable stepWid1 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm stepWid2
			Variable stepWid2 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm stepE1
			Variable stepE1 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm stepE2
			Variable stepE2 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm holdAmps
			Variable holdAmps = V_Value
			
			ControlInfo/W=ClusteringAlgorithm holdPos
			Variable holdPos = V_Value
			
			ControlInfo/W=ClusteringAlgorithm holdWidths
			Variable holdWidths = V_Value
			
			ControlInfo/W=ClusteringAlgorithm refine
			Variable refine = V_Value
			
			ControlInfo/W=ClusteringAlgorithm MaskE1
			Variable maskEnergy1 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm MaskE2
			Variable maskEnergy2 = V_Value
			
			ControlInfo/W=ClusteringAlgorithm PeakToRefine
			Variable pkToRefine = V_Value
			
			Wave refitPWave = $GetrefitPWavePath()
			
			ControlInfo/W=ClusteringAlgorithm holdTensorElems
			Variable holdTensorElems = V_Value
			
			ControlInfo/W=ClusteringAlgorithm expParamSpace
			Variable paramSpaceExplore = V_Value
			
			if(!paramSpaceExplore)
				filterDFT(atomName,nAtoms,osThre,ovpThre,rotSym,LUMOwave,erange=erange,broad1=broad1,broad2=broad2,Eini=Eini,Efin=Efin,gRes=gRes,buildSym=bSym,corrWave=Ecorr,realign=realign,thx=thx,thy=thy,thz=thz,rotOrder=rotOrder,modelAlpha=modAlpha,IPwave=IP,expSpecName=expSpecName,expEnergyName=expEnergyName,expFolderPath=expFolderPath,mol=mol,alpha=alpha,i0=i0,phi=phi,fit=fitExpNEXAFS,rigidShift=rigidShift,thetaList=thetaList,NEXAFStype=spectraType,justModel=justModel,justFit=justFit,anchorStep1=anchorStep1,anchorStep2=anchorStep2,anchorExp1=anchorExp1,anchorExp2=anchorExp2,stepWid1=stepWid1,stepWid2=stepWid2,stepE1=stepE1,stepE2=stepE2,holdAmps=holdAmps,holdWidths=holdWidths,holdPos=holdPos,refinement=refine,maskEnergy1=maskEnergy1,maskEnergy2=maskEnergy2,pkToRefine=pkToRefine,refitPWave=refitPWave,holdTensorElems=holdTensorElems)
				SetDataFolder $startFolder
				break
			else
				Wave osWave  = $GetOSTPath()
				Wave ovpWave = $GetOVPTPath()	
				seqThresholds(atomName,nAtoms,osWave,ovpWave,rotSym,LUMOwave,erange=erange,broad1=broad1,broad2=broad2,Eini=Eini,Efin=Efin,gRes=gRes,buildSym=bSym,corrWave=Ecorr,realign=realign,thx=thx,thy=thy,thz=thz,rotOrder=rotOrder,modelAlpha=modAlpha,IPwave=IP,expSpecName=expSpecName,expEnergyName=expEnergyName,expFolderPath=expFolderPath,mol=mol,alpha=alpha,i0=i0,phi=phi,fit=fitExpNEXAFS,rigidShift=rigidShift,thetaList=thetaList,NEXAFStype=spectraType,justModel=justModel,justFit=justFit,anchorStep1=anchorStep1,anchorStep2=anchorStep2,anchorExp1=anchorExp1,anchorExp2=anchorExp2,stepWid1=stepWid1,stepWid2=stepWid2,stepE1=stepE1,stepE2=stepE2,holdAmps=holdAmps,holdWidths=holdWidths,holdPos=holdPos,refinement=refine,maskEnergy1=maskEnergy1,maskEnergy2=maskEnergy2,pkToRefine=pkToRefine,refitPWave=refitPWave,holdTensorElems=holdTensorElems)				
			endif
			case -1: // control being killed
				break
	endswitch

	return 0
End

Function supportButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
		if(!cmpstr(ba.ctrlname, "DownloadLink"))
			BrowseURL/Z "https://labs.wsu.edu/carbon/xray-analysis-tools/"
		endif
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function loadDFTButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm nFiles
			Variable nFiles = V_Value
			ControlInfo/W=ClusteringAlgorithm atomName
			String atomName = S_Value
			Wave Ecorr = $GetEcorrPath()
			if(WaveExists(Ecorr))	//Energy correction wave was made prior to loading the data
				PRINT "Energy correction wave provided."
				DFTwrapper1("",mol,atomName,nFiles,d=1,Ecorr=Ecorr)
			else //Energy correction wave will be made as part of the data loading process
				PRINT "Energy correction wave will be made."
				DFTwrapper1("",mol,atomName,nFiles,d=1)
			endif
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function loadMolStructure(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
		
			DoWindow/W=moleculeHost ClusteringAlgorithm
			
			if(V_Flag)
				KillWindow/Z moleculeHost
				DoWindow ClusteringAlgorithm
				Button LoadStruct title="Load Molecule Structure"
			else
				DoWindow ClusteringAlgorithm
				Button LoadStruct title="Hide Molecule Structure"
				
				DoWindow/W=moleculeHost ClusteringAlgorithm
				ControlInfo/W=ClusteringAlgorithm mol
				Wave molStruct = $S_Value
				String mol = S_Value
				ControlInfo/W=ClusteringAlgorithm tval
				Variable tval = V_Value
				ControlInfo/W=ClusteringAlgorithm ovpVal
				Variable ovpMax = V_Value
				String dataFName = "root:Packages:DFTClustering:PolarAngles_" + mol +":TransitionFiltering_" + replacestring(".",num2str(tval),"p") + "OS_" + replacestring(".",num2str(ovpMax),"p") + "OVP:"

				//Check if the molecular structure image exists. If it doesn't, load the image. Else, open the panel
				if(!WaveExists(molStruct))
					ImageLoad/T=jpeg/Q/N=$S_Value/O
				endif
				Variable picSizeX = DimSize(molStruct,0)
				Variable picSizeY = DimSize(molStruct,1)
				//Determine optimal coordinates for display of image structure while preserving aspect ratio
				if(picSizeY > 560)
					Variable yDiff = picSizeY - 560
					picSizeY = 560
					picSizeX = picSizeX - yDiff
				elseif(picSizeY < 560)
					yDiff =  560 - picSizeY
					picSizeY = 560
					picSizeX = picSizeX + yDiff
				endif
				//print picSizeX,picSizeY,yDiff
				NewPanel/HOST=ClusteringAlgorithm/EXT=0/N=$("moleculeHost")/K=1/W=(0,0,picSizeX,picSizeY) as "Molecule Structure"
				Display/W=(0,0,picSizeX,picSizeY)/N=WAT/HOST=#
				AppendImage/T/W=$("moleculeHost#WAT") molStruct //HOST=#/S=0 CGHJ
				ModifyGraph noLabel=2,axThick=0,standoff=0,margin=1
				SetAxis/A/R left
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function makeBareAtomAbs(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			makeBAWrap()
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CheckProcFitExp(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			if(checked)
				print "Please provide the location of experimental data,the base name for the spectra and the base name for the energy waves."
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function notifyEcorrWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:Packages:DFTClustering
	String/G  pathtoEcorr = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyselRefPWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:Packages:DFTClustering
	String/G  pathtoRefPWave = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyIPwaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:Packages:DFTClustering
	String/G  pathtoIP = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyLUMOWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:Packages:DFTClustering
	String/G  pathtoLUMO = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyOSTwaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:Packages:DFTClustering
	String/G  pathtoOST = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyOVPTwaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:Packages:DFTClustering
	String/G  pathtoOVPT = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

function/S GetEcorrPath()
   
   SVAR pathtoEcorr = root:Packages:DFTClustering:pathtoEcorr
   if(SVAR_Exists(pathtoEcorr))
	   	return pathtoEcorr
	else
		print "String pathtoEcorr doesn't exist. Reselect DFT Energy Correction wave."
	endif
end

function/S GetIPPath()
   
   SVAR pathtoIP = root:Packages:DFTClustering:pathtoIP
   if(SVAR_Exists(pathtoIP))
	   	return pathtoIP
	else
		print "String pathtoIP doesn't exist. Reselect IP wave."
	endif
end

function/S GetLUMOPath()
   
   SVAR pathtoLUMO = root:Packages:DFTClustering:pathtoLUMO
   if(SVAR_Exists(pathtoLUMO))
	   	return pathtoLUMO
	else
		print "String pathtoLUMO doesn't exist. Reselect LUMO wave."
	endif
end

function/S GetrefitPWavePath()
   
   SVAR pathtoRefPWave = root:Packages:DFTClustering:pathtoRefPWave
   if(SVAR_Exists(pathtoRefPWave))
	   	return pathtoRefPWave
	else
		print "String pathtoRefPWave doesn't exist. Reselect Refinement Parameter wave."
	endif
end

function/S GetOSTPath()
   
   SVAR pathtoOST = root:Packages:DFTClustering:pathtoOST
   if(SVAR_Exists(pathtoOST))
	   	return pathtoOST
	else
		print "String pathtoOST doesn't exist. Reselect OST wave."
	endif
end

function/S GetOVPTPath()
   
   SVAR pathtoOVPT = root:Packages:DFTClustering:pathtoOVPT
   if(SVAR_Exists(pathtoOVPT))
	   	return pathtoOVPT
	else
		print "String pathtoOVPT doesn't exist. Reselect OVPT wave."
	endif
end

Function DisplayWSULogo()

	String LogoPath = SpecialDirPath("Igor Pro User files",0,0,0) + "User Procedures:DFT Clustering:"
	NewPath/q/z/o Path2Logo, LogoPath
	if(!waveexists(root:Packages:DFTClustering:WSUlogo))
		ImageLoad/Z/P=Path2Logo/T=any/Q/N=WSUlogo/O "wsu_2.svg"//"wsu_black.jpg"//"index.png"//"wsu-signature-default.svg"
		if(!V_Flag)
			return 0 //No logo found, no hard feelings
		endif
	endif
	
	Display/W=(559,120,710,225)/N=Logo/HOST=ClusteringAlgorithm
	//Display/W=(9,435,175,555)/N=Logo/HOST=ClusteringAlgorithm
	AppendImage/T/W=$("ClusteringAlgorithm#Logo") root:Packages:DFTClustering:WSUlogo
	ModifyImage/W=$("ClusteringAlgorithm#Logo") WSUlogo ctab= {*,*,Grays,0}
	ModifyGraph/W=$("ClusteringAlgorithm#Logo") margin(left)=1,margin(bottom)=1,margin(top)=1,margin(right)=1
	ModifyGraph/W=$("ClusteringAlgorithm#Logo") mirror=0,nticks=0,standoff=0,axthick=0,btlen=3
	ModifyGraph/W=$("ClusteringAlgorithm#Logo") tkLblRot(left)=90
	SetAxis/W=$("ClusteringAlgorithm#Logo")/A/R left
	
End

Function load3Dstruct(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=clusteringAlgorithm mol
			String molName = S_Value
			SetDataFolder root: 
			loadAtomCoords(molName,"")
			Chem3Dmodule#Chem3Dstart()
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function loadAtomCoords(molName,pathName)

	String molName	//Desired name of file
	String pathName	//Symbolic path where desired file is present
	String columnInfoStr = " "									//Contains set of names for each column in the .out file
	columnInfoStr += "C=1,F=-2,W=7,N=Element;"
	columnInfoStr += "C=1,F=0,W=16,N=atomX;"
	columnInfoStr += "C=1,F=0,W=16,N=atomY;"
	columnInfoStr += "C=1,F=0,W=16,N=atomZ;"
	LoadWave/F={4,9,0}/B=columnInfoStr/D/K=0/O/Q/A/L={0, 2, 0, 0, 4 }/P=$pathName molName
	
	//Make2D wave with x,y,z coordinates for molecule
	String wName = "atomicCoords_" + molName
	Wave atomicCoordinates = $wName
	if(WaveExists(atomicCoordinates))
		KillWaves/Z atomicCoordinates
	endif
	
	String aCoordList = WaveLIst("atom*",";","")	
	Concatenate/O aCoordList, $wName
	
	String wName2 = "elementNames_" + molName
	Wave/T Element
	Duplicate/O Element, $wName2
	KillWaves/Z Element,atomX,atomY,atomZ
	//endif	
End

Function displayClusterInfo(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=ClusteringAlgorithm mol
			String molName = S_Value
			ControlInfo/W=ClusteringAlgorithm clusterPlotName
			String mName = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable OS = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable OVP = V_Value
			ControlInfo/W=ClusteringAlgorithm clusterID
			Variable cl = V_Value
			ControlInfo/W=ClusteringAlgorithm rigidShift
			Variable corr = V_Value
			ControlInfo/W=ClusteringAlgorithm alphaFitVal
			Variable alpha = V_Value
			tdmInfo(mName,OS,OVP,alpha,molName,corr=corr,cl=cl)
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ChangeClusterView(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			ControlInfo/W=ClusteringAlgorithm mol
			String molName = S_Value
			ControlInfo/W=ClusteringAlgorithm clusterPlotName
			String mName = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable OS = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable OVP = V_Value
			ControlInfo/W=ClusteringAlgorithm clusterID
			Variable cl = V_Value
			ControlInfo/W=ClusteringAlgorithm rigidShift
			Variable corr = V_Value
			ControlInfo/W=ClusteringAlgorithm alphaFitVal
			Variable alpha = V_Value
			tdmInfo(mName,OS,OVP,alpha,molName,corr=corr,cl=cl)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function ChangeTransitionView(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			ControlInfo/W=ClusteringAlgorithm mol
			String molName = S_Value
			ControlInfo/W=ClusteringAlgorithm clusterPlotName
			String mName = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable OS = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable OVP = V_Value
			ControlInfo/W=ClusteringAlgorithm transitionID
			Variable tr = V_Value
			ControlInfo/W=ClusteringAlgorithm rigidShift
			Variable corr = V_Value
			ControlInfo/W=ClusteringAlgorithm alphaFitVal
			Variable alpha = V_Value
			tdmInfo(mName,OS,OVP,alpha,molName,corr=corr,tr=tr)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function makeOVPMatPlots(tval,ovpMax,mol,stage)

	Variable tval
	Variable ovpMax
	String mol
	Variable stage
	
	//Define initial clustering run folder
	String iniFName = "root:Packages:DFTClustering:PolarAngles_" + mol +":TransitionFiltering_" + replacestring(".",num2str(tval),"p") + "OS_" + replacestring(".",num2str(ovpMax),"p") + "OVP:"
	String ovpFolder, folderLoc , ovpWaveName
	String stageName
	
	if(stage == 0)
		ovpFolder = "originalPeaks:"
		folderLoc = iniFName + ovpFolder
		ovpWaveName = "ovpWaveAll"
		stageName = "\\f01Unfiltered"
	elseif(stage == 1)
		ovpFolder = "firstFilter:"
		folderLoc = iniFName + ovpFolder
		ovpWaveName = "ovpOriginal"
		stageName = "\\f01" + num2str(tval) + "% OS Filtered"
	elseif(stage >1)
		ovpFolder = "firstFilter:"
		folderLoc = iniFName + ovpFolder
		ovpWaveName = "ovpWaveALL" + num2str(stage-1)
		stageName = "\\f01Merge Stage #" + num2str(stage-1)
	endif
	
	if(!DataFolderExists(folderLoc))
		Abort "Folder containing overlap matrix doesn't exist for specified clustering parameters."
	endif
	
	SetDataFolder $folderLoc
	Wave ovpWave = $ovpWaveName 
	
	if(!WaveExists(ovpWave))
		Abort "This system does not have an overlap matrix information for this symmetry at this stage."
	endif
	
	//Plot clusters
	String plotName = "Overlap_Matrix_Display"
	String tableName = "Overlap_Matrix_Values"
	DoWindow/F $plotName
	Variable n = 200,nTrans = DimSize(ovpWave,0)
	
	if(!V_Flag)
		NewImage/N=$plotName/K=1 ovpWave
		ModifyGraph/W=$plotName width  = n , height = n , margin(right)=108, margin(left)=36, margin(top)=36
		ModifyImage $ovpWaveName  ctab= {*,*,ColdWarm,0}
		ColorScale/C/N=text0/A=RC/E image=$ovpWaveName,fsize=14,fstyle=1,nticks=10
		ColorScale/C/N=text0 "\\Z14\\f01Peak Overlap[%]"
		ModifyGraph fStyle=1,fSize=14,lblPosMode=1 , nticks=10,minor=0 
		if(nTrans <=10)
			ModifyGraph manTick(left)={0,1,0,0},manMinor(left)={0,0},manTick(top)={0,1,0,0},manMinor(top)={0,0}
		endif
		Label left "Peak ID"
		Label top "Peak ID"
		TextBox/C/N=text1/A=RB/E "\\f01\\JC" + num2str(nTrans) + "\rTransitions"
		MoveWindow/W=$plotName 5,0,200,200
		TextBox/C/N=text2/E stageName
		//Edit/N=$tableName ovpWave
	else
		KillWindow/Z $plotName
		KillWindow/Z $tableName
		NewImage/N=$plotName/K=1 ovpWave
		ModifyGraph/W=$plotName width  = n , height = n , margin(right)=108, margin(left)=36, margin(top)=36
		ModifyImage $ovpWaveName  ctab= {*,*,ColdWarm,0}
		ColorScale/C/N=text0/A=RC/E image=$ovpWaveName,fsize=14,fstyle=1,nticks=10
		ColorScale/C/N=text0 "\\Z14\\f01Peak Overlap[%]"
		ModifyGraph fStyle=1,fSize=14,lblPosMode=1 , nticks=10,minor=0
		if(nTrans <=10)
			ModifyGraph manTick(left)={0,1,0,0},manMinor(left)={0,0},manTick(top)={0,1,0,0},manMinor(top)={0,0}
		endif
		Label left "Peak ID"
		Label top "Peak ID"
		TextBox/C/N=text1/A=RB/E "\\f01\\JC" + num2str(nTrans) + "\rTransitions"
		MoveWindow/W=$plotName 5,0,200,200
		TextBox/C/N=text2/E stageName
		//Edit/N=$tableName ovpWave
	endif
End

Function displayOVPMatButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			ControlInfo/W=ClusteringAlgorithm stage
			Variable stage = V_Value
			makeOVPMatPlots(osThre,ovpThre,mol,stage)
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc2(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			ControlInfo/W=ClusteringAlgorithm stage
			Variable stage = V_Value
			makeOVPMatPlots(osThre,ovpThre,mol,stage)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function PopMenuProc2(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			ControlInfo/W=ClusteringAlgorithm stage
			Variable stage = V_Value
			makeOVPMatPlots(osThre,ovpThre,mol,stage)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function makeOSPlotComparison(tval,ovpMax,mol)
	
	Variable tval
	Variable ovpMax
	String mol
	
	//Find data folders containing relevant waves
	String baseFolderName = "root:Packages:DFTClustering:PolarAngles_" + mol + ":TransitionFiltering_" + replacestring(".",num2str(tval),"p") + "OS_" + replacestring(".",num2str(ovpMax),"p") + "OVP:"
	SetDataFolder $baseFolderName
	String pathToOriOS   = baseFolderName + "originalPeaks:"
	String pathToFFOS    = baseFolderName + "firstFilter:"
	String pathToMergeOS = baseFolderName + "AmplitudeFitting:ALL:"
	//Plot the oscillator strengths for each of the filtering stages
	
	String oriName = pathToOriOS   + "allParamsSorted"
	String ffName  = pathToFFOS    + "allParamsF1Sorted"
	String mName   = pathToMergeOS + "pw2DOriginal"
	
	Wave allParamsSorted   = $oriName
	Wave allParamsF1Sorted = $ffName
	Wave mergedPWave2D     = $mName
	
	Variable n = DimSize(allParamsSorted,0),n1 = DimSize(allParamsF1Sorted,0),n2 = DimSize(mergedPWave2D,0)
	String pathToMinOS = baseFolderName + "minOSThreshold"
	NVAR minOSthreshold = $pathToMinOS
	Make/O/N=(n) OSthreshold = minOSthreshold//0.00017686
	SetScale/i x,280,320,OSthreshold
	String mName2 = "pw2DOriginal"
	String plotName = "OscillatorStrengthsComparison"
	Variable col=1
					
	DoWindow $plotName
	if(!V_Flag)
		Display/N=$plotName/K=1/W=(0,0,400,524) allParamsSorted[*][col] vs allParamsSorted[*][0]
		AppendToGraph/W=$plotName OSthreshold
		NewFreeAxis/L f1
		NewFreeAxis/L f2
		AppendtoGraph/L=f1/W=$plotName allParamsF1Sorted[*][col] vs allParamsF1Sorted[*][0]
		AppendtoGraph/L=f1/W=$plotName OSthreshold
		AppendtoGraph/L=f2/W=$plotName mergedPWave2D[*][1] vs mergedPWave2D[*][0]
		AppendtoGraph/L=f2/W=$plotName OSthreshold
		Label left "ƒ\\U"
		Label bottom "Energy[eV]"
		Label f1 "ƒ\\U"
		Label f2 "ƒ\\U"
		ModifyGraph/W=$plotName mode=1,tick(left)=2,tick(f1)=2,zero(f1)=1,mirror=1,axisEnab(left)={0,0.30},axisEnab(f1)={0.34,0.64},freePos(f1)=0,lblPosMode(f1)=2
		ModifyGraph/W=$plotName zero(f2)=1,axisEnab(f2)={0.7,1},freePos(f2)=0,lblPosMode(f2)=2,lsize=1.5
		ModifyGraph/W=$plotName rgb(allParamsF1Sorted)=(0,0,0)
		ModifyGraph/W=$plotName rgb($mName2)=(1,39321,19939),margin(top)=90,tick(bottom)=1
		ModifyGraph mode(OSthreshold)=0,lstyle(OSthreshold)=3,lsize(OSthreshold)=2,rgb(OSthreshold)=(0,0,0)
		ModifyGraph mode(OSthreshold#1)=0,lstyle(OSthreshold#1)=3,lsize(OSthreshold#1)=2,rgb(OSthreshold#1)=(0,0,0)
		ModifyGraph mode(OSthreshold#2)=0,lstyle(OSthreshold#2)=3,lsize(OSthreshold#2)=2,rgb(OSthreshold#2)=(0,0,0)
		ModifyGraph rgb(allParamsF1Sorted)=(1,26221,39321)		
		ModifyGraph grid(left)=2,grid(f1)=2,grid(f2)=2,log(left)=1,log(f1)=1,log(f2)=1,tick(f2)=2,fStyle=1,fSize=11,lblPosMode(left)=1,lblPosMode(f1)=1,lblPosMode(f2)=1
		SetAxis left 0.0001,1
		SetAxis f1 0.0001,1
		SetAxis f2 0.0001,1
		SetAxis bottom 284,320
		String legendText = "\\JC\\JL\\s(allParamsSorted) Unfiltered	["+ num2str(n) +"]\r\\s(allParamsF1Sorted) Filtered		["+ num2str(n1) +"]\r\\s("+mName2+") Clustered	["+ num2str(n2) +"]"
		Legend/C/N=text0/A=MC/E=2/X=0/Y=41 legendText
	else
		Variable i		
		KillWindow/Z $plotName
		Display/N=$plotName/K=1/W=(0,0,400,524) allParamsSorted[*][col] vs allParamsSorted[*][0]
		AppendToGraph/W=$plotName OSthreshold
		NewFreeAxis/L f1
		NewFreeAxis/L f2
		AppendtoGraph/L=f1/W=$plotName allParamsF1Sorted[*][col] vs allParamsF1Sorted[*][0]
		AppendtoGraph/L=f1/W=$plotName OSthreshold
		AppendtoGraph/L=f2/W=$plotName mergedPWave2D[*][1] vs mergedPWave2D[*][0]
		AppendtoGraph/L=f2/W=$plotName OSthreshold
		Label left "ƒ\\U"
		Label bottom "Energy[eV]"
		Label f1 "ƒ\\U"
		Label f2 "ƒ\\U"
		ModifyGraph/W=$plotName mode=1,tick(left)=2,tick(f1)=2,zero(f1)=1,mirror=1,axisEnab(left)={0,0.30},axisEnab(f1)={0.34,0.64},freePos(f1)=0,lblPosMode(f1)=2
		ModifyGraph/W=$plotName zero(f2)=1,axisEnab(f2)={0.7,1},freePos(f2)=0,lblPosMode(f2)=2,lsize=1.5
		ModifyGraph/W=$plotName rgb(allParamsF1Sorted)=(0,0,0)
		ModifyGraph/W=$plotName rgb($mName2)=(1,39321,19939),margin(top)=90,tick(bottom)=1
		ModifyGraph mode(OSthreshold)=0,lstyle(OSthreshold)=3,lsize(OSthreshold)=2,rgb(OSthreshold)=(0,0,0)
		ModifyGraph mode(OSthreshold#1)=0,lstyle(OSthreshold#1)=3,lsize(OSthreshold#1)=2,rgb(OSthreshold#1)=(0,0,0)
		ModifyGraph mode(OSthreshold#2)=0,lstyle(OSthreshold#2)=3,lsize(OSthreshold#2)=2,rgb(OSthreshold#2)=(0,0,0)
		ModifyGraph rgb(allParamsF1Sorted)=(1,26221,39321)		
		ModifyGraph grid(left)=2,grid(f1)=2,grid(f2)=2,log(left)=1,log(f1)=1,log(f2)=1,tick(f2)=2,fStyle=1,fSize=11,lblPosMode(left)=1,lblPosMode(f1)=1,lblPosMode(f2)=1
		SetAxis left 0.0001,1
		SetAxis f1 0.0001,1
		SetAxis f2 0.0001,1
		SetAxis bottom 284,320
		legendText = "\\JC\\JL\\s(allParamsSorted) Unfiltered	["+ num2str(n) +"]\r\\s(allParamsF1Sorted) Filtered		["+ num2str(n1) +"]\r\\s("+mName2+") Clustered	["+ num2str(n2) +"]"		
		Legend/C/N=text0/A=MC/E=2/X=0/Y=41 legendText
	endif
End

Function displayOSPlot(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			makeOSPlotComparison(osThre,ovpThre,mol)
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc3(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			makeOSPlotComparison(osThre,ovpThre,mol)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetVarProc4(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			ControlInfo/W=ClusteringAlgorithm alpha
			Variable alpha = V_Value
			ControlInfo/W=ClusteringAlgorithm transSym2
			String transSym2 = S_Value
			ControlInfo/W=ClusteringAlgorithm stage3
			Variable stage3 = V_Value
			ControlInfo/W=ClusteringAlgorithm thetaList
			String thetaList = S_Value
			ControlInfo/W=ClusteringAlgorithm BBpeak
			Variable whichPk = V_Value
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			displayDFTBBmodel(osThre,ovpThre,alpha,stage3,thetaList,whichPk,mol)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function PopMenuProc3(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			makeOSPlotComparison(osThre,ovpThre,mol)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CheckProc1(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Variable justModel,justFit,modelAlpha,fitARNEXAFS
			ControlInfo/W=ClusteringAlgorithm fitARNEXAFS
			fitARNEXAFS = V_Value
			if(fitARNEXAFS==1)
				CheckBox modelAlpha win=ClusteringAlgorithm,value=1
				CheckBox justModel  win=ClusteringAlgorithm,value=0
				CheckBox justFit    win=ClusteringAlgorithm,value=0
				print "Please provide the location of experimental data,the base name for the spectra and the base name for the energy waves."
			endif

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CheckProc2(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Variable justModel,justFit,modelAlpha,fitARNEXAFS
			ControlInfo/W=ClusteringAlgorithm justFit
			justFit = V_Value
			if(justFit)
				CheckBox justModel   win=ClusteringAlgorithm,value=0
				CheckBox modelAlpha  win=ClusteringAlgorithm,value=0
				CheckBox fitARNEXAFS win=ClusteringAlgorithm,value=0
			endif

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CheckProc3(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Variable justModel,justFit,modelAlpha,fitARNEXAFS
			ControlInfo/W=ClusteringAlgorithm justModel
			justModel = V_Value
			if(justModel)
				CheckBox justFit     win=ClusteringAlgorithm,value=0
				CheckBox modelAlpha  win=ClusteringAlgorithm,value=0
				CheckBox fitARNEXAFS win=ClusteringAlgorithm,value=0
			endif

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function CheckProc4(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			Variable justModel,justFit,modelAlpha,fitARNEXAFS
			ControlInfo/W=ClusteringAlgorithm modelAlpha
			modelAlpha = V_Value
			if(modelAlpha)
				CheckBox justFit     win=ClusteringAlgorithm,value=0
				CheckBox justModel   win=ClusteringAlgorithm,value=0
				//CheckBox fitARNEXAFS win=ClusteringAlgorithm,value=0
			endif

			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function displayDFTBBmodel(tval,ovpVal,alpha,whichMerge,thetaList,whichPk,mol)
	
	Variable tval,ovpVal,alpha,whichMerge,whichPk
	String thetaList,mol
	
	String filterFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":TransitionFiltering_"+replaceString(".",num2str(tval),"p")+"OS_"+replaceString(".",num2str(ovpVal),"p")+"OVP:"
	String ampFitFolder = "AmplitudeFitting:"
	String alphaFolder  = "Alpha_" + replaceString(".",num2str(alpha),"p") + ":"
	String baseFolder = filterFolder + ampFitFolder + alphaFolder
	String miscFolder = baseFolder + "TensorMisc:"
	
	Variable nItems = ItemsInList(thetaList,";"),i,j,magicAngleSpec = findMagicAngle(thetaList)
	SetDataFolder $miscFolder
	DoWindow DFT_BB_FIT
	if(!V_Flag)
		Display/K=1/N=DFT_BB_FIT/W=(0,0,600,400) 
		//Append Experimental NEXAFS and fits to Graph
		for(i=0;i<nItems;i+=1)
			String expSpecList = WaveList("expSpec*",";","")
			String cSpecWave = StringFromList(i,expSpecList,";")
			Wave w = $cSpecWave
			String expEList = WaveList("expEnergy*",";","")
			String cEWave = StringFromList(i,expEList,";")
			Wave y = $cEWave
			String fitList = WaveList("fitresults*",";","")
			String cFitWave = StringFromList(i,fitList,";")
			Wave z = $cFitWave
			AppendToGraph w,z vs y 
			ModifyGraph lstyle($cFitWave)=3
		endfor
		ApplyColorTableToTopGraph("ColdWarm")
		
		//Append Step Edge to Graph
		String stepList = WaveList("dftStep*",";","")
		String cStepWave = StringFromList(magicAngleSpec,stepList,";")
		Wave step = $cStepWave
		AppendToGraph step vs y 
		ModifyGraph lsize=1.5
		//Append BB model to Graph
		String pkList = WaveList("pk*"+num2str(magicAngleSpec),";","") 
		Variable nPks = ItemsInList(pkList)
		NewFreeAxis/L BB
		ModifyGraph axisEnab(left)={0.35,1},axisEnab(BB)={0,0.3},freePos(BB)=0,lblPosMode(BB)=1
		for(i=0;i<nPks;i+=1)
			String cPk = StringFromList(i,pkList,";")
			Wave pk = $cPk
			AppendToGraph/L=BB pk
			ModifyGraph rgb($cPk)=(52428,52428,52428)
		endfor
		SetAxis bottom 280,320
		SetAxis BB 0,50000
		ModifyGraph tick=2,mirror=1,minor=1,fStyle=1,fSize=14,margin(left)=50
		Label left "Mass Absorbance [cm\\S2\\M/g\\U"
		Label BB "BB\\U"
		Label bottom "Photon Energy[eV]"
	endif
	
	//Get the 3D tensor info for sym elements
	Wave normTensor3D,fitEns
	if(whichPk >= numpnts(fitEns))
		Abort "No such peak exists in BB model."
	endif
	Variable pkEn = fitEns[whichPk]
	String symElemList = ""
	for(i=0;i<3;i+=1)
		for(j=0;j<3;j+=1)
			Variable cSymElem = truncate2(normTensor3D[i][j][whichPk],2)
			if(cSymElem != 0)
				if(i==0 && j==0)
					symElemList = AddListItem("XX",symElemList)
				elseif(i==0 && j==1)
					symElemList = AddListItem("XY",symElemList)
				elseif(i==0 && j==2)
					symElemList = AddListItem("XZ",symElemList)
				elseif(i==1 && j==1)
					symElemList = AddListItem("YY",symElemList)
				elseif(i==1 && j==2)
					symElemList = AddListItem("YZ",symElemList)
				elseif(i==2 && j==2)
					symElemList = AddListItem("ZZ",symElemList)
				endif
			endif
		endfor
	endfor
	
	symElemList =  SortList(symElemList)
	//Recolor BB model pks to gray
	pkList = SortList(WaveList("pk*"+num2str(magicAngleSpec),";",""),";",16) 
	nPks = ItemsInList(pkList)
	for(i=0;i<nPks;i+=1)
		cPk = StringFromList(i,pkList,";")
		ModifyGraph/W=DFT_BB_FIT rgb($cPk)=(52428,52428,52428)
	endfor
	
	//Highlight current BB peak and show sym elements
	String pkToShow = StringFromList(whichPk,pkList)
//	DoWindow DFT_BB_FIT
	ModifyGraph/W=DFT_BB_Fit lsize($pkToShow)=2,rgb($pkToShow)=(0,0,0)
	TextBox/K/N=text0
	Tag/W=DFT_BB_Fit/C/N=text0/A=MT/L=0/X=0/Y=10 $pkToShow, pkEn,symElemList
End

Function displayfitResults(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			ControlInfo/W=ClusteringAlgorithm osThre
			Variable osThre = V_Value
			ControlInfo/W=ClusteringAlgorithm ovpThre
			Variable ovpThre = V_Value
			ControlInfo/W=ClusteringAlgorithm resAlpha
			Variable alpha = V_Value
			ControlInfo/W=ClusteringAlgorithm transSym2
			String transSym2 = S_Value
			ControlInfo/W=ClusteringAlgorithm stage3
			Variable stage3 = V_Value
			ControlInfo/W=ClusteringAlgorithm thetaList
			String thetaList = S_Value
			ControlInfo/W=ClusteringAlgorithm BBpeak
			Variable whichPk = V_Value
			ControlInfo/W=ClusteringAlgorithm mol
			String mol = S_Value
			displayDFTBBmodel(osThre,ovpThre,alpha,stage3,thetaList,whichPk,mol)
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function breakPWaveIntoClusters(pw,mName)
	Wave pw
	String mName//Model Name representative of clustering/broadening method
	
	Variable nTrans = DimSize(pw,0),i,j=0,cc=0,tc
	Make/O/N=1 clusterList = 1 
	for(i=0;i<nTrans;i+=1)
		
		tc = pw[i][16]	//current transition cluster
		
		String cweName   = "cEn_" + num2str(tc) + "_" + mName	
		String cwOSName  = "cOS_" + num2str(tc) + "_" + mName
		String cwWiName  = "cWi_" + num2str(tc) + "_" + mName
		String cwThiName = "cTh_" + num2str(tc) + "_" + mName
		Wave en = $cweName,os = $cwOSName,wi = $cwWiName,th = $cwThiName
		if((tc > numpnts(clusterList)-1) || i==0)//This is a new cluster
			Redimension/N=(tc+1) clusterList
			clusterList[tc] = 1
			Make/O/N=1 $cweName,$cwOSName,$cwWiName,$cwThiName
			Wave en = $cweName,os = $cwOSName,wi = $cwWiName,th = $cwThiName
		else
			clusterList[tc] += 1	//Update number of transitions in cluster
			Redimension/N=(clusterList[tc]) $cweName,$cwOSName,$cwWiName,$cwThiName
		endif
		
		en[clusterList[tc]-1] = pw[i][0]	//Populate energy wave
		os[clusterList[tc]-1] = pw[i][1]	
		wi[clusterList[tc]-1] = pw[i][2]	
		th[clusterList[tc]-1] = pw[i][6]				
	endfor
End

Function tdmInfo(mName,OS,OVP,alpha,molName,[corr,tr,cl])

	String mName //Model Name
	Variable OS,OVP	//Threshold values to look up clustering folder
	Variable alpha //Molecular tilt angle
	String molName	//Name of molecule
	Variable corr //Energy correction
	Variable tr //Which DFT transition to get info for
	Variable cl	//Which cluster to focus on
	
	String baseFolder = GetDataFolder(1)
	Variable p1,p2,i,j=0,k=0
	
	//Check that only the transition or the cluster is being specified. Reduce bulk in visualizer
	if(!ParamIsDefault(tr) && !ParamIsDefault(cl))
		Abort "Please only specify either a cluster or a transition to get information for."
	endif 
	
	//Add functionality to fetch the clusters made at different iterations
	//Use it to include clustered peak in tdm plot
	String mFolder = "root:Packages:DFTClustering:PolarAngles_" + molName +":"//Base folder
	String aFolder = "TransitionFiltering_"+replaceString(".",num2str(OS),"p")+"OS_"+replaceString(".",num2str(OVP),"p")+"OVP:"//Start folder of current filtering/clustering run
	
	String fFolder   = mFolder + aFolder +"AmplitudeFitting:Alpha_"+replacestring(".",num2str(alpha),"p")+":TensorMisc:"//Folder with results of DFT BB Model fit to experiment
	String iniFolder = mFolder + aFolder +"firstFilter:"
	SetDataFolder $iniFolder	
	
	Wave pw = AllParamsF1Sorted
	//Make waves containing the tdm component values within desired energy range 
	breakPWaveIntoClusters(pw,mName)
	
	//Make color wave
	Wave M_colors = makeColorWave(0.5)
	Variable nColors = DimSize(M_colors,0),row
	
	//Information about the transition of interest
	Variable atom = pw[tr][14] //What atom the transition derives from
	Variable trEn = pw[tr][0]  //Transition energy
	Variable trWi = pw[tr][2]  //Transition width
	Variable mo   = pw[tr][15] //Transition MO
	Variable cltr = pw[tr][16] //To what cluster does this transition belong to?
	SetDataFolder $fFolder	//Need to fetch the cluster parameter wave
	Wave pwI0Fit
	Wave pw2dFin = make2Dfrom1DPwave(pwI0Fit)
	Variable i0 = pwI0Fit[1]
	Wave pw2dini = pw2dOriginal
	//Information about the clusters
	Variable nClusters = DimSize(pw2dini,0)
	if(cl >= nClusters)
		return -1
	endif
	Variable clEn = pw2dini[cl][0]//Cluster Energy [This comes from the paramter wave use to fit the experiment]
	Variable clWi = truncate2(2.355*pw2dini[cl][2],2)//Cluster Width
	Variable clOS = truncate2(pw2dini[cl][1],3)//Cluster Amplitude
	Variable clTh = truncate2(pw2dini[cl][6],1)//Cluster TDM Theta***
	Variable clMT = truncate2(pwI0Fit[pwI0Fit[0] + 4 + cl*11+10],1)//Cluster modTheta
	Variable clXX = truncate2(pw2dini[cl][8],2)
	Variable clYY = truncate2(pw2dini[cl][9],2)
	Variable clZZ = truncate2(pw2dini[cl][10],2)
	//Filtered Transition List -> Trace back to transitions
	//These are variables and strings used for the cluster label
	Variable loE,hiE,loOS,hiOS,loW,hiW,maxOS 	//Range of Energies,OS and Widths.Transition w/ max OS
	String cenLabel,camLabel,cwiLabel,cthLabel,cmtLabel,tenLabel,twiLabel,tosLabel,tensorLabel,ntcLabel,trMaxLabel,clPCClabel
	
	//These are variables and strings used for the transition label
	Variable trE  = pw[tr][0]	 //Transition energy
	Variable trW  = truncate2(2.355*pw[tr][2],2)//2.355*pw[tr][2]  //Transition width
	Variable trOS = truncate2(pw[tr][1],3)  //Transition OS
	Variable trTh = truncate2(pw[tr][6],1)  //Transition TDM theta
	Variable trCl = pw[tr][16] //Cluster transition belongs to
	Variable trxx = truncate2(pw[tr][8],3)  //Transition xx tensor element
	Variable tryy = truncate2(pw[tr][9],3)  //Transition yy tensor element
	Variable trzz = truncate2(pw[tr][10],3) //Transition zz tensor element
	String trLabel
	
	SetDataFolder $iniFolder
	Wave percentDiff
	if(!WaveExists(percentDiff))
		Make/O/N=(nClusters) percentDiff
	endif
	Make/O/N=(2000) dSumPk=0,dClPk=0,dummySumPk=0,dummyClPk=0
	SetScale/I x,280,320,dSumPk,dClPk,dummySumPk,dummyClPk
	dClPk = pw2dini[cl][1] * gauss(x,clEn,pw2dini[cl][2])
	String pltName = "tdmInfoPlot_" + mName
	DoWindow $pltName
	if(!V_Flag)
		Display/W=(0,0,1200,400)/N=$pltName/K=1
		String enList = WaveList("cEn*"+mName+"",";","")
		String osList = WaveList("cOS*"+mName+"",";","")
		String wiList = WaveList("cWi*"+mName+"",";","")
		Variable nw = ItemsInList(enList)
		for(i=0;i<nw;i+=1)
			String enn = StringFromList(i,enList)
			String osn = StringFromList(i,osList)
			String win = StringFromList(i,wiList)
			Wave enw = $enn,osw = $osn,wiw = $win
			if(i==cl)
				Variable ntc = numpnts(enw)	//Number of Transitions in cluster
				loE  = enw[0];hiE  = enw[ntc-1]	//Range of energies
				WaveStats/Q osw
				loOS = truncate2(V_min,3);hiOS = truncate2(V_max,3);maxOS = V_maxloc	//Range of oscillator strengths
				WaveStats/Q wiw
				loW  = 2.355*V_min;hiW  = 2.355*V_max	//Range of peak widths
				//Determine percent difference between summed and merged peak 
				Variable transInCl = numpnts(enw)
				for(k=0;k<transInCl;k+=1)
					dSumPk += osw[k] * gauss(x,enw[k],wiw[k])//sqrt(2*Pi) * wiw[k] * osw[k] * gauss(x,enw[k],wiw[k])
				endfor
				Variable pd = calcPercentDiff(dSumPk,dClPk)//StatsCorrelation(dSumPk,dClPk)
				percentDiff[i] = pd
			else
				transInCl = numpnts(enw)
				dummySumPk = 0;dummyClPk = 0
				for(k=0;k<transInCl;k+=1)
					dummySumPk += osw[k] * gauss(x,enw[k],wiw[k])//sqrt(2*Pi) * wiw[k] * osw[k] * gauss(x,enw[k],wiw[k])
				endfor
				dummyClPk = sqrt(2*Pi) * pw2dini[i][2] * pw2dini[i][1] * gauss(x,pw2dini[i][0],pw2dini[i][2])
				Variable dummypd = calcPercentDiff(dummySumPk,dummyClPk)
				percentDiff[i] = dummypd
			endif
			AppendToGraph/W=$pltName osw vs enw
			ModifyGraph rgb($osn)=(M_colors[j][0],M_colors[j][1],M_colors[j][2],M_colors[j][3]),mode($osn)=1,offset($osn)={corr,0}
			j+=1
			if(j>=nColors)
				j = 0
			endif
		endfor 
		
		//Determine color to use for info box
		j=0
		for(i=0;i<nw;i+=1)
			if(!ParamIsDefault(cl))
				if(i==cl)
					row = j
					break
				else
					j+=1
				endif
			elseif(!ParamIsDefault(tr))
				if(i==trCl)
					row = j
					break
				else
					j+=1
				endif
			endif
			if(j>=nColors)
				j=0
			endif
		endfor
		//ModifyGraph mirror(bottom)=1,minor=1,fSize=16,lsize=1.5,fStyle=1
		Label left "ƒ[a.u.]\\U"
		Label bottom "Transition Energy[eV]"
		Variable adj = 0.7
		ModifyGraph lblPosMode(left)=1
		SetAxis bottom *,320
		SetAxis left 0,*
		//Find the wave that contains the transition to get info for
		String whichEWave  = StringFromList(trCl,enList)
		String whichOSWave = StringFromList(trCl,osList)
		Wave ew = $whichEWave
		Variable trPos = ew[BinarySearch(ew,trE)]
		
		//Fetch final clustered peaks
		SetDataFolder $fFolder
		String clusteredPeaks = WaveList("*spec2",";","")
		Variable ni = ItemsInList(clusteredPeaks)
		j=0
		//Plot clusters after fitting amplitudes to experiment
		for(i=0;i<ni;i+=1)
			String fitClPkName = "fitCl_" + num2str(i)
			Make/O/N=(2000) $fitClPkName
			Wave w = $fitClPkName
			SetScale/I x,280,360,w
			w = i0*pw2dFin[i][2]*gauss(x,pw2dFin[i][0],pw2dFin[i][1])
			AppendTograph/R w
			ModifyGraph rgb($fitClPkName)=(M_colors[j][0],M_colors[j][1],M_colors[j][2]),lsize($fitClPkName)=3//,offset($fitClPkName)={corr,0}//,lstyle($fitClPkName)=3
			j+=1
			if(j>=nColors)
				j = 0
			endif
		endfor
		j=0
		//Plot clusters BEFORE fitting amplitudes to experiment
		for(i=0;i<ni;i+=1)
			String oriClPkName = "oriCl_" + num2str(i)
			Make/O/N=(2000) $oriClPkName
			Wave w = $oriClPkName
			SetScale/I x,280,360,w
			w = i0*pw2dini[i][1]*gauss(x,pw2dini[i][0],pw2dini[i][2])
			AppendTograph/R w
			ModifyGraph rgb($oriClPkName)=(M_colors[j][0],M_colors[j][1],M_colors[j][2]),lsize($oriClPkName)=3,offset($oriClPkName)={corr,0},lstyle($oriClPkName)=3
			j+=1
			if(j>=nColors)
				j = 0
			endif
		endfor
		//Stuff for Cluster Label
		cenLabel = "\rPos="    + num2str(clEn) + "eV"
		camLabel = "  Amp=" + num2str(clOS)
		cwiLabel = "\rWid="  + num2str(clWi) + "eV"
		cthLabel = "  θ="     + num2str(clTh) + " °"
		cmtLabel = "\rθ\Bm\M= " + num2str(clMT) + " °"
		tenLabel = "Tr En Range: "  + num2str(loE)  + " - "  + num2str(hiE) + "eV"
		twiLabel = "\rTr Wid Range: " + num2str(loW)  + " - "  + num2str(hiW) + "eV"
		ntcLabel = "  #Trs: " +  num2str(ntc)
		tosLabel = "\rTr OS Range: "    + num2str(loOS) + " - " + num2str(hiOS) + "eV"
		tensorLabel = "\r\n\\JC\rxx= " + num2str(clXX) +"\ryy= "+ num2str(clYY) + "\rzz="+ num2str(clZZ)
		trMaxLabel = "\rMax OS Transition = "  + num2str(maxOS)// +  " " + num2str(maxOS)
		clPCClabel = "\r%Diff = " + num2str(pd)
		String clusterToTag = StringFromList(cl,osList)
		
		//Stuff for Transition Label
		if(!ParamIsDefault(cl))
			Tag/C/W=$pltName/N=text0/S=1/G=(M_colors[row][0],M_colors[row][1],M_colors[row][2])/D=1.25/TL={lineRGB=(M_colors[row][0],M_colors[row][1],M_colors[row][2]),lThick=2} $clusterToTag, 0,"\\K(0,0,0)\\JCCluster " + num2str(cl) + " Info" + cenLabel + camLabel + cwiLabel + cthLabel + cmtLabel + ntcLabel + trMaxLabel + clPCClabel
			AppendText/N=text0 tenLabel + twiLabel + tosLabel + "\rCluster Tensor:" + tensorLabel	
			Tag/K/W=$pltName/N=text1
		elseif(!ParamIsDefault(tr))
			Tag/C/W=$pltName/N=text1/S=1/G=(M_colors[row][0],M_colors[row][1],M_colors[row][2])/D=1.25/TL={lineRGB=(M_colors[row][0],M_colors[row][1],M_colors[row][2]),lThick=2}/TL={lThick=2} $whichOSWave, trPos,"\\JCTransition "+num2str(tr)+" Info\rEn = "+num2str(trE)+" eV W = "+num2str(trW)+" eV\rOS = "+num2str(tros)+" θ = "+num2str(trTh)+"°\rCluster: "+num2str(trCl)+"\rTransition Tensor\rxx ="+num2str(trxx)+"\ryy ="+num2str(tryy)+"\rzz ="+num2str(trzz)
			Tag/K/W=$pltName/N=text0
		endif
		Label right "Mass Absorbance [cm\S2\M/g] \U"
		ModifyGraph mirror(bottom)=1,minor=1,fSize=16,fStyle=1
		Legend/C/W=$pltName/N=text2/J/S=1/A=RT "\\s(cOS_0_"+mName+") OS\r\\s(fitCl_0) Fit BB\r\\s(oriCl_0) Ini BB"
	else
		SetDataFolder $iniFolder
		enList = WaveList("cEn*"+mName+"",";","")
		osList = WaveList("cOS*"+mName+"",";","")
		wiList = WaveList("cWi*"+mName+"",";","")
		nw = ItemsInList(enList)
		for(i=0;i<nw;i+=1)
			enn = StringFromList(i,enList)
			osn = StringFromList(i,osList)
			win = StringFromList(i,wiList)
			Wave enw = $enn,osw = $osn,wiw = $win
			if(i==cl)
				ntc = numpnts(enw)	//Number of Transitions in cluster
				loE  = enw[0];hiE  = enw[ntc-1]	//Range of energies
				WaveStats/Q osw
				loOS = truncate2(V_min,3);hiOS = truncate2(V_max,3);maxOS = V_maxloc	//Range of oscillator strengths
				WaveStats/Q wiw
				loW  = 2.355*V_min;hiW  = 2.355*V_max	//Range of peak widths
				//Determine PCC 
				transInCl = numpnts(enw)
				for(k=0;k<transInCl;k+=1)
					dSumPk += osw[k] * gauss(x,enw[k],wiw[k])//sqrt(2*Pi) * wiw[k] * osw[k] * gauss(x,enw[k],wiw[k])
				endfor
				pd = calcPercentDiff(dSumPk,dClPk)//PCC = StatsCorrelation(dSumPk,dClPk)
			else
				transInCl = numpnts(enw)
				dummySumPk = 0;dummyClPk = 0
				for(k=0;k<transInCl;k+=1)
					dummySumPk += osw[k] * gauss(x,enw[k],wiw[k])//sqrt(2*Pi) * wiw[k] * osw[k] * gauss(x,enw[k],wiw[k])
				endfor
				dummyClPk = sqrt(2*Pi) * pw2dini[i][2] * pw2dini[i][1] * gauss(x,pw2dini[i][0],pw2dini[i][2])
				dummypd = calcPercentDiff(dummySumPk,dummyClPk)
				percentDiff[i] = dummypd
			endif
		endfor
		
		//Determine color to use for info box
		j=0
		for(i=0;i<nw;i+=1)
			if(!ParamIsDefault(cl))
				if(i==cl)
					row = j
					break
				else
					j+=1
				endif
			elseif(!ParamIsDefault(tr))
				if(i==trCl)
					row = j
					break
				else
					j+=1
				endif
			endif
			if(j>=nColors)
				j=0
			endif
		endfor
		
		//Find the wave that contains the transition to get info for
		whichEWave  = StringFromList(trCl,enList)
		whichOSWave = StringFromList(trCl,osList)
		Wave ew = $whichEWave
		trPos = ew[BinarySearch(ew,trE)]
		cenLabel = "\rPos="    + num2str(clEn) + "eV"
		camLabel = "  Amp=" + num2str(clOS)
		cwiLabel = "\rWid="  + num2str(clWi) + "eV"
		cthLabel = "  θ="     + num2str(clTh) + " °"
		cmtLabel = "\rθ\Bm\M= " + num2str(clMT) + " °"
		tenLabel = "Tr En Range: "  + num2str(loE)  + " - "  + num2str(hiE) + "eV"
		twiLabel = "\rTr Wid Range: " + num2str(loW)  + " - "  + num2str(hiW) + "eV"
		ntcLabel = "  #Trs: " +  num2str(ntc)
		tosLabel = "\rTr OS Range: "    + num2str(loOS) + " - " + num2str(hiOS) + "eV"
		tensorLabel = "\r\n\\JC\rxx= " + num2str(clXX) +"\ryy= "+ num2str(clYY) + "\rzz="+ num2str(clZZ)	
		trMaxLabel = "\rMax OS Transition = "  + num2str(maxOS) 
		clPCClabel = "\r%Diff = " + num2str(pd)
		clusterToTag = StringFromList(cl,osList)
		//Stuff for Transition Label
		if(!ParamIsDefault(cl))
			Tag/C/W=$pltName/N=text0/S=1/G=(M_colors[row][0],M_colors[row][1],M_colors[row][2])/D=1.25/TL={lineRGB=(M_colors[row][0],M_colors[row][1],M_colors[row][2]),lThick=2} $clusterToTag, 0,"\\K(0,0,0)\\JCCluster " + num2str(cl) + " Info" + cenLabel + camLabel + cwiLabel + cthLabel + cmtLabel + ntcLabel + trMaxLabel + clPCClabel
			AppendText/W=$pltName/N=text0 tenLabel + twiLabel + tosLabel + "\rCluster Tensor:" + tensorLabel	
			Tag/K/W=$pltName/N=text1
		elseif(!ParamIsDefault(tr))
			Tag/C/W=$pltName/N=text1/S=1/G=(M_colors[row][0],M_colors[row][1],M_colors[row][2])/D=1.25/TL={lineRGB=(M_colors[row][0],M_colors[row][1],M_colors[row][2]),lThick=2} $whichOSWave, trPos,"\\JCTransition "+num2str(tr)+" Info\rEn = "+num2str(trE)+" eV W = "+num2str(trW)+" eV\rOS = "+num2str(tros)+" θ = "+num2str(trTh)+"°\rCluster: "+num2str(trCl)+"\rTransition Tensor\rxx ="+num2str(trxx)+"\ryy ="+num2str(tryy)+"\rzz ="+num2str(trzz)
			Tag/K/W=$pltName/N=text0
		endif 
	endif
	
	//Plot the percent difference between the summed peak and the merged peak for each cluster
	DoWindow perDiffPlot
	if(!V_Flag)
		Display/K=1/N=perDiffPlot percentDiff
		ModifyGraph/W=perDiffPlot mode=4,marker=19,grid=2,mirror=1,fStyle=1,fSize=12
		Label/W=perDiffPlot left "% Difference"
		Label/W=perDiffPlot bottom "Cluster ID"
		SetAxis/W=perDiffPlot left 0,*
	endif
	
	//Plot the summed Peak and the clustered peak along the corresponding PCC
	DoWindow sumVSCluster
	if(!V_Flag)
		Display/N=sumVSCluster/K=1 dSumPk,dClPk
		Legend/C/N=text2/J/W=sumVSCluster "\\JCCluster "+num2str(cl)+"\r\\s(dClPk) Clustered Pk\r\\s(dSumPk) Summed Pk\r%Diff = " + num2str(pd)
		ModifyGraph/W=sumVSCluster grid=2,mirror=1,nticks=10,minor=1,fStyle=1
		Label/W=sumVSCluster left "ƒ[A.U.]\U"
		Label/W=sumVSCluster bottom "Transition Enerrgy [eV]"
		ModifyGraph/W=sumVSCluster lsize=2,lstyle(dClPk)=3,rgb(dClPk)=(0,0,0)
		SetAxis/W=sumVSCluster bottom clEn-2*clWi,clEn+2*clWi
	else
		Legend/C/N=text2/J/W=sumVSCluster "\\JCCluster "+num2str(cl)+"\r\\s(dClPk) Clustered Pk\r\\s(dSumPk) Summed Pk\r%Diff = " + num2str(pd)
		SetAxis/W=sumVSCluster bottom clEn-2*clWi,clEn+2*clWi
	endif
	
	DoWindow tdmOSPlot
	SetDataFolder $iniFolder
	enn = StringFromList(cl,enList)
	osn = StringFromList(cl,osList)
	String thetaList = WaveList("cTh*"+mName+"",";","")
	String thn = StringFromList(cl,thetaList)
	Wave enw = $enn,osw = $osn,th = $thn
	WaveStats/Q enw
	if(!V_Flag)
		Display/N=tdmOSPlot/K=1 osw vs enw
		NewFreeAxis/L Theta
		AppendToGraph/W=tdmOSPlot/L=Theta th vs enw 
		ModifyGraph/W=tdmOSPlot grid=2,mirror=1,nticks=10,minor=1,fStyle=1
		Label/W=tdmOSPlot left  "ƒ[A.U.]\U"
		Label/W=tdmOSPlot Theta "θ[°]\U"
		Label/W=tdmOSPlot bottom "Transition Energy [eV]"
		ModifyGraph/W=tdmOSPlot lblPosMode(left)=1,lblPosMode(Theta)=1,axisEnab(left)={0,0.46},axisEnab(Theta)={0.5,1},freePos(Theta)=0,mirror=1,nticks=10,minor=1,fStyle=1,grid=2
		ModifyGraph/W=tdmOSPlot mode($osn)=1,mode($thn)=3,marker($thn)=19,rgb($thn)=(0,0,0),lsize($osn)=1.5
		SetAxis/W=tdmOSPlot bottom V_min,V_max
		SetAxis/W=tdmOSPlot Theta  0,90
		Legend/W=tdmOSPlot/C/N=text0/J/F=0/D=1.2/A=RC "\\JCCluster "+num2str(cl)+"\r\\s("+osn+") OS\r\\s("+thn+") θ"
	else
		String tlist = TraceNameList("tdmOSPlot",";",1)
		Variable n = ItemsInList(tlist)
		for(i=0;i<n;i+=1)
			String cTrace = StringFromList(i,tlist)
			RemoveFromGraph/W=tdmOSPlot $cTrace
		endfor
		AppendToGraph/W=tdmOSPlot osw vs enw 
		AppendToGraph/W=tdmOSPlot/L=Theta th vs enw 
		ModifyGraph/W=tdmOSPlot lblPosMode(left)=1,lblPosMode(Theta)=1,axisEnab(left)={0,0.46},axisEnab(Theta)={0.5,1},freePos(Theta)=0,mirror=1,nticks=10,minor=1,fStyle=1,grid=2
		ModifyGraph/W=tdmOSPlot mode($osn)=1,mode($thn)=3,marker($thn)=19,rgb($thn)=(0,0,0),lsize($osn)=1.5
		SetAxis/W=tdmOSPlot bottom V_min,V_max
		SetAxis/W=tdmOSPlot Theta  0,90
		Label/W=tdmOSPlot left  "ƒ[A.U.]\U"
		Label/W=tdmOSPlot Theta "θ[°]\U"
		Label/W=tdmOSPlot bottom "Transition Energy [eV]"
		Legend/W=tdmOSPlot/C/N=text0/J/F=0/D=1.2/A=RC "\\JCCluster "+num2str(cl)+"\r\\s("+osn+") OS\r\\s("+thn+") θ"	
	endif
	//Plot the Oscillators Strengths for just that one cluster and TDM Thetas
	SetDataFolder $baseFolder
End

Function/WAVE makeColorWave(alpha)
	Variable alpha
	
	Variable opacity = 65535 * alpha
	Make/O/N=(9,4) M_colors = 0 //Make color wave
	M_colors[0][0] = 230; M_colors[0][1] = 25 ; M_colors[0][2] = 75 
	M_colors[1][0] = 60 ; M_colors[1][1] = 180; M_colors[1][2] = 75 
	M_colors[2][0] = 0  ; M_colors[2][1] = 130; M_colors[2][2] = 200
	M_colors[3][0] = 145; M_colors[3][1] = 30 ; M_colors[3][2] = 180
	M_colors[4][0] = 70 ; M_colors[4][1] = 240; M_colors[4][2] = 240
	M_colors[5][0] = 250; M_colors[5][1] = 190; M_colors[5][2] = 212
	M_colors[6][0] = 0  ; M_colors[6][1] = 128; M_colors[6][2] = 128
	M_colors[7][0] = 170; M_colors[7][1] = 110; M_colors[7][2] = 40 
	M_colors[8][0] = 128; M_colors[8][1] = 0  ; M_colors[8][2] = 0  
	
	M_colors *= 256
	
	M_colors[0][3] = opacity 
	M_colors[1][3] = opacity 
	M_colors[2][3] = opacity  
	M_colors[3][3] = opacity  
	M_colors[4][3] = opacity  
	M_colors[5][3] = opacity  
	M_colors[6][3] = opacity  
	M_colors[7][3] = opacity  
	M_colors[8][3] = opacity  
	return M_colors
End

Function BIC_Panel(OS)
	Wave OS
	
	Variable nx = numpnts(OS),i
	Make/O/N=(nx)/T ost2
	Make/O/N=(nx) osScale2
	for(i=0;i<nx;i+=1)
		ost2[i] = num2str(OS[i])
		osScale2[i] = i
	endfor
	
	NewPanel/N=BICvsRCS/K=1/W=(0,0,300,250)
	Slider slider0 vert=0,side=2,limits={0,nx-1,1},pos={38,183},size={235,57},userTicks={osScale2,ost2},proc=BIC_Slider
	TitleBox SliderName title="OS[%]",pos={129,159}
	
	//Select BIC wave
	Button selBICwave                                                  ,pos={142,16} ,size={150,20}
	MakeButtonIntoWSPopupButton("BICvsRCS", "selBICwave", "notifyBICWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selBICwaveTitle                                              ,pos={36,16} ,size={115,12},title="Select BIC Wave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtoBIC

	//Select RCS wave
	Button selRCSwave                                                  ,pos={142,39} ,size={150,20}
	MakeButtonIntoWSPopupButton("BICvsRCS", "selRCSwave", "notifyRCSWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selRCSwaveTitle                                           ,pos={36,39} ,size={115,12},title="Select RCS Wave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtoRCS
			
	//Select nPeaks wave
	Button selnPeakswave                                                  ,pos={142,62} ,size={150,20}
	MakeButtonIntoWSPopupButton("BICvsRCS", "selnPeakswave", "notifynPeaksWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selnPeakwaveTitle                                             ,pos={36,62} ,size={115,12},title="Select nPeaks Wave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtonPeaks

	//Select pDiff wave
	Button selpDiffwave                                                  ,pos={142,85} ,size={150,20}
	MakeButtonIntoWSPopupButton("BICvsRCS", "selpDiffwave", "notifypDiffWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selpDiffwaveTitle                                           ,pos={36,85} ,size={115,12},title="Select pDiff Wave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtopDiff
	
	//Select OVP wave
	Button selOVPwave                                                  ,pos={142,108} ,size={150,20}
	MakeButtonIntoWSPopupButton("BICvsRCS", "selOVPwave", "notifyOVPWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selOVPwaveTitle                                           ,pos={36,108} ,size={115,12},title="Select OVP Wave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtoOVP
	
	//Select OS wave
	Button selOSwave                                                  ,pos={142,131} ,size={150,20}
	MakeButtonIntoWSPopupButton("BICvsRCS", "selOSwave", "notifyOSWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selOSwaveTitle                                           ,pos={36,131} ,size={115,12},title="Select OS Wave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtoOS

End

Function notifynPeaksWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtonPeaks = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifypDiffWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtopDiff = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyBICWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtoBIC = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyRCSWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtoRCS = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyOVPWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtoOVP = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

Function notifyOSWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtoOS = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

function/S GetnPeaksPath()
   
   SVAR pathtonPeaks = root:pathtonPeaks
   if(SVAR_Exists(pathtonPeaks))
	   	return pathtonPeaks
	else
		print "String pathtonPeaks doesn't exist. Reselect 2D nPeaks wave."
	endif
end

function/S GetpDiffPath()
   
   SVAR pathtopDiff = root:pathtopDiff
   if(SVAR_Exists(pathtopDiff))
	   	return pathtopDiff
	else
		print "String pathtopDiff doesn't exist. Reselect percent Difference wave."
	endif
end

function/S GetRCSPath()
   
   SVAR pathtoRCS = root:pathtoRCS
   if(SVAR_Exists(pathtoRCS))
	   	return pathtoRCS
	else
		print "String pathtoRCS doesn't exist. Reselect 2D reduced Chi Squared wave."
	endif
end

function/S GetBICPath()
   
   SVAR pathtoBIC = root:pathtoBIC
   if(SVAR_Exists(pathtoBIC))
	   	return pathtoBIC
	else
		print "String pathtoBIC doesn't exist. Reselect BIC wave."
	endif
end

function/S GetOVPPath()
   
   SVAR pathtoOVP = root:pathtoOVP
   if(SVAR_Exists(pathtoOVP))
	   	return pathtoOVP
	else
		print "String pathtoOVP doesn't exist. Reselect OVP wave."
	endif
end

function/S GetOSPath()
   
   SVAR pathtoOS = root:pathtoOS
   if(SVAR_Exists(pathtoOS))
	   	return pathtoOS
	else
		print "String pathtoOS doesn't exist. Reselect OS wave."
	endif
end

Function BIC_Slider(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
			endif
			Wave RCS    = $GetRCSPath()
			Wave BIC    = $GetBICPath()
			Wave nPeaks = $GetnPeaksPath()
			Wave pDiff  = $GetpDiffPath()
			Wave OVP    = $GetOVPPath()
			Wave OS     = $GetOSPath()
			ControlInfo/W=BICvsRCS slider0
			Variable row = V_Value
			dispBICandRCS(RCS,BIC,nPeaks,pDiff,row,OVP,OS)
			break
	endswitch

	return 0
End

Function dispBICandRCS(RCS,BIC,nPeaks,pDiff,row,OVP,OS)
	
	Wave RCS,BIC,nPeaks,pDiff,OVP,OS
	Variable row
	
	String name1 = NameOfWave(RCS)
	String name2 = NameOfWave(BIC)
	String name3 = NameOfWave(nPeaks)
	String name4 = NameOfWave(pDiff)
	
	if(row >= numpnts(OS))
		Abort "No further data."
	endif
	WaveStats/q RCS
	Variable maxRCS = V_max,minRCS = V_min
	WaveStats/q BIC
	Variable maxBIC = V_max,minBIC = V_min
	WaveStats/q nPeaks
	Variable maxnPeaks = V_max,minnPeaks = V_min
	WaveStats/q pDiff
	Variable maxpDiff = V_max,minpDiff = V_min
	Variable OSval = OS[row]
	DoWindow/F BICvsRCSPlot
	if(!V_Flag)
		Display/W=(0,0,400,500)/N=BICvsRCSPlot/K=1 
		AppendToGraph/W=BICvsRCSPlot RCS[row][*] VS OVP
		AppendToGraph/R/W=BICvsRCSPlot BIC[row][*] VS OVP
		ModifyGraph mode=4,marker=19,rgb($name2)=(0,0,0)
		ModifyGraph mirror(bottom)=1,minor(bottom)=1,fStyle=1
		NewFreeAxis/L nPeaksAxis
		NewFreeAxis/R pDiffAxis
		SetAxis left minRCS,maxRCS
		SetAxis right minBIC,maxBIC
		AppendToGraph/W=BICvsRCSPlot/L=nPeaksAxis nPeaks[row][*] VS OVP
		AppendToGraph/W=BICvsRCSPlot/R=pDiffAxis pDiff[row][*] VS OVP
		SetAxis nPeaksAxis 0,maxnPeaks
		SetAxis pDiffAxis 0,maxpDiff
		ModifyGraph mode=4,marker=19,rgb(nPeaks)=(1,9611,39321),rgb(pDiff)=(19729,1,39321)
		ModifyGraph axisEnab(nPeaksAxis)={0,0.49},axisEnab(pDiffAxis)={0,0.49},axisEnab(left)={0.51,1},axisEnab(right)={0.51,1},freePos(nPeaksAxis)=0,freePos(pDiffAxis)=0
		Label/W=BICvsRCSPlot right "BIC"
		Label/W=BICvsRCSPlot bottom "OVP[%]"
		Label/W=BICvsRCSPlot left "Χ\\S2\\M\\BRED\\M\\U"
		Label/W=BICvsRCSPlot nPeaksAxis "Number of Peaks"
		Label/W=BICvsRCSPlot pDiffAxis  "% Difference"
		ModifyGraph lblPosMode(nPeaksAxis)=1,lblPosMode(pDiffAxis)=1,lblPosMode(left)=1,lblPosMode(right)=1
		Legend/C/N=text0/J/S=1/A=MC "\\JCOS = "+num2str(osval)+"%\r\\JL\\s(rChiSq) Χ\\S2\\M\\BRED\\M\r\\s("+name2+") BIC\r\\s("+name3+") nPeaks\r\\s("+name4+") %Diff"
	else
		DoWindow/F BICvsRCSPlot
		String tlist = TraceNameList("BICvsRCSPlot",";",1)
		Variable n = ItemsInList(tlist),i
		for(i=0;i<n;i+=1)
			String cTrace = StringFromList(i,tlist)
			RemoveFromGraph/W=BICvsRCSPlot $cTrace
		endfor
		AppendToGraph/W=BICvsRCSPlot RCS[row][*] VS OVP
		AppendToGraph/R/W=BICvsRCSPlot BIC[row][*] VS OVP
		ModifyGraph mode=4,marker=19,rgb($name2)=(0,0,0)
		ModifyGraph mirror(bottom)=1,minor(bottom)=1,fStyle=1
		SetAxis left 0,maxRCS
		SetAxis right 0,maxBIC
		AppendToGraph/W=BICvsRCSPlot/L=nPeaksAxis nPeaks[row][*] VS OVP
		AppendToGraph/W=BICvsRCSPlot/R=pDiffAxis pDiff[row][*] VS OVP
		SetAxis nPeaksAxis 0,maxnPeaks
		SetAxis pDiffAxis 0,maxpDiff
		ModifyGraph mode=4,marker=19,rgb(nPeaks)=(1,9611,39321),rgb(pDiff)=(19729,1,39321)
		ModifyGraph axisEnab(nPeaksAxis)={0,0.49},axisEnab(pDiffAxis)={0,0.49},axisEnab(left)={0.51,1},axisEnab(right)={0.51,1},freePos(nPeaksAxis)=0,freePos(pDiffAxis)=0
		Label/W=BICvsRCSPlot right "BIC"
		Label/W=BICvsRCSPlot bottom "OVP[%]"
		Label/W=BICvsRCSPlot left "Χ\\S2\\M\\BRED\\M\\U"
		Label/W=BICvsRCSPlot nPeaksAxis "Number of Peaks"
		Label/W=BICvsRCSPlot pDiffAxis  "% Difference"
		ModifyGraph lblPosMode(nPeaksAxis)=1,lblPosMode(pDiffAxis)=1,lblPosMode(left)=1,lblPosMode(right)=1
		Legend/C/N=text0/J/S=1/A=MC "\\JCOS = "+num2str(osval)+"%\r\\JL\\s(rChiSq) Χ\\S2\\M\\BRED\\M\r\\s("+name2+") BIC\r\\s("+name3+") nPeaks\r\\s("+name4+") %Diff"
	endif
End

function truncate2(inValue,targetDP)
    Variable inValue
    Variable targetDP
    targetDP = round(targetDP)
    inValue = round(inValue * (10^targetDP)) / (10^targetDP)
    return invalue
end

Function PlottingPanel(OS,OVP)
	Wave OS,OVP
	
	Variable nx = numpnts(OS),ny=numpnts(OVP),i
	Make/O/N=(nx)/T ost2
	Make/O/N=(ny)/T ovpt2
	Make/O/N=(nx) osScale2
	Make/O/N=(ny) ovpScale2
	for(i=0;i<nx;i+=1)
		ost2[i] = num2str(OS[i])
		osScale2[i] = i
	endfor
	for(i=0;i<ny;i+=1)
		ovpt2[i] = num2str(OVP[i])
		ovpScale2[i] = i
	endfor
	
	NewPanel/N=ClusterResultPlotter/K=1/W=(0,0,850,250)
	Slider slider0 vert=0,side=2,limits={0,nx-1,1},pos={275,160},size={235,57},userTicks={osScale2,ost2},proc=DFTPlotter_Slider
	Slider slider1 vert=0,side=2,limits={0,ny-1,1},pos={15,57},size={803,57},userTicks={ovpScale2,ovpt2},proc=DFTPlotter_Slider
	TitleBox Slider1Name title="OS[%]",pos={379,130}
	TitleBox Slider2Name title="OVP[%]",pos={379,20}
End

Function DFTPlotter_Slider(sa) : SliderControl
	STRUCT WMSliderAction &sa

	switch( sa.eventCode )
		case -3: // Control received keyboard focus
		case -2: // Control lost keyboard focus
		case -1: // Control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
			endif
			Wave OVP = root:OVP
			Wave OS = root:OS
			ControlInfo/W=ClusterResultPlotter slider0
			Variable osVal = OS[V_Value]
			ControlInfo/W=ClusterResultPlotter slider1
			Variable ovpVal = OVP[V_Value]
			ControlInfo/W=ClusteringAlgorithm mol
			String molName = S_Value
			ControlInfo/W=ClusteringAlgorithm thetaList
			String thetas = S_Value
			ControlInfo/W=ClusteringAlgorithm alpha
			Variable alphaVal = V_Value
			plotRawDFTvsClusterDFT(osVal,ovpVal,molName)
			plotParamChanges(osVal,ovpVal,alphaVal,molName)
			plotDFTBBvsEXP(thetas,osVal,ovpVal,alphaVal,molName)
			break
	endswitch

	return 0
End

Function plotRawDFTvsClusterDFT(osVal,ovpVal,mol)
	Variable osVal,ovpVal
	String mol
	
	String iniFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":"//GetDataFolder(1)
	String ampFolder = iniFolder + "TransitionFiltering_" + replaceString(".",num2str(osVal),"p")+"OS_" +replaceString(".",num2str(ovpVal),"p")+"OVP:AmplitudeFitting:ALL:"
	SetDataFolder $ampFolder
	Wave Total_Specf0,fitSpecAll,mergedOriAll,res,OS,energy,En,percentDiff_All
	Variable pDiff = percentDiff_All[0]
	DoWindow RawDFTvsClusterDFT
	if(!V_Flag)
		Display/N=RawDFTvsClusterDFT/W=(0,0,400,500)/K=1 Total_Specf0,fitSpecAll,mergedOriAll vs energy
		Label/W=RawDFTvsClusterDFT left "Transition Intensity [a.u.]\\U"
		Label/W=RawDFTvsClusterDFT bottom "Transition Energy[eV]"
		NewFreeAxis/W=RawDFTvsClusterDFT/L residual
		NewFreeAxis/W=RawDFTvsClusterDFT/L oss
		AppendToGraph/W=RawDFTvsClusterDFT/L=residual res
		AppendToGraph/W=RawDFTvsClusterDFT/L=oss OS vs En
		ModifyGraph/W=RawDFTvsClusterDFT axisEnab(oss)={0,0.2},axisEnab(left)={0.25,0.75},axisEnab(residual)={0.8,1},freePos(residual)=0, freePos(oss)=0,lblPosMode=1
		ModifyGraph/W=RawDFTvsClusterDFT lsize=1.5,lstyle(fitSpecAll)=3,rgb(fitSpecAll)=(0,0,0),rgb(res)=(52428,34958,1),rgb(mergedOriAll)=(2,39321,1),grid=2,mirror=1,nticks=5,minor=1,fStyle=1,lsize(fitSpecAll)= 2.5
		ModifyGraph/W=RawDFTvsClusterDFT mode(OS)=1,rgb(OS)=(1,9611,39321),mirror(bottom)=2,log(oss)=1
		WaveStats/Q res
		SetAxis/W=RawDFTvsClusterDFT residual V_min,V_max
		Label/W=RawDFTvsClusterDFT residual "Residual\\U"
		WaveStats/Q OS
		SetAxis/W=RawDFTvsClusterDFT oss 0.001,1
		Label/W=RawDFTvsClusterDFT oss "OS [a.u.]\\U"
		Legend/C/N=text0/J/A=RC/X=4.11/Y=-1.88/W=RawDFTvsClusterDFT "OS = " + num2str(osVal) + "% OVP=" + num2str(ovpVal) + "%\r\\s(Total_Specf0) Unfiltered\r\\s(fitSpecAll) Fit\r\\s(mergedOriAll) Merged\r\\s(res) Residual\r% Diff = " + num2str(pDiff) +"%\r\\s(OS) OS"
	else
		DoWindow/F RawDFTvsClusterDFT
		String tlist = TraceNameList("RawDFTvsClusterDFT",";",1)
		Variable n = ItemsInList(tlist),i
		for(i=0;i<n;i+=1)
			String cTrace = StringFromList(i,tlist)
			RemoveFromGraph/W=RawDFTvsClusterDFT $cTrace
		endfor
		AppendToGraph/W=RawDFTvsClusterDFT Total_Specf0,fitSpecAll,mergedOriAll vs energy		
		AppendToGraph/W=RawDFTvsClusterDFT/L=residual res
		AppendToGraph/W=RawDFTvsClusterDFT/L=oss OS vs En
		ModifyGraph/W=RawDFTvsClusterDFT lsize=1.5,lstyle(fitSpecAll)=3,rgb(fitSpecAll)=(0,0,0),rgb(res)=(52428,34958,1),rgb(mergedOriAll)=(2,39321,1),grid=2,mirror=1,nticks=5,minor=1,fStyle=1,lsize(fitSpecAll)= 2.5
		ModifyGraph/W=RawDFTvsClusterDFT mode(OS)=1,rgb(OS)=(1,9611,39321),mirror(bottom)=2,log(oss)=1
		ModifyGraph/W=RawDFTvsClusterDFT axisEnab(oss)={0,0.2},axisEnab(left)={0.25,0.75},axisEnab(residual)={0.8,1},freePos(residual)=0, freePos(oss)=0,lblPosMode=1		
		WaveStats/Q res
		SetAxis/W=RawDFTvsClusterDFT residual V_min,V_max
		Label/W=RawDFTvsClusterDFT residual "Residual\\U"
		WaveStats/Q OS
		SetAxis/W=RawDFTvsClusterDFT oss 0.001,1
		Label/W=RawDFTvsClusterDFT oss "OS [a.u.]\\U"
		Label/W=RawDFTvsClusterDFT left "Transition Intensity [a.u.]\\U"
		Label/W=RawDFTvsClusterDFT bottom "Transition Energy[eV]"
		Legend/C/N=text0/J/A=RC/W=RawDFTvsClusterDFT "OS = " + num2str(osVal) + "% OVP=" + num2str(ovpVal) + "%\r\\s(Total_Specf0) Unfiltered\r\\s(fitSpecAll) Fit\r\\s(mergedOriAll) Merged\r\\s(res) Residual\r% Diff = " + num2str(pDiff) +"%\r\\s(OS) OS"
	endif
	
	SetDataFolder $iniFolder

End

Function plotParamChanges(osVal,ovpVal,alpha,mol)

	Variable osVal,ovpVal,alpha
	String mol
	
	String iniFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":"//GetDataFolder(1)
	String ampFolder = iniFolder + "TransitionFiltering_" + replaceString(".",num2str(osVal),"p")+"OS_" +replaceString(".",num2str(ovpVal),"p")+"OVP:AmplitudeFitting:Alpha_"+ replaceString(".",num2str(alpha),"p")+":TensorMisc:"
	SetDataFolder $ampFolder
	Wave ampChange,enChange,modThetaChange,fitTDMTheta,iniTDMTheta

	DoWindow ParamChangePlot
	if(!V_Flag)
		Display/N=ParamChangePlot/W=(405,0,705,500)/K=1 ampChange
		NewFreeAxis/L deltaEn 
		NewFreeAxis/L deltaTh
		NewFreeAxis/L deltaTDM
		AppendToGraph/W=ParamChangePlot/L=deltaEn enChange
		AppendToGraph/W=ParamChangePlot/L=deltaTh modThetaChange
		AppendToGraph/W=ParamChangePlot/L=deltaTDM fitTDMTheta,iniTDMTheta//tdmThetaChange
		ModifyGraph mirror=1,minor=1,fStyle=1,fSize=14,lblPosMode(left)=1,lblPosMode(deltaEn)=1,lblPosMode(deltaTh)=1,lblPosMode(deltaTDM)=1
		ModifyGraph axisEnab(left)={0,0.22},axisEnab(deltaEn)={0.25,0.47},axisEnab(deltaTh)={0.5,0.72},axisEnab(deltaTDM)={0.75,1},freePos(deltaEn)=0,freePos(deltaTh)=0,freePos(deltaTDM)=0
		ModifyGraph mode=4,marker=19,rgb(enChange)=(1,34817,52428),rgb(modThetaChange)=(1,39321,19939),rgb(iniTDMTheta)=(0,0,0)//,rgb(tdmThetaChange)=(19729,1,39321)
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
	else
		DoWindow/F ParamChangePlot
		String tlist = TraceNameList("ParamChangePlot",";",1)
		Variable n = ItemsInList(tlist),i
		for(i=0;i<n;i+=1)
			String cTrace = StringFromList(i,tlist)
			RemoveFromGraph/W=ParamChangePlot $cTrace
		endfor
		AppendToGraph/W=ParamChangePlot/L ampChange
		AppendToGraph/W=ParamChangePlot/L=deltaEn enChange
		AppendToGraph/W=ParamChangePlot/L=deltaTh modThetaChange
		AppendToGraph/W=ParamChangePlot/L=deltaTDM fitTDMTheta,iniTDMTheta//tdmThetaChange
		ModifyGraph mirror=1,minor=1,fStyle=1,fSize=14,lblPosMode(left)=1,lblPosMode(deltaEn)=1,lblPosMode(deltaTh)=1,lblPosMode(deltaTDM)=1
		ModifyGraph axisEnab(left)={0,0.22},axisEnab(deltaEn)={0.25,0.47},axisEnab(deltaTh)={0.5,0.72},axisEnab(deltaTDM)={0.75,1},freePos(deltaEn)=0,freePos(deltaTh)=0,freePos(deltaTDM)=0
		ModifyGraph mode=4,marker=19,rgb(enChange)=(1,34817,52428),rgb(modThetaChange)=(1,39321,19939),rgb(iniTDMTheta)=(0,0,0)//,rgb(tdmThetaChange)=(19729,1,39321)
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
	endif
	
	SetDataFolder $iniFolder
End

Function plotDFTBBvsEXP(thetaList,osVal,ovpVal,alpha,mol)

	String thetaList,mol
	Variable osVal,ovpVal,alpha
	
	String iniFolder = "root:Packages:DFTClustering:PolarAngles_"+mol+":"//GetDataFolder(1)
	String fitFolder = iniFolder + "TransitionFiltering_" + replaceString(".",num2str(osVal),"p")+"OS_" +replaceString(".",num2str(ovpVal),"p")+"OVP:AmplitudeFitting:Alpha_"+ replaceString(".",num2str(alpha),"p")+":TensorMisc:"
	
	Variable n = ItemsInList(thetaList),i,pDiff
	SetDataFolder $fitFolder
	Make/O/N=(n) percentDiff
	String fitNXFS = WaveList("fitResult*",";","")
	String pks = WaveList("pk*"+num2str(2),";","")
	Variable npks = ItemsInList(pks)
	String energies = WaveList("expEnergy*",";",""),enWaveName = StringFromList(0,energies)
	Wave enWave = $enWaveName 
	String expNXFS = WaveList("expSpec*",";","")	
	String steps = WaveList("dftStep*",";","")
	Wave step = $StringFromList(0,steps)
	DoWindow ModelvsExpPlot
	
	String Rlist = "0;1;39321;65535",Glist = "3204;34817;1;43690",Blist = "13107;52428;1;0"
	DoWindow ModelvsExpPlot
	if(!V_Flag)
		Display/N=ModelvsExpPlot/K=1/W=(720,0,1430,500)
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
			
			String fitName = StringFromList(i,fitNXFS)
			Wave fitW  = $fitName
			AppendToGraph/W=ModelvsExpPlot fitW vs enWave
			ModifyGraph lstyle($fitName)=3
			res =(( nxfs - fitW)/nxfs)*100	
			pDiff = calcPercentDiff(nxfs,fitW)
			ModifyGraph rgb($expName)=(r,g,b),rgb($fitName)=(r,g,b),lsize($expName)=3,lsize($fitName)=3
			
			AppendToGraph/W=ModelvsExpPlot/L=residuals res vs enWave
			ModifyGraph rgb($resName)=(r,g,b)
			percentDiff[i] = pDiff
		endfor
				
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
		SetAxis bottom 283,320
		SetAxis residuals -100,100
		SetAxis peaks *,90000
		//Make the legend
		String totalLegend = "θ[°]  EXP  DFT  %Diff\r",legendPortion = "",whichList
		whichList = fitNXFS
		for(i=0;i<n;i+=1)
			String theta = StringFromList(i,thetaList)
			String dftSpec = StringFromList(i,whichList)
			String expSpec = StringFromList(i,expNXFS)
			if(i != (n-1))
				legendPortion += theta + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) +"\r\n"
			else
				legendPortion += theta + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) 
			endif
		endfor
		Legend/C/A=RT/N=text0/J/W=ModelvsExpPlot totalLegend + legendPortion	 + "\r\\JC# of Peaks = "+num2str(npks)+"\rΧ\\S2\\M="//++"
	else
		DoWindow/F ModelvsExpPlot
		//Remove old traces
		String tlist = TraceNameList("ModelvsExpPlot",";",1)
		Variable nTraces = ItemsInList(tlist)
		for(i=0;i<nTraces;i+=1)
			String cTrace = StringFromList(i,tlist)
			RemoveFromGraph/W=ModelvsExpPlot $cTrace
		endfor
		
		for(i=0;i<n;i+=1)
			r = str2num(StringFromList(i,Rlist))
			g = str2num(StringFromList(i,Glist))
			b = str2num(StringFromList(i,Blist))		
			SetDataFolder $fitFolder
			expName = StringFromList(i,expNXFS)
			Wave nxfs = $expName
			AppendToGraph/W=ModelvsExpPlot nxfs vs enWave
			
			x = numpnts(nxfs)
			resName = "res" + num2str(i)
			Make/O/N=(x) $resName
			Wave res = $resName
			
			fitName = StringFromList(i,fitNXFS)
			Wave fitW  = $fitName
			res =(( nxfs - fitW)/nxfs)*100	
			pDiff = calcPercentDiff(nxfs,fitW)
			percentDiff[i] = pDiff
			AppendToGraph/W=ModelvsExpPlot fitW vs enWave
			ModifyGraph lstyle($fitName)=3
			res =(( nxfs - fitW)/nxfs)*100	
			pDiff = calcPercentDiff(nxfs,fitW)
			ModifyGraph rgb($expName)=(r,g,b),rgb($fitName)=(r,g,b),lsize($expName)=3,lsize($fitName)=3
			
			AppendToGraph/W=ModelvsExpPlot/L=residuals res vs enWave
			ModifyGraph rgb($resName)=(r,g,b)
			percentDiff[i] = pDiff
		endfor
		
		//Append the DFT step edge onto graph
		AppendToGraph/W=ModelvsExpPlot step vs enWave
		ModifyGraph lsize($StringFromList(0,steps))=2,rgb($StringFromList(0,steps))=(0,0,0)
		
		for(i=0;i<npks;i+=1) 	
			pkName = StringFromList(i,pks)
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
		SetAxis bottom 283,320
		SetAxis residuals -100,100
		SetAxis peaks *,90000
		
		totalLegend = "θ[°]  EXP  DFT  %Diff\r"
		legendPortion = ""
		whichList = fitNXFS
		for(i=0;i<n;i+=1)
			theta = StringFromList(i,thetaList)
			dftSpec = StringFromList(i,whichList)
			expSpec = StringFromList(i,expNXFS)
			if(i != (n-1))
				legendPortion += theta + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) +"\r\n"
			else
				legendPortion += theta + "° \\s('"+ expSpec +"') \\s('" + dftSpec +"') " + num2str(percentDiff[i]) 
			endif
		endfor
		Legend/C/A=RT/N=text0/J/W=ModelvsExpPlot totalLegend + legendPortion	 + "\r\\JC# of Peaks = "+num2str(npks)+"\rΧ\\S2\\M="//++"
	endif
	SetDataFolder $iniFolder
End

Function overlapVisualizerPanel()
	NewPanel/N=overlapVisualizer/K=1/W=(0,0,400,100)
	SetVariable pk1 vert=0,side=2,limits={0,10000,1},pos={60,20},size={100,18},value=_NUM:0,proc=ovpVisualizerSV
	SetVariable pk2 vert=0,side=2,limits={0,10000,1},pos={60,50},size={100,18},value=_NUM:0,proc=ovpVisualizerSV
	TitleBox pk1Name title="Pk1",pos={20,20}
	TitleBox pk2Name title="Pk2",pos={20,50}
	
	//Select pWave
	Button selpWave                                                  ,pos={215,50} ,size={150,20}
	MakeButtonIntoWSPopupButton("overlapVisualizer", "selpWave", "notifypWaveSelection", options=PopupWS_OptionFloat)
	TitleBox selpWaveTitle                                              ,pos={256,20} ,size={115,12},title="Select pWave",frame=0//,fColor=(fc1,fc2,fc3)
	String/G pathtopWave
End

Function notifypWaveSelection(event, wavepath, windowName, ctrlName)
	Variable event
	String wavepath
	String windowName
	String ctrlName
	String startFolder = GetDataFolder(1)
	SetDataFolder root:
	String/G  pathtopWave = wavePath
	SetDataFolder $startFolder
	print "Selected wave:",wavepath, " using control", ctrlName
end

function/S GetpWavePath()
   
   SVAR pathtopWave = root:pathtopWave
   if(SVAR_Exists(pathtopWave))
	   	return pathtopWave
	else
		print "String pathtopWave doesn't exist. Reselect pWave."
	endif
end

Function ovpVisualizerSV(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
						
			ControlInfo/W=overlapVisualizer pk1
			Variable pk1 = V_Value
			ControlInfo/W=overlapVisualizer pk2
			Variable pk2 = V_Value
			
			Wave pWave = $GetpWavePath()
			Variable mu1  = pWave[pk1][0]  , mu2  = pWave[pk2][0]
			Variable sd1  = pWave[pk1][2]  , sd2  = pWave[pk2][2]
			Variable amp1 = pWave[pk1][1]  , amp2 = pWave[pk2][1]
			compGaussErf2(mu1,sd1,amp1,mu2,sd2,amp2,2000,0)//,d=1,pk1=pk1,pk2=pk2)
			DoUpdate
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End