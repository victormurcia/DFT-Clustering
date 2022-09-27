#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Function makeBAWrap()
	DFT_Init()
	DFT_atomsPanel()
End

Function DFT_Init()

	print "------------------------------------------------- Analysis of NEXAFS via DFT calculations ---------------------------------------------------"
	print "Victor Manuel Murcia Ruiz ----------------------------- Washington State University-------------------------------------victor.murcia@wsu.edu"
	
	//Setup Folders
	String CurrentFolder=GetDataFolder(1)
	print "Initiallizing NEXAFS processing..."
	NewDataFolder/O/S root:Packages
	NewDataFolder/O/S root:Packages:NXA
	NewDataFolder/O/S root:Packages:NXA:ElementLibrary
	
	//Add new elements here (also, change create molecule panel and mol structure)
	NewPath/Z/Q/O ElementPath SpecialDirPath("Igor Pro User files",0,0,0)+"User Procedures:Element Library:"
	LoadWave/a/J/o/b="N=Elements,F=-2;N=AtomicWeight,F=-1;"/p=ElementPath "AtomicWeight.txt"
	WAVE/T Elements
	WAVE AtomicWeight
	SetDataFolder root:Packages:NXA
	
	String/G molName="", wList //list of molecules for the popupMenu
	Make/o/n=(1,2) ChemFormSelWave
	Make/o/t/n=(1,2) ChemFormTxtWave
	
	//Load Elemental Data data
	Variable/G Errors=0
	PathInfo ElementPath
	If( V_Flag==0 || !WaveExists(:ElementLibrary:H_f1) )
		DFT_LoadBareAtomData()//Load the Element Scattering Factor Data
	endif
	
	SetDataFolder $CurrentFolder
	If( Errors==0 )
		Print "DFT Analysis code successfully initialized."
	elseIf( Errors==5 )
		Print "User Cancel: DFT Analysis code NOT initialized."
	else
		Print "ERROR: DFT Analysis code NOT initiailzed properly."
	endif
End

Function DFT_LoadBareAtomData()
	WAVE/T Elements=:ElementLibrary:Elements
	WAVE AtomicWeight=:ElementLibrary:AtomicWeight
	Variable i, numE=numpnts(Elements)
	NewPath/Z/Q/O ElementPath SpecialDirPath("Igor Pro User files",0,0,0)+"User Procedures:ElementLibrary"
	String fileList=IndexedFile( ElementPath, -1, ".nff") //see if files are actually in here
	For( ; itemsInLIst(fileList) < numE; ) //Loop, letting user select another folder
		NewPath/M="Can't find all the data here.  Select another folder."/o/q ElementPath
		If( V_Flag==-1 ) //user pressed cancel
			KillDataFolder root:Packages:NXA
			KillPath/Z ElementPath
			return -1
		endif
		fileList=IndexedFile( ElementPath, -1, ".nff")
	endfor
	//Start with Carbon element as template
	SetDataFolder root:Packages:NXA:ElementLibrary
	LoadWave/a/d/g/o/b="N='C_energy';N=C_f1raw;N=C_f2raw;"/p=ElementPath "c.nff"
	Wave C_energy,C_f1raw,C_f2raw
	Duplicate/d/o C_energy, ::mu_energy //this goes in the parent NXA folder
	Note/K/NOCR C_f1raw, "W:"+num2str(AtomicWeight[12])+";"
	Note/K/NOCR C_f2raw, "W:"+num2str(AtomicWeight[12])+";"
	For( i=0; i<numE; i+=1 )
		If( stringmatch(Elements[i],"C") )
			continue //Already did Carbon
		else
			print "Loading "+Elements[i]
			DFT_LoadHenke(i) //Loads rest of files and adds extra energies from each element into Mu_energy
		endif
	Endfor
	Sort ::mu_energy, ::mu_energy
	For( i=0; i<numE; i+=1 )
		DFT_interpHenke(Elements[i]) // interpolates all element scattering factors based on Mu_energy
	Endfor
	SetDataFolder root:Packages:NXA
End

//loads the Henke data
Function DFT_LoadHenke(index)
	Variable index
	WAVE/T Elements
	WAVE AtomicWeight
	String element=Elements[index]
	String loadTxt, fileName
	Sprintf loadTxt, "N=%s_energy;N=%s_f1raw;N=%s_f2raw;", element, element, element
	fileName= element+".nff"
	LoadWave/a/d/g/o/b=loadTxt/p=ElementPath fileName
	WAVE/Z energy=$(element+"_energy"), f1raw=$(element+"_f1raw"), f2raw=$(element+"_f2raw")
	WAVE/Z mu_energy=::mu_energy
	Note/K/NOCR f1raw, "W:"+num2str(AtomicWeight[index])+";"
	Note/K/NOCR f2raw, "W:"+num2str(AtomicWeight[index])+";"
	Duplicate/d/o f2raw, derivative, logEnergy
	Duplicate/d/o energy energyCopy
	logEnergy=log(energy)
	Differentiate derivative /X=logEnergy
	Variable i, npts
	For( i=0;  ; i+=1 ) //add extra energies at the edges into the general energy wave
		WaveStats/Q/Z derivative
		If( V_max>500 )
			Insertpoints inf, 1, mu_energy
			mu_energy[inf]=energyCopy[V_maxLoc]
			DeletePoints V_maxLoc, 1, derivative, energyCopy //don't find this point again
		else
			break
		endif
	Endfor
End

//interpolates the Heneke Data with the higher data density
Function DFT_InterpHenke(element)
	String element
	Duplicate/d/o ::mu_energy $(element+"_f2")
	WAVE/Z energy=::mu_energy, f2=$(element+"_f2"), f2raw=$(element+"_f2raw"), ElementE=$(element+"_energy")
	f2=interp(energy,ElementE,f2raw)
End

//Struct for constructing a molecule bare atom absorption spectrum
Structure DFT_molStruct
	SVAR molName
	WAVE mu, energy, f1, f2, SelWave, A, W
	WAVE/T TxtWave, Elements
EndStructure

Function DFT_fillMolStruct(m)
	Struct DFT_molStruct &m
	String CurrentFolder=GetDataFolder(1)
	SetDataFolder root:Packages:NXA
	SVAR/Z m.molName
	WAVE/Z m.f1, m.f2, m.SelWave=ChemFormSelWave
	WAVE/Z/T m.TxtWave=ChemFormTxtWave
	WAVE/Z m.A=AtomsWave
	WAVE/Z m.mu, m.energy=mu_energy
	SetDataFolder root:Packages:NXA:ElementLibrary
	WAVE/Z/T m.Elements
	WAVE/Z m.W=AtomicWeight
	SetDataFolder $currentFolder
End

//Makes a panel for constructing a molecule from constituent atoms to tabulate the mass abs coef
Function DFT_atomsPanel() : Panel
	
	ControlInfo/W=ClusteringAlgorithm mol
	String molName = S_Value
	
	NewPanel /W=(557,66,940,267)/N=atomsPanel as "Chemical Fomula"
	Button H ,pos={12,41}  ,size={25,20},proc=DFT_Button,title="H",fSize=12
	Button He,pos={210,41} ,size={25,20},proc=DFT_Button,title="He",fSize=12
	Button Li,pos={12,62}  ,size={25,20},proc=DFT_Button,title="Li",fSize=12
	Button Be,pos={38,62}  ,size={25,20},proc=DFT_Button,title="Be",fSize=12
	Button B ,pos={80,62}  ,size={25,20},proc=DFT_Button,title="B",fSize=12
	Button C ,pos={106,62} ,size={25,20},proc=DFT_Button,title="C",fSize=12
	Button N ,pos={132,62} ,size={25,20},proc=DFT_Button,title="N",fSize=12
	Button O ,pos={158,62} ,size={25,20},proc=DFT_Button,title="O",fSize=12
	Button F ,pos={184,62} ,size={25,20},proc=DFT_Button,title="F",fSize=12
	Button Ne,pos={210,62} ,size={25,20},proc=DFT_Button,title="Ne",fSize=12
	Button Na,pos={12,83}  ,size={25,20},proc=DFT_Button,title="Na",fSize=12
	Button Mg,pos={38,83}  ,size={25,20},proc=DFT_Button,title="Mg",fSize=12
	Button Al,pos={80,83}  ,size={25,20},proc=DFT_Button,title="Al",fSize=12
	Button Si,pos={106,83} ,size={25,20},proc=DFT_Button,title="Si",fSize=12
	Button P ,pos={132,83} ,size={25,20},proc=DFT_Button,title="P",fSize=12
	Button S ,pos={158,83} ,size={25,20},proc=DFT_Button,title="S",fSize=12
	Button Cl,pos={184,83} ,size={25,20},proc=DFT_Button,title="Cl",fSize=12
	Button Ar,pos={210,83} ,size={25,20},proc=DFT_Button,title="Ar",fSize=12
	Button K ,pos={12,104} ,size={25,20},proc=DFT_Button,title="K",fSize=12
	Button Ca,pos={38,104} ,size={25,20},proc=DFT_Button,title="Ca",fSize=12
	Button Sc,pos={12,135} ,size={25,20},proc=DFT_Button,title="Sc",fSize=12
	Button Ti,pos={38,135} ,size={25,20},proc=DFT_Button,title="Ti",fSize=12
	Button V ,pos={64,135} ,size={25,20},proc=DFT_Button,title="V",fSize=12
	Button Cr,pos={90,135} ,size={25,20},proc=DFT_Button,title="Cr",fSize=12
	Button Mn,pos={116,135},size={25,20},proc=DFT_Button,title="Mn",fSize=12
	Button Fe,pos={142,135},size={25,20},proc=DFT_Button,title="Fe",fSize=12
	Button Co,pos={168,135},size={25,20},proc=DFT_Button,title="Co",fSize=12
	Button Ni,pos={194,135},size={25,20},proc=DFT_Button,title="Ni",fSize=12
	Button Cu,pos={220,135},size={25,20},proc=DFT_Button,title="Cu",fSize=12
	Button Zn,pos={246,135},size={25,20},proc=DFT_Button,title="Zn",fSize=12
	Button Ga,pos={80,104} ,size={25,20},proc=DFT_Button,title="Ga",fSize=12
	Button Ge,pos={106,104},size={25,20},proc=DFT_Button,title="Ge",fSize=12
	Button As,pos={132,104},size={25,20},proc=DFT_Button,title="As",fSize=12
	Button Se,pos={158,104},size={25,20},proc=DFT_Button,title="Se",fSize=12
	Button Br,pos={184,104},size={25,20},proc=DFT_Button,title="Br",fSize=12
	Button Kr,pos={210,104},size={25,20},proc=DFT_Button,title="Kr",fSize=12
	ListBox ElementList,pos={280,40},size={85,116},proc=DFT_ElementListBox
	ListBox ElementList,listWave=root:Packages:NXA:ChemFormTxtWave,selWave=root:Packages:NXA:ChemFormSelWave,mode= 5
	SetVariable Name,pos={16,12},size={112,16},bodyWidth=80,title="Name"
	SetVariable Name,value=_STR:molName//value= _STR:""
	SetVariable Formula,pos={156,13},size={191,16},bodyWidth=150,title="Formula"
	SetVariable Formula,value= _STR:"",noedit= 1
	Button AtomsOkay,pos={145,170},size={50,20},proc=DFT_Button,title="Okay"
	Button AtomsOkay,fSize=12,fColor=(0,52224,0)
	Button AtomsCancel,pos={202,170},size={50,20},proc=DFT_Button,title="Cancel"
	Button AtomsCancel,fSize=12,fColor=(65280,0,0)
	TitleBox ElementTitle,pos={95,38},size={43,13},title="Elements",frame=0
End

Function DFT_Button(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	If( ba.eventCode!=2 )
		return 0
	endif
	String dn, rFolder, dFolder
	String CurrentFolder=GetDataFolder(1)
	STRUCT DFT_molStruct m; DFT_fillMolStruct(m)
	STRUCT DFT_dStruct d

	variable Num
	If( stringmatch(ba.ctrlName,"newMol") )
		redimension/n=(0,2) m.SelWave, m.TxtWave
		DFT_atomsPanel()
	elseif( stringMatch(ba.win,"atomsPanel") )
		If( cmpstr( ba.ctrlname, "atomsOkay")==0 )
			DFT_calcMu(m)
			DFT_fill_dStruct(d)
			DoWindow/k atomsPanel
		elseIf( cmpstr( ba.ctrlname, "atomsCancel")==0 )
			DoWindow/k atomsPanel		
		else
			DFT_AddElement(ba.ctrlname)
		endif
	
	endif
End

//calculate per-mol absorption for an isolated carbon atom
Function DFT_calcMu(m)//, [molN])
	Struct DFT_molStruct &m
//	String molN
	String newName //new molecule name without dialog box
	String currentFolder=getDataFolder(1), titleStr
	SetDataFolder root:Packages:NXA
	WAVE m.energy = Mu_energy
	Wave Mu_Energy
	Duplicate/d/o Mu_Energy, lambda, m.mu, fSum
	fSum=0
	variable wSum=0, zStar=0
	Variable Na=6.0221415e23 //Avogadro's number [at/mol]
	Variable re=2.81794e-13 //classical electron radius [cm]
	lambda=1.23984e-4/m.energy //wavelength [cm]
	Controlinfo/W=atomsPanel Name
	If( StrLen(S_Value) )  //get molecular information from the atoms panel
		Duplicate/d/o m.W $(S_Value+"_Atoms")
		WAVE m.A=$(S_Value+"_Atoms")
		m.molName=S_Value
//		molN=m.molName
		DFT_ReadMolList(m)
//	If( !ParamIsDefault(molN) ) // Not being called from the Atoms Panel
//		Duplicate/d/o m.W $(molN+"_Atoms")
	else // called from another function (Atom list & molName must be made already)
		WAVE m.A=$(m.molName+"_Atoms")
	endif
	variable i, j, z
	For( i=0; i<numpnts(m.Elements); i+=1 )
		WAVE f2=$(":ElementLibrary:"+m.Elements[i]+"_f2") //imaginary scattering factor of each atom
		If( !waveExists(f2) )
			print "Error: Couldn't find "+m.Elements[i]+"_f2 wave!!!"
			continue
		endif
		fSum+=m.A[i]*f2
		wSum+=m.A[i]*m.W[i] //molecular weight
		z=i+1 //atomic number
		For( j=0; j<m.A[i]; j+=1 )
			zStar+=z-((z)/82.5)^2.37 //relativistic correction to Z (See X-ray Dat Booklet by Thompson et al.)
		endfor
	EndFor
	
	m.mu=2*re*lambda*Na*fSum/wSum  //Actual calculation of mu [cm^2/g]
	
	Note/NOCR m.mu, ("Mol:"+m.molName+";Mw:"+num2str(wSum)+";Zstar:"+num2str(zStar)+";")
	Duplicate/d/o m.mu $(m.MolName+"_mu")
	DoWindow/F MuWindow //display results
	If( V_Flag==0 && !stringmatch(m.molName,"Carbon") )
		display/k=1/W=(252.75,272.75,542.25,461) m.mu vs m.energy
		DoWindow/C MuWindow
		ModifyGraph margin(left)=36,margin(bottom)=29,margin(top)=14,margin(right)=14
		ModifyGraph log=1,grid=2,tick=2,mirror=1
		Label Bottom "Energy [eV]"
		Label Left "Mass Abs [cm\\S2\\M/g]"
		titleStr=m.MolName+" Bare Atom Absorption"
		TextBox/C/N=text0/H={0,0.5,10}/A=MC/X=-13.68/Y=-28.57 titleStr
	elseif(!stringmatch(m.molName,"Carbon" ) )
		titleStr=m.MolName+" Bare Atom Absorption"
		TextBox/W=MuWindow/C/N=text0/H={0,0.5,10}/A=MC/X=-13.68/Y=-28.57 titleStr
	endif
	SetDataFolder $CurrentFolder
End

Function DFT_Fill_dStruct(d,[dn]) //dn is the base name of the data set (stripped of "_OD" for example)
	Struct DFT_dStruct &d
	String dn //use current folder and given wave name, skip finding it in the raw data panel
	String dFolder
	String CurrentFolder=GetDataFolder(1)
	DoWindow NEXAFSrawData
	If( ParamIsDefault(dn) && V_Flag==1 ) //if no data name given, look it up in the current graph (Raw data panel)
		dn=StringFromList(0,ImageNameList("NEXAFSrawData",";")) //name of the data
		dFolder=RemoveEnding(StringByKey( "ZWAVEDF",ImageInfo("NEXAFSrawData",dn,0) ),"Raw:")
		If( strLen(dn)==0 || strLen(dFolder)==0 )
			return 0
		endif
	elseif( ParamIsDefault(dn) ) //find the data info from the top graph
		String wn0=StringFromList(0,tracenamelist("",";",1))
		dn=removeEnding(ParseFilePath(1,wn0,"_",1,0),"_") //remove extra tags like "_ma" or "_OD" or "_delta" to get the base name
		WAVE wv=TraceNameToWaveRef("",wn0)
		String wvFolder=GetWavesDataFolder(wv,1)
		dFolder=removeEnding(removeEnding(removeEnding(wvFolder,"Raw:"),"MassAbs:"),"KramersKronig:") //parent folder for the data tree
	else
		dFolder=removeEnding(removeEnding(removeEnding(CurrentFolder,"Raw:"),"MassAbs:"),"KramersKronig:") //parent folder for the data tree
	endif
	SetDataFolder root:Packages:NXA
	SVAR/Z d.mol=MolName
	SVAR/Z d.SLmol, d.SCmol
	WAVE/Z d.muE=mu_energy
	
	SetDataFolder $CurrentFolder
End

//Data Structure
Structure DFT_dStruct
	SVAR mol
	SVAR SLmol, SCmol
	Variable Dwell
	WAVE ma, maE, maU
	WAVE muE 
	//not filled by "Fill_dStruct"
	WAVE OD, uOD, mu, muU
	Variable bkg, corrDetSat
EndStructure

Function DFT_AddElement(element)
	String Element
	Struct DFT_molStruct m; DFT_fillMolStruct(m)
	Variable nElements=DimSize(m.TxtWave,0), Enum=DFT_WhichWaveItem(element,m.Elements), i, InList=0, listEnum
	m.SelWave[0,nElements-1][1]=2 //deselect all
	For( i=0; i<nElements; i+=1 ) //search for that element in the current elemental list
		If( stringMatch(m.TxtWave[i][0], element) )
			m.TxtWave[i][1]= num2str(str2num(m.TxtWave[i][1])+1) //increment that element
			m.SelWave[i][1]= 3 //highlight the cell
			DFT_MakeFormula()
			return 0
		endif
	endFor
	//not in the current list, so add it
	For( i=0; i<nElements; i+=1 ) //find position to insert the new element
		listEnum=DFT_WhichWaveItem(m.TxtWave[i][0],m.Elements)
		If( listEnum>Enum )
			break
		endif
	Endfor
	InsertPoints i, 1, m.txtWave, m.SelWave
	m.TxtWave[i][0]=Element
	m.TxtWave[i][1]="1"
	m.SelWave[i][0]=0
	m.SelWave[i][1]=3
	DFT_MakeFormula()
end

//Returns the index of a text wave that points to the match to the given string (returns -1 if no match)
Function DFT_WhichWaveItem(str,txtW)
	String str
	WAVE/T txtW
	variable i, index=-1
	For( i=0; i<numpnts(txtW); i+=1 )
		If( stringMatch(str, txtW[i]) )
			index=i
			break
		endif
	endfor
	return index
end

//Display the chemical formula for the user to see in the atoms panel
Function DFT_MakeFormula()
	WAVE/T List=root:Packages:NXA:ChemFormTxtWave
	String str=""
	Variable i, nElements=DimSize(List,0)
	For( i=0; i<nElements; i+=1 ) //make the string from the list
		str+=List[i][0]
		If( str2num(List[i][1])<1 || str2num(List[i][1])>1 )
			str+=List[i][1]
		Endif
		If( i<nElements-1 )
			str+=","
		endif
	endfor
	SetVariable Formula value=_STR:str, win=atomsPanel
End

Function DFT_ElementListBox(ctrlName,row,col,event) : ListBoxControl
	String ctrlName
	Variable row
	Variable col
	Variable event	//1=mouse down, 2=up, 3=dbl click, 4=cell select with mouse or keys
					//5=cell select with shift key, 6=begin edit, 7=end
	Struct DFT_molStruct m; DFT_fillMolStruct(m)
	If( stringmatch(ctrlName,"ElementList") )
		If( event==7)
			IF( stringMatch(m.txtWave[row][col],"0") ) //user wants to get rid of element
				DeletePoints row, 1, m.txtWave, m.selwave
			endif
			DFT_MakeFormula()
		endif
	endif
	return 0
End

//Reads the Listbox Listwave in the atoms panel to make the elements list "m.A"
Function DFT_ReadMolList(m)
	Struct DFT_molStruct &m
	Variable i, j
	m.A=0
	For( i=0; i<numpnts(m.Elements); i+=1 )
		For( j=0; j<DimSize(m.TxtWave,0); j+=1 )
			If( StringMatch(m.Elements[i],m.TxtWave[j][0]) )
				m.A[i]=str2num(m.TxtWave[j][1])
				break
			endif
		Endfor
	Endfor
End
