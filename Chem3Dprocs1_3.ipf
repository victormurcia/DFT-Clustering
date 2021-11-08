//#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma version = 1.3
#pragma IgorVersion = 7.0.3
#pragma independentmodule=Chem3Dmodule 

// CHEM3D
// A molecular visualization package for Igor 7
//
// Copyright (c) 2017 Richard Knochenmuss
// See end of this file for license 

constant bondtol=1.15 // note tolerance factor for making bonds, based on covalent radii
constant defaultatomscale=1.0 // initial atom sphere sizes   1=vdW radii
constant minatom=0.05 // Atom size slider min value
constant maxatom=1.0 // Atom size slider max value,  1=vdW radii
constant bondthick=0.06 // base bond thickness
constant bondthickmax=8 // bond thickness slider max value

constant GCdefaultpanelmode=1 // <0:minimized, 1:object, 2:axes, 3:ranges, 4:rotation
								 	// valid minimized values are -1, -2, -3, -4
constant BgR=56797  // panel background R,G,B
constant BgG=56797
constant BgB=56797

constant macpaneldefault=1.5  // scale factors for Gizmo Control panels
constant winpaneldefault=1

constant defX=100	// initial gizmo window position and size
constant defY=100
constant defdXwin=300
constant defdYwin=300
constant defdXmac=500
constant defdYmac=500


//---------------------------------------------------------------//
Menu "Chem 3D"
	"Make new Chem3D/F4",/Q, Chem3Dmodule#Chem3Dstart()
	"Chem 3D Notes",/Q, Chem3Dmodule#Chem3Dnotes()
	submenu  "Advanced Settings"
		"Molecule scatter plot",/Q,Chem3Dmodule#Chem3Dadvanced(1) 
		"Axes",/Q,Chem3Dmodule#Chem3Dadvanced(2) 
		"Gizmo Info Window",/Q,Chem3Dmodule#Chem3Dadvanced(3) 
	end
	"Igor 3D help",/Q,Displayhelptopic "3D Graphics"
end		
								 	
//---------------------------------------------------------------//
Function Chem3Dstart()
	string wlist, str
	variable ii, defPres
	
	defPres=panelresolution("") // save existing panel setting
	str="SetIgorOption PanelResolution = "+num2str(screenresolution) // use full resln
	execute str
	
	Newpanel/K=1/W=(97,89,347,291) 
	ModifyPanel cbRGB=(BGR, BGG, BGB)
	Dowindow/C Chem3DstartPanel
	
	SetDrawLayer UserBack
	SetDrawEnv fsize= 14,fstyle= 1
	DrawText 58,26,"Chem3D Start Panel"
	DrawText 21,53,"Select wave with coordinates (Angstrom)"
	DrawText 50,67,"(2D, N atoms x 3 columns)"
	PopupMenu popxyz,pos={36.00,69.00},size={168.00,23.00},title="XYZ wave "
	PopupMenu popxyz,font="Tahoma",fSize=12
	PopupMenu popxyz,mode=1,value=wavelist("*",";","MINCOLS:3,MAXCOLS:3")
	
	DrawText 47,123,"Select text wave with atoms"
	PopupMenu popatoms,pos={39.00,127.00},size={168.00,23.00},title="Atom list "
	PopupMenu popatoms,font="Tahoma",fSize=12

	PopupMenu popatoms,mode=1,value=wavelist("*",";","TEXT:1")
	
	Button Make3dbutton,pos={99.00,159.00},size={50.00,20.00},title="Go"
	Button Make3dbutton,proc=Chem3DgoButtonProc
	
	str="SetIgorOption PanelResolution = "+num2str(defPres) // restore panel resln
	execute str
end	

//---------------------------------------------------------------//
Function Chem3DgoButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			controlinfo/W=Chem3DstartPanel popxyz
			wave xyz=$s_value
			controlinfo/W=Chem3DstartPanel popatoms
			wave atoms=$s_value
			Chem3D(atoms, xyz, 1) 
			KillWindow Chem3DstartPanel
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End			

//---------------------------------------------------------------//
Function Chem3Dadvanced(mode)  //1=scatter (molecule), 2=axes, 3 info window
	variable mode 
	string wstr, gizname
	variable ii, jj
	wstr=winrecreation("",0)
	ii=strsearch(wstr, "Chem3D_",0)
	if (ii<0)
		DoAlert/T="Chem3D", 0,"Top window is not a Chem3D Gizmo"
	else
		jj=strsearch(wstr,"_GCpanel",0,2)
		if (jj>=0)
			gizname=wstr[ii, jj-1]
			//print gizname
			dowindow/F $gizname
		else
			jj=strsearch(wstr,"(",0,2)
			gizname=wstr[ii, jj-1]
			//print gizname
		endif
		
		if (mode==1)
			modifygizmo/N=$gizname edit={object,scatter0}
		elseif(mode==2)
			modifygizmo/N=$gizname edit={object,axes0}
		elseif(mode==3)
			modifygizmo/N=$gizname showinfo
		endif
	endif
End				 	

//---------------------------------------------------------------//
Function Chem3D(atoms, xyz, bonds01) // element list, triplet with x,y,z in Angstrom
	wave/T atoms
	wave xyz
	variable bonds01
	variable npts, ii, jj, r1, r2, dd, dispnr

	string colwv, sizewv, winnam, bondnam, str, str2
	String savedDataFolder = GetDataFolder(1)

	npts=dimsize(xyz,0)
	
	for (ii=0; ii<50; ii+=1)
		winnam="Chem3D_"+num2str(ii)
		dowindow $winnam
		//print winnam, v_flag
		if (v_flag==0)
			break
		endif
	endfor

	newdatafolder/O/S $winnam
	
	// make color wave
	colwv= "AtomColors"//nameofwave(xyz)+"_C"
	make/O/N=(npts,4) $colwv  // RGBA
	wave cwv=$colwv
	
	// Make size wave
	sizewv= "AtomSizes" //nameofwave(xyz)+"_S"
	make/O/N=(npts,3) $sizewv
	variable/G AtomScale=1
	variable/G AtomMode=0 // 0=atoms, 1=at nrs
	wave siz=$sizewv
	for (ii=0; ii<npts; ii+=1)  
		r1=2*Chem3DvdwRadius(atoms[ii])
		siz[ii][]=r1*defaultatomscale*AtomScale // note scale factors
	endfor
	
	// atom numbers
	Duplicate/O/T atoms, AtomNumbers
	AtomNumbers=num2str(p)
	
	// set atom colors
	Chem3DatomColors(atoms, cwv)
	
	str2=igorinfo(2) // OS win or mac?
	variable/G GCpscl
	if (cmpstr(str2,"Windows")==0)
		GCpscl=winpaneldefault
	else
		GCpscl=macpaneldefault
	endif
	
	// Gizmo control variables
	variable GCx, GCy, GCdx, GCdy
	GCx=defX //100
	GCy=defY //100
	if (cmpstr(str2,"Windows")==0)
		GCdx=defdXwin //300
		GCdy=defdYwin //300
	else	
		GCdx=defdXmac //500
		GCdy=defdXmac //500
	endif
	variable/G  GCzoom, GCaxes, GCpanX, GCpanY, GCrefresh
	GCzoom=1
	GCaxes=1
	NewGizmo/K=1/W=(GCx, GCy, GCx+GCdx, GCy+GCdy)
	DoWindow/C $winnam
	ModifyGizmo scalingOption=63
	
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo modifyObject=light0,objectType=light,property={position,-0.2515,-0.3591,0.8988,0.0000}
	ModifyGizmo modifyObject=light0,objectType=light,property={direction,-0.2515,-0.3591,0.8988}
	ModifyGizmo modifyObject=light0,objectType=light,property={ambient,0.6667,0.6667,0.6667,1.0000}
	ModifyGizmo modifyObject=light0,objectType=light,property={specular,0.8,0.8,0.8,1.000000}
	ModifyGizmo setDisplayList=0, object=light0
	
	AppendToGizmo attribute shininess={35,1032},name=shininess0
	AppendToGizmo attribute specular={1,1,1,1,1032},name=specular0
	AppendToGizmo attribute blendFunc={770,771},name=blendFunc0	
	// antialiasing, see help topic Printing Gizmo Windows:
	ModifyGizmo setDisplayList=1, opName=enable0, operation=enable, data=2848
	ModifyGizmo setDisplayList=2, opName=enable1, operation=enable, data=2832
	ModifyGizmo setDisplayList=3, attribute=blendFunc0
	
	ModifyGizmo setDisplayList=4, attribute=shininess0
	ModifyGizmo setDisplayList=5, attribute=specular0
	
	AppendToGizmo sphere={1,20,20},name=sphere0
	ModifyGizmo modifyObject=sphere0,objectType=sphere,property={colorType,0}
//	ModifyGizmo setObjectAttribute={sphere0,shininess0}
//	ModifyGizmo setObjectAttribute={sphere0,specular0}
		
	AppendToGizmo Scatter=xyz,name=scatter0
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ markerType,0}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ rotationType,0}
	//ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,2} // std spheres
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ Shape,7}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter, property={ objectName,sphere0} // custom spheres
	//ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ size,1} // single size
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ sizeType,1}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ sizeWave,siz}
	
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ scatterColorType,1}
	ModifyGizmo ModifyObject=scatter0,objectType=scatter,property={ colorWave,cwv}
	ModifyGizmo setDisplayList=6, object=scatter0
	
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""

	// Axes
	make/O/N=(12) Ax, AxT, AxN
	make/O/T/N=(12) AxL
	Ax=1
	AxT=0
	AxN=0
	AxL=""
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo setDisplayList=7, object=axes0
	// set properties of all axes:
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisRange,-1,-1,-1,1,-1,-1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisColor,0,0,0,1}
	//ModifyGizmo modifyObject=axes0,objectType=Axes,property={-1,Clipped,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisLabel,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisLabelDistance,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisLabelRGBA,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,labelBillboarding,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,tickScaling,3}
	if (cmpstr(str2,"Windows")==0)
		ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,fontScaleFactor,2}
	else
		ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,fontScaleFactor,1}
		ModifyGizmo ModifyObject=axes0,objectType=Axes,property={-1,axisLabelScale,1.5}
	endif
	
	// modify axes with ticks and labels:
	AxT[0]=1
	AxT[1]=1
	AxT[3]=1
	AxN[0]=1
	AxN[1]=1
	AxN[3]=1
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,ticks,3}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={3,ticks,3}
	AxL[0]="X"
	AxL[1]="Y"
	AxL[3]="Z"
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={3,axisLabel,1}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={0,axisLabelText,"X"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={1,axisLabelText,"Y"}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={3,axisLabelText,"Z"}
	// turn off axes in front of object, with std orientation:
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 11,visible,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 6,visible,0}
	ModifyGizmo ModifyObject=axes0,objectType=Axes,property={ 2,visible,0}
	Ax[11]=0
	Ax[6]=0
	Ax[2]=0
	
	ModifyGizmo aspectratio=1
	
	ModifyGizmo namedHookStr={GChook,"Chem3Dmodule#GChook"}  // check for resize and kill
	SetWindow $winnam, hook(GCHook2) = Chem3Dmodule#GChook2  // check if moved
	
	ModifyGizmo SETQUATERNION={0.490507,-0.123205,-0.243853,0.827501}
	Getgizmo currotation
	//ModifyGizmo home={-GizmoEulerA, -GizmoEulerB,  -GizmoEulerC} was igor bug
	ModifyGizmo home={GizmoEulerA, GizmoEulerB,  GizmoEulerC}
	
	String/G atomlist=GetWavesDataFolder(atoms,2)
	String/G xyzlist=GetWavesDataFolder(xyz,2)
	
	dispnr+=1
		
	// make bonds as paths
	variable/G BondsOnOff=bonds01
	variable/G BondMult=1
	Chem3DmakeBonds(atoms, xyz)

	//ModifyGizmo showInfo
	//ModifyGizmo infoWindow={540,100,1000,300}
	DoUpdate
	
	GCniceAxislimits() // also makes variables storing coords
	GCmakePanel(winnam) // initial panel, see constants above for default 
	
	GCrangeAllData(winnam) // set to all data + margin for symbols
	GCrangeCubeData(winnam) // set to isotropic axes

	str="ShowTools/W="+winnam
	execute/P/Q/Z str
	
	// set dependency to trigger redraw if coords or elements changed
	variable/G chg
	str="Chem3Dmodule#Chem3DdataChange("+GetWavesDataFolder(xyz,2)+","+GetWavesDataFolder(atoms,2)+")"
	str="SetFormula root:"+winnam+":chg,\""+str+"\""
	execute/P/Q/Z str
	
	if (! bonds01)
		str="DoWindow/F "+winnam
		str+=" ; Chem3Dmodule#Chem3DremoveBonds("+xyzlist+")"
		print str
		execute/P/Q/Z str
	endif
	
	str="DoWindow/F "+winnam
	Execute/P/Q str
	
	GCscatterSaveCoords(winnam) // save the incoming x,y,z coords, set rotation variables
	
	SetDataFolder savedDataFolder
end

//---------------------------------------------------------------//
Function GCniceAxislimits()
	//string str
	variable x1, x2, y1, y2, z1, z2, ii, jj
	jj=0.075
	
//	ModifyGizmo autoscale
//	doupdate
	
	Getgizmo datalimits
	//print GizmoXmin, GizmoXmax, GizmoYmin, GizmoYmax, GizmoZmin, GizmoZmax 
	ii=GizmoXmax - GizmoXmin
	x1=GizmoXmin - jj*ii
	x2=GizmoXmax + jj*ii
	ii=GizmoYmax - GizmoYmin
	y1=GizmoYmin - jj*ii
	y2=GizmoYmax + jj*ii
	ii=GizmoZmax - GizmoZmin
	z1=GizmoZmin - jj*ii
	z2=GizmoZmax + jj*ii
	
//	Getgizmo userboxlimits  
//	ii=GizmoBoxXmax - GizmoBoxXmin
//	x1=GizmoBoxXmin - jj*ii
//	x2=GizmoBoxXmax + jj*ii
//	ii=GizmoBoxYmax - GizmoBoxYmin
//	y1=GizmoBoxYmin - jj*ii
//	y2=GizmoBoxYmax + jj*ii
//	ii=GizmoBoxZmax - GizmoBoxZmin
//	z1=GizmoBoxZmin - jj*ii
//	z2=GizmoBoxZmax + jj*ii
	
	ModifyGizmo scalingoption=0
	ModifyGizmo setOuterBox={x1,x2,y1,y2,z1,z2}
	//print x1,x2,y1,y2,z1,z2
	variable/G GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
	GCxmin=x1
	GCxmax=x2
	GCymin=y1
	GCymax=y2
	GCzmin=z1
	GCzmax=z2
end

// background task and control for checking if gizmos have moved
//---------------------------------------------------------------//
//Function GCcheckMoved(s)		// This is the function that will be called periodically
//	STRUCT WMBackgroundStruct &s
//	
//	GCmovePanel() 
//	//Printf "Task %s called, ticks=%d\r", s.name, s.curRunTicks
//	return 0	// Continue background task
//End
//Function GCStartCheckTask()
//	Variable numTicks = 6		// 60 ticks/sec
//	CtrlNamedBackground GCcheck, period=numTicks, proc=GCcheckMoved
//	CtrlNamedBackground GCcheck, start
//End
//Function GCStopCheckTask()
//	CtrlNamedBackground GCcheck, stop
//End

////---------------------------------------------------------------//
//Function GCmovePanel() // check if gizmo moved, if so update panel position
//	string gwinlist, winstr, wvstr, fullwvstr, recstr, str, topwin
//	variable nwin, ii, jj, kk
//	
//	gwinlist=winlist("*",";","WIN:65536")  // gizmo list, IP7 only
//	nwin=itemsinlist(gwinlist)
//	for (ii=0; ii<nwin; ii+=1) // check all gizmos
//		winstr=stringfromlist(ii, gwinlist)
//		str="root:"+winstr
//		if (datafolderexists(str) ) // it has a folder?
//			String savedDataFolder = GetDataFolder(1)
//			setdatafolder $str
//			if (exists("GCx"))
//				nvar GCx, GCy, GCdx, GCdy
//				getgizmo winpixels
//				if ((GCx!=v_left) ||(GCy!=v_top))
//					//print "moved", winstr, GCx,v_left, GCy, v_top
//					str=winstr+"_GCpanel"
//					killWindow/Z $str
//					//execute/P/Q str
//					GCx=v_left
//					GCy=v_top
//					GCdx=v_right-v_left
//					GCdy=v_top-v_bottom
//					GCmakePanel(winstr)
//				endif
//			endif
//			setdatafolder savedDataFolder
//		endif
//	endfor
//end


//---------------------------------------------------------------//
Function GCmarkAllRefresh() // called by BeforeExperimentSaveHook hook, marks windows as needing refresh
	string gwinlist, winstr, wvstr, fullwvstr, recstr, str, topwin
	variable nwin, ii, jj, kk
	
	//gwinlist=winlist("*",";","WIN:4096")  // XOP windows, IP6
	gwinlist=winlist("*",";","WIN:65536")  // gizmo list, IP7 only
	nwin=itemsinlist(gwinlist)
	for (ii=0; ii<nwin; ii+=1) // check all gizmos
		winstr=stringfromlist(ii, gwinlist)
		str="root:"+winstr
		if (datafolderexists(str) ) // it has a folder?
			String savedDataFolder = GetDataFolder(1)
			setdatafolder $str
			if (exists("GCrefresh"))
				nvar GCrefresh
				GCrefresh=1
			endif
			setdatafolder savedDataFolder
		endif 
	endfor
end

//---------------------------------------------------------------//
Function GCcheckAllRefresh() // called by after compiled hook, mainly to reset windows after expt reloaded
	string gwinlist, winstr, wvstr, fullwvstr, recstr, str, topwin
	variable nwin, ii, jj, kk
	
	//gwinlist=winlist("*",";","WIN:4096")  // XOP windows, IP6
	gwinlist=winlist("*",";","WIN:65536")  // gizmo list, IP7 only
	nwin=itemsinlist(gwinlist)
	for (ii=0; ii<nwin; ii+=1) // check all gizmos
		winstr=stringfromlist(ii, gwinlist)
		str="root:"+winstr
		if (datafolderexists(str) ) // it has a folder?
			String savedDataFolder = GetDataFolder(1)
			setdatafolder $str
			if (exists("GCrefresh"))
				nvar GCrefresh
				if (GCrefresh)
					Dowindow/F $winstr
					str=winstr+"_GCpanel"
					killWindow/Z $str
					GCmakePanel(winstr)
					GCrefresh=0
					ModifyGizmo/N=$winstr SETQUATERNION={0.490507,-0.123205,-0.243853,0.827501}
					Getgizmo/N=$winstr currotation
					ModifyGizmo/N=$winstr home={GizmoEulerA, GizmoEulerB,  GizmoEulerC}
					nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
					ModifyGizmo/N=$winstr scalingoption=0
					ModifyGizmo/N=$winstr setOuterBox={GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax}
					Dowindow/hide=1 $str
				endif
			endif
			setdatafolder savedDataFolder
		endif 
	endfor
end

//---------------------------------------------------------------//
Function GCmovePanel(gname) // check if gizmo moved, if so update panel position
	string gname

	string  winstr, str
	variable xx, yy, dx

	getgizmo/N=$gname winpixels
	xx=v_left
	yy=v_top
	winstr=gname+"_GCpanel"
	str="root:"+gname+":GCpanelMode"
	nvar GCpanelMode=$str
	str="root:"+gname+":GCpdx"
	nvar GCpdx=$str
	str="root:"+gname+":GCpdy"
	nvar GCpdy=$str
	str="root:"+gname+":GCpscl"
	nvar GCpscl=$str
	
//	str=igorinfo(2)
//	if (cmpstr(str,"Windows")==0)
//		//dx=screenresolution/Panelresolution(str)
//		xx=(v_right+8) //*dx
//	else
//		xx=v_right
//	endif
	GetWindow $gname, wsizeouter // this includes frames on windows
	xx=xx+(v_right-v_left)

//print "move panel", xx, yy, xx+GCpdx, yy+GCpdy
	MoveWindow/W=$winstr xx, yy, xx+GCpdx*GCpscl, yy+GCpdy*GCpscl 
end

//---------------------------------------------------------------//
Function GCmakePanel(GizmoName) 
	string GizmoName
	string str, str2
	variable xx, yy, dx, dy, defPres
	
	String savedDataFolder = GetDataFolder(1)
	
	str="root:"+GizmoName
	setdatafolder $str
	
	defPres=panelresolution("") // save existing panel setting
	nvar GCpscl
	xx=round(GCpscl*screenresolution)
	str="SetIgorOption PanelResolution = "+num2str(xx) // use full resln
	execute str
	
	if (!exists("GCpdx"))
		variable/G GCpdx, GCpdy, GCpanelMode
		GCpanelMode=GCdefaultpanelmode // <0:minimized, 1:default, 2:axes, 3:ranges, 4:rotation
	else
		nvar GCpdx, GCpdy, GCpanelMode
	endif
	nvar GCzoom, GCaxes
	
	str=GizmoName+"_GCpanel"
	Getgizmo/N=$GizmoName winpixels

	str2=igorinfo(2) // also needed below
//	if (cmpstr(str2,"Windows")==0)
//		xx=(v_right+8) //*dx
//	else
//		xx=v_right
//	endif
	xx=v_left
	yy=v_top
	if (GCpanelMode<0)
			GCpdx=30 
			GCpdy=240
	else
			GCpdx=200  // a default size, may be resized by individual panels
			GCpdy=240  
	endif
	
	GetWindow $GizmoName, wsizeouter // this includes frames on windows
	xx=xx+(v_right-v_left)

	//print "new panel", xx, yy, xx+GCpdx, yy+GCpdy
	newpanel/K=1/N=$str/W=(xx, yy, xx+GCpdx*GCpscl, yy+GCpdy*GCpscl )
	ModifyPanel cbRGB=(BGR, BGG, BGB)
	DefaultGUIControls/W=$str native
	
	// permanent controls
	Button minimizebutton,pos={5 ,2 },size={20 ,20 },proc=GCMinimizeButtonProc,title=" < "
	Button minimizebutton,help={"Shrink/Expand"},font="Tahoma",fStyle=1
	Button minimizebutton,fColor=(65535,49151,49151), thumbColor=(65535,21845,0)
	
	if (GCpanelMode<=0)
		Button minimizebutton,title=" > "
	else
		Button minimizebutton,title=" < "
	endif
	
	SetDrawLayer UserBack
	SetDrawEnv fsize= 9
	DrawText 3,193,"Zoom"
	
	Slider zoomslider,pos={6.00 ,28.0 },size={13.00 ,150.00 },proc=GCzoomSliderProc
	Slider zoomslider,help={"Zoom"},limits={1.75,0.2,0},value= GCzoom, ticks= -1
	Slider zoomslider,labelBack=(65535,65535,65535),  focusRing=0
	if (cmpstr(str2,"Windows")==0)
		Slider zoomslider, side= 1
	else
		Slider zoomslider, side= 0
	endif
	Slider zoomslider labelBack=(49151,65535,49151)
	
	CheckBox Axescheck,pos={7.00,206.00 },size={16.00 ,15.00},proc=GCAxesCheckProc,title=""
	CheckBox Axescheck,font="Tahoma",variable=GCaxes
	SetDrawEnv fsize= 9
	DrawText 4,233,"Axes"
	
	SetWindow kwtopwin, hook(GCpanelHook) = Chem3Dmodule#GCpanelhook
	
	// page-dependent controls
	if (GCpanelMode==1) // object mode
		GCobjectMode(GizmoName) 
	elseif (GCpanelMode==2) // axes mode
		GCaxesMode(GizmoName) 
	elseif (GCpanelMode==3) // range mode	
		GCrangeMode(GizmoName)
	elseif (GCpanelMode==4) // rotate mode	
		GCrotateMode(GizmoName)
	else
		GCmovepanel(GizmoName)
	endif
	
	setdatafolder savedDataFolder
	
	str="SetIgorOption PanelResolution = "+num2str(defPres) // restore panel resln
	execute str
end

//------------------------------//
Function GCrangeButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			str="KillWindow/Z "+str // qeue panel to be killed
			Execute/P/Q str
			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			str="root:"+gizname+":GCpanelMode"
			nvar GCpanelMode=$str
			GCpanelMode=3 // <0:minimized, 1:default, 2:axes, 3:ranges

			str=""
			str="Chem3Dmodule#GCmakePanel(\""+GizName+"\") " // qeue remake of panel
			Execute/P/Q str
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrangeMode(GizmoName)  
	string GizmoName
	
	String savedDataFolder = GetDataFolder(1)
	string str, ctrl, str2
	variable ypos, dy, ii, jj
	
	str="root:"+GizmoName
	setdatafolder $str
	
	nvar GCpdx, GCpdy
	nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
	
	GCpdx=150
	GCpdy=240
	GCmovepanel(GizmoName)
	
	Button ObjectButton,pos={30,2.50},size={50.00,15},title="Chem3D",font="Tahoma"
	Button ObjectButton,fSize=11 , proc=GCobjectButtonProc
	
	Button AxesButton,pos={80,2.50},size={30.00,15},title="Axes",font="Tahoma"
	Button AxesButton,fSize=11, proc=GCAxesButtonProc
	
	Button RotateButton,pos={110,2.50},size={37,15},title="Rotate",font="Tahoma"
	Button RotateButton,fSize=11, proc=GCRotateButtonProc
	
	TitleBox Chem3Dpaneltitle,pos={61,22},size={61.00,12.50},title="Range Panel"
	TitleBox Chem3Dpaneltitle,frame=0, fstyle=1
		
	ypos=36
	
	SetVariable setvarXL,pos={30.75,ypos},size={57.75,13.50},title="X:",fSize=9,  focusRing=0
	SetVariable setvarXL,limits={-inf,inf,0}, value=GCxmin, fstyle=0, proc=GCrangeSetVarProc
	SetVariable setvarXH,pos={91.50,ypos},size={54.00,13.50},title="-",fSize=9,  focusRing=0
	SetVariable setvarXH,limits={-inf,inf,0}, value=GCxmax, proc=GCrangeSetVarProc
	ypos+=18
	Slider Xrangeslider,pos={30.50,ypos},size={115.00,10.50},  focusRing=0
	Slider Xrangeslider,labelBack=(65535,65535,65535),proc=GCrangeSliderProc
	Slider Xrangeslider,limits={-1,1,0},vert= 0,ticks= 0,live=1, value=0
	str2=igorinfo(2)
	if (cmpstr(str2,"Windows")==0)
		Slider Xrangeslider, side=1
	else
		Slider Xrangeslider, side=0
	endif
	ypos+=15
	TitleBox Xrangetitle,pos={49.50,ypos},size={61.00,12.50},title="Range (shift=offset)"
	TitleBox Xrangetitle,frame=0
	
	ypos+=20	
	SetVariable setvarYL,pos={30.75,ypos},size={57.75,13.50},title="Y:",fSize=9,  focusRing=0
	SetVariable setvarYL,limits={-inf,inf,0}, value=GCymin, fstyle=0, proc=GCrangeSetVarProc
	SetVariable setvarYH,pos={91.50,ypos},size={54.00,13.50},title="-",fSize=9,  focusRing=0
	SetVariable setvarYH,limits={-inf,inf,0}, value=GCymax, proc=GCrangeSetVarProc
	ypos+=18
	Slider Yrangeslider,pos={30.50,ypos},size={115.00,10.50},  focusRing=0
	Slider Yrangeslider,labelBack=(65535,65535,65535),proc=GCrangeSliderProc
	Slider Yrangeslider,limits={-1,1,0},vert= 0,ticks= 0,live=1, value=0
	if (cmpstr(str2,"Windows")==0)
		Slider Yrangeslider, side=1
	else
		Slider Yrangeslider, side=0
	endif
	ypos+=15
	TitleBox Yrangetitle,pos={49.50,ypos},size={61.00,12.50},title="Range (shift=offset)"
	TitleBox Yrangetitle,frame=0
	
	ypos+=20	
	SetVariable setvarZL,pos={30.75,ypos},size={57.75,13.50},title="Z:",fSize=9,  focusRing=0
	SetVariable setvarZL,limits={-inf,inf,0}, value=GCzmin, fstyle=0, proc=GCrangeSetVarProc
	SetVariable setvarZH,pos={91.50,ypos},size={54.00,13.50},title="-",fSize=9,  focusRing=0
	SetVariable setvarZH,limits={-inf,inf,0}, value=GCzmax, proc=GCrangeSetVarProc
	ypos+=18
	Slider Zrangeslider,pos={30.50,ypos},size={115.00,10.50},  focusRing=0
	Slider Zrangeslider,labelBack=(65535,65535,65535),proc=GCrangeSliderProc
	Slider Zrangeslider,limits={-1,1,0},vert= 0,ticks= 0,live=1, value=0
	if (cmpstr(str2,"Windows")==0)
		Slider Zrangeslider, side=1
	else
		Slider Zrangeslider, side=0
	endif
	ypos+=15
	TitleBox Zrangetitle,pos={49.50,ypos},size={61.00,12.50},title="Range (shift=offset)"
	TitleBox Zrangetitle,frame=0
	
	ypos+=20
	Button RangeAllButton,pos={32,ypos},size={50.00,15},title="All Data",font="Tahoma"
	Button RangeAllButton,fSize=11, proc=GCrangeAllButtonProc
	
	Button RangeCubeButton,pos={90,ypos},size={50.00,15},title="Cube",font="Tahoma"
	Button RangeCubeButton,fSize=11, proc=GCrangeCubeButtonProc
	
	ypos+=25
//	Button RangeNotesButton,pos={60,ypos},size={50.00,15},title="Help",font="Tahoma"
//	Button RangeNotesButton,fSize=11 , proc=GCrangeNotesButtonProc

	setdatafolder savedDataFolder
end

//------------------------------//
Function GCrangeCubeButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	variable ii
	string str, gname

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			ii=strsearch(str,"_GCpanel",0,2)
			gname=str[0,ii-1]
			GCrangeCubeData(gname)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrangeCubeData(gname) 
	string gname
	string str
	variable ii, jj, midx, midy, midz, dx, dy, dz, dmax
				
	String savedDataFolder = GetDataFolder(1)
	str="root:"+GName
	setdatafolder $str
	nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
	svar xyzlist
	wave xyz=$xyzlist
	
	dx=GCxmax-GCxmin
	midx=GCxmin+dx/2	

	dy=GCymax-GCymin
	midy=GCymin+dy/2	

	dz=GCzmax-GCzmin
	midz=GCzmin+dz/2	
	
	dmax=max(dx, dy)
	dmax=max(dmax, dz)
	
	GCxmin=midx-dmax/2
	GCxmax=midx+dmax/2
	
	GCymin=midy-dmax/2
	GCymax=midy+dmax/2
	
	GCzmin=midz-dmax/2
	GCzmax=midz+dmax/2
	
	ModifyGizmo/N=$gname setOuterBox={GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax}
	
	//GCrangeAfterChangeProc(gname)  now in gizmo hook
	
	setdatafolder savedDataFolder
End

//------------------------------//
Function GCrangeAllButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	variable ii
	string str, gname

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			ii=strsearch(str,"_GCpanel",0,2)
			gname=str[0,ii-1]
			GCrangeAllData(gname)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrangeAllData(gname) 
	string gname
	string str
	variable ii, jj, delta, mid, mx, mn
	variable margin 
	variable mindelta
				
	String savedDataFolder = GetDataFolder(1)
	str="root:"+GName
	setdatafolder $str
	nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax, AtomScale
	svar xyzlist
	wave xyz=$xyzlist
	
	margin= 1 + 0.35*AtomScale  // extra for symbols
	mindelta=2.5
	
	// X axis
	mx=xyz[0][0]
	mn=xyz[0][0]
	for (ii=1; ii<dimsize(xyz,0); ii+=1)
		mx= max(xyz[ii][0],mx)
		mn= min(xyz[ii][0],mn)
	endfor
	delta=mx-mn
	mid=mn+delta/2
	delta=max(delta, mindelta)
	delta*=margin  
	GCxmin=mid-delta/2
	GCxmax=mid+delta/2
	
	// Y axis
	mx=xyz[0][1]
	mn=xyz[0][1]
	for (ii=1; ii<dimsize(xyz,0); ii+=1)
		mx= max(xyz[ii][1],mx)
		mn= min(xyz[ii][1],mn)
	endfor
	delta=mx-mn
	mid=mn+delta/2
	delta=max(delta, mindelta)
	delta*=margin  
	GCymin=mid-delta/2
	GCymax=mid+delta/2
	
	// Z axis
	mx=xyz[0][2]
	mn=xyz[0][2]
	for (ii=1; ii<dimsize(xyz,0); ii+=1)
		mx= max(xyz[ii][2],mx)
		mn= min(xyz[ii][2],mn)
	endfor
	delta=mx-mn
	mid=mn+delta/2
	delta=max(delta, mindelta)
	delta*=margin  
	GCzmin=mid-delta/2
	GCzmax=mid+delta/2
	
	ModifyGizmo/N=$gname setOuterBox={GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax}
	
	//GCrangeAfterChangeProc(gname)  now in gizmo hook
	
	setdatafolder savedDataFolder
End

//------------------------------//
Function GCresetAxesLimits(gname)
	string gname
	string str, wrec
	variable ii, jj
	
	String savedDataFolder = GetDataFolder(1)
	str="root:"+gname
	setdatafolder $str
	nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
	
	wrec=winrecreation(gname,0)
	ii=strsearch(wrec,"setOuterBox={",0,2)
	if (ii>=0)
		ii+=12
		jj=strsearch(wrec,"}",ii,2)
		str=wrec[ii+1, jj-1]
		//print str
		//print itemsinlist(str, ",")
		GCxmin=str2num(stringfromlist(0, str,","))
		GCxmax=str2num(stringfromlist(1, str,","))
		GCymax=str2num(stringfromlist(2, str,","))
		GCymax=str2num(stringfromlist(3, str,","))
		GCzmax=str2num(stringfromlist(4, str,","))
		GCzmax=str2num(stringfromlist(5, str,","))
	endif
	setdatafolder savedDataFolder
end

//------------------------------//
Function GCrangeSetVarProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva
	string str, gizname
	variable ii

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			str=sva.win
			ii=strsearch(str,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			String savedDataFolder = GetDataFolder(1)
			str="root:"+GizName
			setdatafolder $str
			nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
			ModifyGizmo/N=$gizname scalingoption=0
			ModifyGizmo/N=$gizname setOuterBox={GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax}
			setdatafolder savedDataFolder
			
			//GCrangeAfterChangeProc(gizname)  now in gizmo hook
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCallRangeSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname
	variable ii, curval, delta, midx, midy, midz

	ii=sa.eventCode
	curval = sa.curval
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if (( ii & 1 ) || ( ii & 4 ) || ( ii & 2 ) ) // mouse up or click in scale
				//print "released at", curval
				str=sa.win
				ii=strsearch(sa.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				
				String savedDataFolder = GetDataFolder(1)
				str="root:"+GizName
				setdatafolder $str
				nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
					
				delta=GCxmax-GCxmin
				midx=GCxmin+delta/2
				
				ii=GCymax-GCymin
				midy=GCymin+ii/2
				delta=max(delta, ii)
				
				ii=GCzmax-GCzmin
				midz=GCzmin+ii/2
				delta=max(delta, ii)
				delta=delta/2 *1.2^curval  // scaling here
				
				GCxmin=midx-delta
				GCxmax=midx+delta
				
				GCymin=midy-delta
				GCymax=midy+delta
				
				GCzmin=midz-delta
				GCzmax=midz+delta
					
				ModifyGizmo/N=$gizname scalingoption=0
				ModifyGizmo/N=$gizname setOuterBox={GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax}
	
				str=sa.ctrlName
				Slider $str,value=0 // spring back
				
				setdatafolder savedDataFolder
			endif
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrangeSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname
	variable ii, curval, delta, mid
	
	ii=sa.eventCode
	curval = sa.curval
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if (( ii & 1 ) ||( ii & 4 ) || ( ii & 2 ) ) // mouse up or click in scale
				//print "released at", curval
				str=sa.win
				ii=strsearch(sa.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				
				String savedDataFolder = GetDataFolder(1)
				str="root:"+GizName
				setdatafolder $str
				nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
				
				strswitch(sa.ctrlName)	
					case "Xrangeslider":	
						delta=GCxmax-GCxmin
						if (sa.eventmod & 2) // shift
							GCxmin+=delta*curval/4
							GCxmax=GCxmin+delta
						else
							mid=GCxmin+delta/2
							delta*=2^curval
							GCxmin=mid-delta/2
							GCxmax=mid+delta/2
						endif
						break
					case "Yrangeslider":	
						delta=GCymax-GCymin
						if (sa.eventmod & 2) // shift
							GCymin+=delta*curval/4
							GCymax=GCymin+delta
						else
							mid=GCymin+delta/2
							delta*=2^curval
							GCymin=mid-delta/2
							GCymax=mid+delta/2
						endif
						break
					case "Zrangeslider":	
						delta=GCzmax-GCzmin
						if (sa.eventmod & 2) // shift
							GCzmin+=delta*curval/4
							GCzmax=GCzmin+delta
						else
							mid=GCzmin+delta/2
							delta*=2^curval
							GCzmin=mid-delta/2
							GCzmax=mid+delta/2
						endif
						break
				endswitch
				ModifyGizmo/N=$gizname scalingoption=0
				ModifyGizmo/N=$gizname setOuterBox={GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax}
	
				str=sa.ctrlName
				Slider $str,value=0 // spring back
				
				//GCrangeAfterChangeProc(gizname)  now in gizmo hook
				
				setdatafolder savedDataFolder
			endif
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrangeAfterChangeProc(gname) // perform actions after a range change
	string gname								// called from gizmo hook 

	Chem3DatomScale(gname)
End
	

//------------------------------//
Function GCobjectButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			str="KillWindow/Z "+str // qeue panel to be killed
			Execute/P/Q str
			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			str="root:"+gizname+":GCpanelMode"
			nvar GCpanelMode=$str
			GCpanelMode=1 // <0:minimized, 1:default, 2:axes, 3:ranges

			str=""
			str="Chem3Dmodule#GCmakePanel(\""+GizName+"\") " // qeue remake of panel
			Execute/P/Q str
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
// this proc part of GizmoControl, but must be adapted to objects displayed
// => is not as general as other parts of GC
Function GCobjectMode(GizmoName)  
	string GizmoName
	
	String savedDataFolder = GetDataFolder(1)
	string str, str2,  ctrl
	variable ypos, dy, ii, jj
	
	str="root:"+GizmoName
	setdatafolder $str
	
	nvar GCpdx, GCpdy
	GCpdx=135
	GCpdy=240
	GCmovepanel(GizmoName) // ensures positioning wrt gizmo window
	
	Button AxesButton,pos={30,2.50},size={30.00,15},title="Axes",font="Tahoma"
	Button AxesButton,fSize=11, proc=GCAxesButtonProc
	
	Button RangeButton,pos={61,2.50},size={46.00,15},title="Ranges",font="Tahoma"
	Button RangeButton,fSize=11, proc=GCrangeButtonProc
	
	Button RotateButton,pos={108,2.50},size={25,15},title="Rot",font="Tahoma"
	Button RotateButton,fSize=11, proc=GCRotateButtonProc
	
	TitleBox Chem3Dpaneltitle,pos={45,22},size={61.00,12.50},title="Chem3D Panel"
	TitleBox Chem3Dpaneltitle,frame=0, fstyle=1
		
	ypos=36
	
	Button BondToggleButton,pos={48.5,ypos},size={50.00,15.00},title="Bonds 0/1",font="Tahoma"
	Button BondToggleButton,fSize=9, proc=Chem3DbondToggleButtonProc
	
	Slider BondThickslider,pos={30.50,55},size={100.00,10.50}, focusRing=0
	Slider BondThickslider,labelBack=(65535,65535,65535), proc=Chem3DbondThickSliderProc
	Slider BondThickslider,limits={0,bondthickmax,0},variable=BondMult,vert= 0,ticks= 0
	str2=igorinfo(2)
	if (cmpstr(str2,"Windows")==0)
		Slider BondThickslider, side=1
	else
		Slider BondThickslider, side=0
	endif
	TitleBox bondthicktitle,pos={45,68},size={61.00,12.50},title="Bond Thickness"
	TitleBox bondthicktitle,frame=0, fSize=9
	
	Slider AtomScaleSlider,pos={30.50,85},size={100.00,10.50}, focusRing=0
	Slider AtomScaleSlider,labelBack=(65535,65535,65535), proc=Chem3DAtomScaleSliderProc
	Slider AtomScaleSlider,limits={minatom,maxatom,0},variable=AtomScale,vert= 0,ticks= 0
	if (cmpstr(str2,"Windows")==0)
		Slider AtomScaleSlider, side=1
	else
		Slider AtomScaleSlider, side=0
	endif
	SetDrawLayer Overlay
	DrawLine 125,82.5,125,99  // marks vdW position
	SetDrawLayer UserBack
	
	TitleBox atomscaletitle,pos={57,98},size={61.00,12.50},title="Atom Size"
	TitleBox atomscaletitle,frame=0, fSize=9
	
	TitleBox vdWtitle,pos={117.00,100},size={15.00,9.00},title="vdW"
	TitleBox vdWtitle,fSize=8,frame=0
	
	ypos=115
	Button AtomNumbersToggleButton,pos={34,ypos},size={85.00,17.00},font="Tahoma"
	Button AtomNumbersToggleButton,fSize=9, proc=Chem3DnumbersButtonProc
	nvar AtomMode
	if (AtomMode==0)
		Button AtomNumbersToggleButton,title="Atom Numbers"  
		TitleBox vdWtitle disable=0 
	else
		Button AtomNumbersToggleButton,title="Show Atoms" 
		TitleBox vdWtitle disable=1
		TitleBox atomscaletitle,pos={48,103},title="Number Size"
	endif
	
	ypos=140
	Slider Rangeslider,pos={30.50,ypos},size={100.00,10.50}, focusRing=0
	Slider Rangeslider,labelBack=(65535,65535,65535),proc=GCallRangeSliderProc
	Slider Rangeslider,limits={-1,1,0},vert= 0,ticks= 0,live=1, value=0
	if (cmpstr(str2,"Windows")==0)
		Slider Rangeslider, side=1
	else
		Slider Rangeslider, side=0
	endif
	
	ypos=153
	TitleBox rangetitle,pos={60.00,ypos},size={55.00,9.00},title="Scale Axes"
	TitleBox rangetitle,fSize=9,frame=0
	
	Slider VpanSlider,pos={37.50,157.00},size={10.50,54.00}
	Slider VpanSlider,labelBack=(65535,65535,65535), proc=GCpanSliderProc
	Slider VpanSlider,limits={0.1,-0.1,0},value= 0,ticks= 0
	
	Slider HpanSlider,pos={53.50,178.00},size={76.50,10.50}
	Slider HpanSlider,labelBack=(65535,65535,65535), proc=GCpanSliderProc
	Slider HpanSlider,limits={0.1,-0.1,0},value= 0,vert= 0,ticks= 0
	
	if (cmpstr(str2,"Windows")==0)
		Slider VpanSlider, side=1
		Slider HpanSlider, side=1
	else
		Slider VpanSlider, side=0
		Slider HpanSlider, side=0
	endif
	Slider VpanSlider size={12,54},labelBack=(49151,65535,65535)
	Slider HpanSlider size={76,10},labelBack=(49151,65535,65535)
	
	TitleBox movetitle,pos={57.00,194.00},size={65.50,12.50},title="Move in window"
	TitleBox movetitle, frame=0
	TitleBox movetitle fColor=(1,52428,52428)

	ypos=215
	Button ExportButton,pos={30,ypos},size={45.00,17.00},title="Save 2D",font="Tahoma"
	Button ExportButton,fSize=9, proc=GCsave2DButtonProc
	
	//ypos=190
	Button GraphButton,pos={79,ypos},size={50.00,17.00},title="Graph 2D",font="Tahoma"
	Button GraphButton,fSize=9, proc=GCgraph2DButtonProc

	setdatafolder savedDataFolder
end

//------------------------------//
Function GCsave2DButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname= str[0,ii-1]
			Dowindow/F $gizname
			ii=2*screenresolution
			SavePICT/WIN=$gizname/E=-5/B=(ii)
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCgraph2DButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string wstr, str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			
			ii=2*screenresolution
			SavePICT/O/WIN=$gizname/E=-5/B=(ii)/P=igoruserfiles as "GCtempImg.png"
			ImageLoad/P=IgorUserFiles/T=rpng/O/Q/N=Chem3DtempImg "GCtempImg.png"
			DeleteFile/Z/P=IgorUserFiles "GCtempImg.png"	
			wave Chem3DtempImg
			for (ii=1; ii<100; ii+=1)
				wstr=gizname+"snapshot"+num2str(ii)
				dowindow $wstr
				if (!v_flag)
					break
				endif
			endfor
			Newimage/K=1 Chem3DtempImg
			Dowindow/C $wstr
			ModifyGraph tick=3,nticks=0,noLabel=2
			ModifyGraph mirror=2
			ModifyGraph margin(left)=17,margin(bottom)=14,margin(top)=15,margin(right)=14
			//ModifyGraph axRGB=(65535,65535,65535)
			SetWindow $wstr, hook(GC2Dhook) = GC2Dhook
			str=gizname+"img"+num2str(ii)
			Rename Chem3DtempImg, $str
			SetWindow $wstr, note=str // save name of wave in window note for hook fn
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//---------------------------------------------------------------//
Function GC2Dhook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	string str

	//print "2dhook", s.eventCode
	switch(s.eventCode)
		case 0:	// Activate    
			break
		case 1:	// Deactivate
			break
		case 2:	// Kill
			str=s.winname
			GetWindow $str, note
			str="killwaves/Z "+s_value // kill wave if window closed
			execute/P/Q str
			break
		case 17: 	// kill vote- must be present as precursor to kill
			break 
	endswitch
	
	return hookResult		// 0 if nothing done, else 1
End

//------------------------------//
Function Chem3DnumbersButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			str="root:"+gizname+":AtomMode"
			nvar AtomMode=$str
			str="root:"+gizname+":AtomScale"
			nvar AtomScale=$str
			str="root:"+gizname+":AtomSizes"
			wave AtomSizes=$str			
			if (AtomMode==0)
				str="root:"+gizname+":AtomNumbers"
				ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ Shape,8}
				ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ TextWave,$str}
				AtomSizes=1
				str="root:"+gizname+":AtomScaleSave"
				variable/G $str=AtomScale
				AtomScale=1
				AtomSizes=1
				//ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ sizeType,0}
				//ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ size,1.5} // single size
				Button AtomNumbersToggleButton,title="Show Atoms",font="Tahoma"//, userdata="on"
				TitleBox atomscaletitle,pos={48,100},size={61.00,12.50},title="Number Size"
				TitleBox vdWtitle disable=1
				Slider AtomScaleSlider,limits={0.2,5,0}
				AtomMode=1
				Chem3DatomScale(gizname)
			else
				str="root:"+gizname+":AtomSizes"
				ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ Shape,7}
				ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter, property={ objectName,sphere0} // custom spheres
				str="root:"+gizname+":AtomScaleSave"
				nvar AtomScaleSave=$str
				AtomScale=AtomScaleSave
				//ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ sizeType,1}
				//ModifyGizmo/N=$gizname ModifyObject=scatter0,objectType=scatter,property={ sizeWave,$str}
				Button AtomNumbersToggleButton,title="Atom Numbers",font="Tahoma"//, userdata="off"
				TitleBox atomscaletitle,pos={58,100},size={61.00,12.50},title="Atom Size"
				TitleBox vdWtitle disable=0
				Slider AtomScaleSlider,limits={minatom,maxatom,0}
				AtomMode=0
				Chem3DatomScale(gizname)
			endif
			DoWindow/F $gizname // these 3 lines cause an update of the slider position
			str=ba.win
			DoWindow/F $str
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrotateMode(GizmoName)  
	string GizmoName
	
	String savedDataFolder = GetDataFolder(1)
	string str, str2,  ctrl
	variable ypos, dy, ii, jj
	
	str="root:"+GizmoName
	setdatafolder $str
	
	nvar GCpdx, GCpdy
	nvar GcScXrot, GcScYrot, GcScZrot
	GCpdx=155
	GCpdy=240
	GCmovepanel(GizmoName) // ensures positioning wrt gizmo window
	
	Button ObjectButton,pos={30,2.5},size={50.00,15},title="Chem3D",font="Tahoma"
	Button ObjectButton,fSize=11 , proc=GCobjectButtonProc
	
	Button AxesButton,pos={80,2.5},size={30.00,15},title="Axes",font="Tahoma"
	Button AxesButton,fSize=11, proc=GCAxesButtonProc
	
	Button RangeButton,pos={110,2.5},size={43.00,15},title="Ranges",font="Tahoma"
	Button RangeButton,fSize=11, proc=GCrangeButtonProc
	
	TitleBox Paneltitle,pos={45,30},size={61.00,12.50},title="Coordinate Rotation"
	TitleBox Paneltitle,frame=0, fstyle=1

	// X rotation
	ypos=55
	Slider RotXslider,pos={30.50,ypos},size={118,10.50}, focusRing=0
	Slider RotXslider,labelBack=(65535,65535,65535), proc=GCscatterRotSliderProc
	Slider RotXslider,limits={-pi, pi,0},variable=GcScXrot,vert= 0,ticks= 0
	str2=igorinfo(2)
	if (cmpstr(str2,"Windows")==0)
		Slider RotXslider, side=1
	else
		Slider RotXslider, side=0
	endif
	TitleBox rotXtitle,pos={60,ypos+15},size={61.00,12.50},title="Rotate around X"
	TitleBox rotXtitle,frame=0, fSize=9
	
	// Y rotation
	ypos=95
	Slider RotYslider,pos={30.50,ypos},size={118,10.50}, focusRing=0
	Slider RotYslider,labelBack=(65535,65535,65535), proc=GCscatterRotSliderProc
	Slider RotYslider, limits={pi,-pi,0},variable=GcScYrot,vert= 0,ticks= 0
	str2=igorinfo(2)
	if (cmpstr(str2,"Windows")==0)
		Slider RotYslider, side=1
	else
		Slider RotYslider, side=0
	endif
	TitleBox rotYtitle,pos={60,ypos+15},size={61.00,12.50},title="Rotate around Y"
	TitleBox rotYtitle,frame=0, fSize=9
	
	// Z rotation
	ypos=135
	Slider RotZslider,pos={30.50,ypos},size={118,10.50}, focusRing=0
	Slider RotZslider,labelBack=(65535,65535,65535), proc=GCscatterRotSliderProc
	Slider RotZslider,limits={pi, -pi,0},variable=GcScZrot,vert= 0,ticks= 0
	str2=igorinfo(2)
	if (cmpstr(str2,"Windows")==0)
		Slider RotZslider, side=1
	else
		Slider RotZslider, side=0
	endif
	TitleBox rotZtitle,pos={60,ypos+15},size={61.00,12.50},title="Rotate around Z"
	TitleBox rotZtitle,frame=0, fSize=9
	
	ypos=175
	Button SaveRefButton,pos={40,ypos},size={100,15},title="Save as Reference",font="Tahoma"
	Button SaveRefButton,fSize=11, proc=GCsaveRefButtonProc
	
	ypos=200
	Button RestoreRefButton,pos={40,ypos},size={100,15},title="Restore Reference",font="Tahoma"
	Button RestoreRefButton,fSize=11, proc=GCrestoreRefButtonProc

	setdatafolder savedDataFolder
end

//------------------------------//
Function GCsaveRefButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]

			GCscatterSaveCoords(gizname)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCrestoreRefButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]

			GCscatterRestoreCoords(gizname)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCscatterRotSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname
	variable ii, npts, r1
	
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				str=sa.win
				ii=strsearch(sa.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				
				GCscatterRotCoords(gizname)

			endif
			break
	endswitch

	return 0
End

//------------------------------//
Function GCscatterRotCoords(gname)
	string gname
	
	string str
	String savedDataFolder = GetDataFolder(1)
	str="root:"+gname
	setdatafolder str
	
	svar XYZlist
	wave wv3d=$XYZlist
	
	str="root:"+gname+":"+parsefilepath(0,XYZlist, ":", 1, 0) + "_save"
	wave wv3dref= $str
	
	variable dx, dy, dz
	
	make/D/free/N=(3,3) rotmat
	
	// slider rotation values are applied sequentially to ref xyz
	nvar GcScXrot, GcScYrot, GcScZrot
	rotmat[][0]={cos(GcScXrot), 0, sin(GcScXrot)}
	rotmat[][1]={0, 1, 0}
	rotmat[][2]={-1*sin(GcScXrot), 0 , cos(GcScXrot)}
	matrixop/O wv3D=(rotmat x wv3dref^t)^t
	
	rotmat[][0]={1, 0, 0}
	rotmat[][1]={0, cos(GcScYrot), -1*sin(GcScYrot)}
	rotmat[][2]={0, sin(GcScYrot),  cos(GcScYrot)}
	matrixop/O wv3D=(rotmat x wv3d^t)^t

	rotmat[][0]={cos(GcScZrot), -1*sin(GcScZrot), 0}
	rotmat[][1]={sin(GcScZrot), cos(GcScZrot), 0}
	rotmat[][2]={0, 0, 1}
	matrixop/O wv3D=(rotmat x wv3d^t)^t

	setdatafolder savedDataFolder
end

//------------------------------//
Function GCscatterSaveCoords(gname)
	string gname
	variable axnr, ang
	
	string str
	String savedDataFolder = GetDataFolder(1)
	str="root:"+gname
	setdatafolder str
	
	svar XYZlist
	str="root:"+gname+":"+parsefilepath(0,XYZlist, ":", 1, 0) + "_save"
	Duplicate/O $XYZlist, $str
	
	variable/G GcScXrot, GcScYrot, GcScZrot
	GcScXrot=0
	GcScYrot=0
	GcScZrot=0
	
	setdatafolder savedDataFolder
End

//------------------------------//
Function GCscatterRestoreCoords(gname)
	string gname
	variable axnr, ang
	
	string str
	String savedDataFolder = GetDataFolder(1)
	str="root:"+gname
	setdatafolder str
	
	svar XYZlist
	str="root:"+gname+":"+parsefilepath(0,XYZlist, ":", 1, 0) + "_save"
	Duplicate/O $str, $XYZlist
	
	variable/G GcScXrot, GcScYrot, GcScZrot
	GcScXrot=0
	GcScYrot=0
	GcScZrot=0
	
	setdatafolder savedDataFolder
End

//------------------------------//
Function GCRotateButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			str="KillWindow/Z "+str // qeue panel to be killed
			Execute/P/Q str
			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			str="root:"+gizname+":GCpanelMode"
			nvar GCpanelMode=$str
			GCpanelMode=4 // <0:minimized, 1:default, 2:axes, 3:ranges, 4:rotate

			str=""
			str="Chem3Dmodule#GCmakePanel(\""+GizName+"\") " // qeue remake of panel
			Execute/P/Q str
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


//------------------------------//
Function GCaxesMode(GizmoName) 
	string GizmoName
	
	String savedDataFolder = GetDataFolder(1)
	string str, ctrl
	variable ypos, dy, ii, jj, dyT,dyC
	
	str="root:"+GizmoName
	setdatafolder $str
	wave Ax, AxT, AxL, AxN
	nvar GCpdx, GCpdy
	
	GCpdx=200
	GCpdy=240
	GCmovepanel(GizmoName)
	
	Button ObjectButton,pos={30,2.50},size={50.00,15},title="Chem3D",font="Tahoma"
	Button ObjectButton,fSize=11 , proc=GCobjectButtonProc
	
	Button RangeButton,pos={80,2.50},size={46.00,15},title="Ranges",font="Tahoma"
	Button RangeButton,fSize=11, proc=GCrangeButtonProc
	
	Button RotateButton,pos={126,2.50},size={37,15},title="Rotate",font="Tahoma"
	Button RotateButton,fSize=11, proc=GCRotateButtonProc
	
	TitleBox Axespaneltitle,pos={78,22},size={61.00,12.50},title="Axes Panel"
	TitleBox Axespaneltitle,frame=0, fstyle=1
	
	ypos=33
	TitleBox AxNr,pos={31.00,ypos},size={24.00,19.00},title="Axis",font="Tahoma"
	TitleBox AxNr,frame=0
	TitleBox AxOnOff,pos={53.00,ypos},size={21.00,19.00},title="0/1",font="Tahoma"
	TitleBox AxOnOff,frame=0
	TitleBox AxTics,pos={75.00,ypos},size={28.00,19.00},title="Ticks",font="Tahoma"
	TitleBox AxTics,frame=0
	TitleBox AxTicNr,pos={102.00,ypos},size={28.00,19.00},title="Num",font="Tahoma"
	TitleBox AxTicNr,frame=0
	TitleBox AxLbL,pos={145.00,ypos},size={29.00,19.00},title="Label",font="Tahoma"
	TitleBox AxLbL,frame=0
	
	ypos=46
	DrawLine 29, ypos, 195, ypos
	
	ypos=48
	dy=15
	
	str=igorinfo(2)
	if (cmpstr(str,"Windows")==0)
		dyT=0
		dyC=2
	else
		dyT=1
		dyC=0
	endif
	
	for (ii=0; ii<12; ii+=1)
		ctrl="AXX"+num2str(ii)
		if (ii<4)
			str="X"+num2str(ii)
		elseif ((ii>=4)&&(ii<8))
			str="Y"+num2str(ii-4)
		else
			str="Z"+num2str(ii-8)
		endif
		jj=GCgizmoAxisKey(ii)
		//print str, jj
		TitleBox $ctrl,pos={33.00, ypos+dyT},size={18.00,19.00},title=str,font="Tahoma"
		TitleBox $ctrl,frame=0, fsize=11
		
		ctrl="AX"+num2str(ii)
		CheckBox $ctrl,pos={56.00, ypos+dyC},size={16.00,15.00},title="",font="Tahoma"
		CheckBox $ctrl, value=Ax[jj], proc=GCaxesOnOffCheckProc
		
		ctrl="AXT"+num2str(ii)
		CheckBox $ctrl ,pos={81.00, ypos+dyC},size={16.00,15.00},title="",font="Tahoma"
		CheckBox $ctrl ,value= AxT[jj], proc=GCaxesTickOnOffCheckProc
		
		ctrl="AXN"+num2str(ii)
		CheckBox $ctrl ,pos={106, ypos+dyC},size={16.00,15.00},title="",font="Tahoma"
		CheckBox $ctrl ,value= AxN[jj], proc=GCaxesTickNrOnOffCheckProc
		
		ctrl="AXL"+num2str(ii)
		SetVariable $ctrl ,pos={126.00, ypos+dyT}, size={70.00,16.00},title=" ",  focusRing=0
		SetVariable $ctrl , value= AxL[jj], font="Tahoma", fsize=9.5, proc=GCaxesSetLabelProc
		SetVariable $ctrl  frame=0, labelBack=(65535,65535,65535)  //(55428,55428,55428) 

		ypos+=dy
	endfor
	GCaxesTitles(GizmoName) 
	
	setdatafolder savedDataFolder
end

//------------------------------//
Function GCaxesTitles(GizmoName) 
	string GizmoName
	
	String savedDataFolder = GetDataFolder(1)
	string str, ctrl
	variable ypos, dy, ii, jj
	
	str="root:"+GizmoName
	setdatafolder $str
	wave Ax, AxT, AxL, AxN
	
	for (ii=0; ii<12; ii+=1)
		ctrl="AXX"+num2str(ii)
		jj=GCgizmoAxisKey(ii)

		if (Ax[jj])  
			if (AxN[jj])
				TitleBox $ctrl fstyle=1,fColor=(39321,13101,1)  // color+bold if numerical ticks shown
			else
				TitleBox $ctrl fstyle=1,fColor=(0,0,0) // bold if axis is turned on
			endif 
		else
			TitleBox $ctrl fstyle=0,fColor=(0,0,0) // not turned on
		endif
	endfor

	setdatafolder savedDataFolder
end

//------------------------------//
Function GCgizmoAxisKey(ii)  // strange axis numbering in gizmos, see ModifyGizmo for Axis Objects
	variable ii 
	make/free/N=(12) fwdkey
	fwdkey={0,9,10,11,1,6,7,8,2,3,4,5} // x0,x1,x2,x3, y,y,y,y, z,z,z,z
	return fwdkey[ii]
end

//------------------------------//
Function GCaxesSetLabelProc(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	string str, str2, df, gwin
	variable ii, jj, kk
	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			str=sva.win
			ii=strsearch(str, "_GCpanel",0,2)
			gwin=str[0,ii-1]
			df="root:"+gwin
			str2=df+":AxL"
			wave/T AxL=$str2
			
			str=sva.ctrlName
			ii=strsearch(str, "AXL",0,2)
			jj=str2num(str[ii+3,inf])
			
			kk=GCgizmoAxisKey(jj)
			AxL[kk]=sval
			
			if (strlen(sval))
				ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={kk,axisLabel,1}
			else
				ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={kk,axisLabel,0}
			endif
			ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={kk, axisLabelText, sval}
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCaxesOnOffCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	string str, str2, df, gwin
	variable ii, jj, kk
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			str=cba.win
			ii=strsearch(str, "_GCpanel",0,2)
			gwin=str[0,ii-1]
			df="root:"+gwin
			str2=df+":AX"
			wave AX=$str2
			
			str=cba.ctrlName
			ii=strsearch(str, "AX",0,2)
			jj=str2num(str[ii+2,inf])
			
			kk=GCgizmoAxisKey(jj)
			AX[kk]=checked
			
			ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={kk, visible, checked}
			
			GCaxesTitles(gwin) 
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCaxesTickNrOnOffCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	string str, str2, df, gwin, ctrl
	variable ii, jj, kk
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			str=cba.win
			ii=strsearch(str, "_GCpanel",0,2)
			gwin=str[0,ii-1]
			df="root:"+gwin
			str2=df+":AxN"
			wave AxN=$str2
			
			str2=df+":AxT"
			wave AxT=$str2
			
			str=cba.ctrlName
			ii=strsearch(str, "AXN",0,2)
			jj=str2num(str[ii+3,inf])
			
			kk=GCgizmoAxisKey(jj)
			AxN[kk]=checked
			
			if (checked)
				ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={ kk,ticks,3}
				AxT[kk]=1
				ctrl="AXT"+num2str(jj)
				CheckBox $ctrl ,value= 1
			else
				ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={ kk,ticks,1}
			endif
			
			GCaxesTitles(gwin) 
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCaxesTickOnOffCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	string str, str2, df, gwin, ctrl
	variable ii, jj, kk
	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			str=cba.win
			ii=strsearch(str, "_GCpanel",0,2)
			gwin=str[0,ii-1]
			df="root:"+gwin
			str2=df+":AxT"
			wave AxT=$str2
			
			str2=df+":AxN"
			wave AxN=$str2
			
			str=cba.ctrlName
			ii=strsearch(str, "AXT",0,2)
			jj=str2num(str[ii+3,inf])
			
			kk=GCgizmoAxisKey(jj)
			AxT[kk]=checked
			
			if (checked)
				ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={ kk,ticks,1}
			else
				ModifyGizmo/N=$gwin ModifyObject=axes0,objectType=Axes,property={ kk,ticks,0}
				AxN[kk]=0
				ctrl="AXN"+num2str(jj)
				CheckBox $ctrl ,value= 0
			endif
			
			GCaxesTitles(gwin) 
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCaxesButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			str="KillWindow/Z "+str // qeue panel to be killed
			Execute/P/Q str
			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			str="root:"+gizname+":GCpanelMode"
			nvar GCpanelMode=$str
			GCpanelMode=2 // <0:minimized, 1:default, 2:axes, 3:ranges

			str=""
			str="Chem3Dmodule#GCmakePanel(\""+GizName+"\") " // qeue remake of panel
			Execute/P/Q str
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCaxesCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba
	string str, gizname
	variable ii

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			str=cba.win
				ii=strsearch(cba.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				if (checked)
					ModifyGizmo setDisplayList=2, object=axes0
				else
					RemovefromGizmo/Z displayItem=axes0
				endif
				str="root:"+gizname+":GCaxes"
				nvar GCaxes=$str
				GCaxes=checked 
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//------------------------------//
Function GCpanSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname, wstr
	variable ii
	variable panlim=1.7
	
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				wstr=sa.win
				ii=strsearch(wstr,"_GCpanel",0,2)
				gizname=wstr[0,ii-1]
				str="root:"+gizname+":GCpanX"
				nvar GCpanX=$str
				str="root:"+gizname+":GCpanY"
				nvar GCpanY=$str
				
				controlinfo/W=$wstr HpanSlider
				GCpanX+=v_value
				GCpanX=max(GCpanX, -panlim)
				GCpanX=min(GCpanX,panlim)
				controlinfo/W=$wstr VpanSlider
				GCpanY+=v_value		
				GCpanY=max(GCpanY, -panlim)
				GCpanY=min(GCpanY,panlim)
				ModifyGizmo/N=$gizname pan={GCpanX,GCpanY}
				
				Slider HpanSlider, win=$wstr, value=0
				Slider VpanSlider, win=$wstr, value=0
			endif
			break
	endswitch

	return 0
End


//------------------------------//
Function GCzoomSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname
	variable ii
	
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				str=sa.win
				ii=strsearch(sa.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				ModifyGizmo/N=$gizname zoomFactor=curval
				str="root:"+gizname+":GCzoom"
				nvar GCzoom=$str
				GCzoom=curval
			endif
			break
	endswitch

	return 0
End

//------------------------------//
Function GCMinimizeButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	variable ii
	string str, gizname
	
	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			str="KillWindow/Z "+str // qeue panel to be killed
			Execute/P/Q str
			
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			str="root:"+gizname+":GCpanelMode"
			nvar GCpanelMode=$str
			GCpanelMode*=-1 // flip mode

			str=""
			str="Chem3Dmodule#GCmakePanel(\""+GizName+"\") " // qeue remake of panel
			Execute/P/Q str
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//---------------------------------------------------------------//
Function GCpanelHook(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	string str

	//print "panel hook", s.eventCode
	switch(s.eventCode)
		case 0:	// Activate    
//			str=s.winName
//			dowindow/hide=0 $str  // gizmo hook hides on deactivate, but if
//			str="dowindow/F  "+str  // focus comes here, we must unhide
//			Execute/P/Q str
			break
		case 1:	// Deactivate
			str=s.winName
			str="Chem3Dmodule#GCpanelShowHide(\""+str+"\")"
			execute/P/Q/Z str
			break
		case 12: 	// move
			break
		case 15: 	// hide
			break
		case 16: 	// unhide
			break
	endswitch
	
	return hookResult		// 0 if nothing done, else 1
End

//---------------------------------------------------------------//
Function GCpanelShowHide(pname)
	string pname
	variable ii, pactive, gactive
	string gizname

	GetWindow $pname, active  // panel active?
	pactive=v_value
	
	ii=strsearch(pname,"_GCpanel",0,2)
	gizname=pname[0,ii-1]
	GetWindow/Z $gizname, active   // gizmo active?
	gactive=v_value
	
	if (!pactive && !gactive) // if neither panel nor gizmo are active, hide panel:
		dowindow/hide=1 $pname
	endif
end

//---------------------------------------------------------------//
Function GChook2(s)
	STRUCT WMWinHookStruct &s

	Variable hookResult = 0
	string str

	//print "hook2", s.eventCode
	
	switch(s.eventCode)
		case 0:	// Activate    
			//print "activate"
			// idea is to ensure panel is visilble when gizmo is top window
			// but may cause various problems, such as apparent wrong name in 
			// WIndow Control panel (ctrl-Y)
			str=s.winName+"_GCpanel"
			dowindow/hide=0 $str
			dowindow/F $str
			doupdate
			str=s.winName
			dowindow/F $str
			doupdate
			break
		case 1:	// Deactivate  possibly hide panel
			str=s.winName
			str=s.winName+"_GCpanel" 
			str="Chem3Dmodule#GCpanelShowHide(\""+str+"\")"
			execute/P/Q/Z str
			break
		case 12: 	// move
			str=s.winName
			GCmovePanel(str)
			break
		case 15: 	// hide
			str=s.winName+"_GCpanel"
			dowindow/hide=1 $str
			break
		case 16: 	// unhide
			str=s.winName+"_GCpanel"
			dowindow/hide=0 $str
			str=s.winName
			GCmovePanel(str)
			break
	endswitch
	
	return hookResult		// 0 if nothing done, else 1
End

//---------------------------------------------------------------//
Function GChook(s)
	STRUCT WMGizmoHookStruct &s
	string str, wrec
	variable ii, jj

	//print s.eventName
	strswitch(s.eventName)
		case "mouseDown":
			break
		case "mouseMoved":
			break
		case "rotation":
			break
		case "kill":
			str=s.winName+"_GCpanel"
			str="KillWindow/Z "+str
			execute/P/Q str
			
			str=s.winname
			str="Killdatafolder/Z "+str
			execute/P/Q str
			break
		case "scaling":
			break
		case "scale":  // axes changed
			str=s.winName
			GCrangeAfterChangeProc(str)
			GCresetAxesLimits(str)
			break
		case "wheel":
			break
		case "resize":
			str=s.winName
			GCmovePanel(str)
			str=s.winName+"_GCpanel"
			dowindow/F $str
			break
	endswitch

	return 0 
End

//---------------------------------------------------------------//
Function Chem3DdataChange(xyz, atoms)
	wave xyz
	wave/T atoms
	string gwinlist, winstr, wvstr, fullwvstr, recstr, str, topwin
	variable nwin, ii, jj, kk, mm, r1
	
	fullwvstr=GetWavesDataFolder(xyz,2)
	wvstr=nameofwave(xyz) // xyz wave that has been changed
	
	gwinlist=winlist("*",";","WIN:65536")  // gizmo list, IP7 only
	nwin=itemsinlist(gwinlist)
	for (ii=0; ii<nwin; ii+=1) // check all gizmos
		winstr=stringfromlist(ii, gwinlist)
		if (strsearch(winstr,"Chem3d",0,2)>=0) // a Chem3D gizmo?
			str="root:"+winstr+":BondsOnOff"
			nvar bonds01=$str
			if (bonds01)
				recstr=winrecreation(winstr,0)
				kk=0
				do
					kk=strsearch(recstr,"Scatter=",kk,2) // check all scatter plots
					if (kk>=0)
						kk+=8
						jj=strsearch(recstr,"\r",kk)
						str=recstr[kk, jj]
						jj=strsearch(str,wvstr,0) // our xyz wave?
						if (jj>0) // yes, redraw bonds
							topwin=removeending(winlist("*",";","win:"),";") // save our top win
							dowindow/F winstr // activate relvant gizmo
							Chem3DremoveBonds(xyz)
							String savedDataFolder = GetDataFolder(1)
							str="root:"+winstr
							setdatafolder str
							svar atomlist
							wave/T atoms=$atomlist
							wave coords=$fullwvstr
							Chem3DmakeBonds(atoms, coords)
							wave AtomColors
							Chem3DatomColors(atoms, AtomColors)
							Chem3DatomScale(winstr)	
							setdatafolder savedDataFolder
							dowindow/F $topwin
							break
						endif // is our changed wave
					else
						break
					endif 
				while(1) // scatter plots
			endif // has bonds
		endif // is Chem3D gizmo
	endfor // all gizmos
	return 0
end

//------------------------------//
Function Chem3DBondToggleButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	string str, gizname
	variable ii

	switch( ba.eventCode )
		case 2: // mouse up
			str=ba.win
			ii=strsearch(ba.win,"_GCpanel",0,2)
			gizname=str[0,ii-1]
			dowindow/F gizname // activate relvant gizmo
			String savedDataFolder = GetDataFolder(1)
			str="root:"+gizname
			setdatafolder str
			svar atomlist, xyzlist
			nvar BondsOnOff
			wave/T atoms=$atomlist
			wave xyz=$xyzlist
			if (BondsOnOff)
				Chem3DremoveBonds(xyz)
				BondsOnOff=0
			else
				Chem3DmakeBonds(atoms, xyz)
				BondsOnOff=1
			endif
			setdatafolder savedDataFolder
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

//---------------------------------------------------------------//
Function Chem3DremoveBonds(xyz) // remove bond paths
	wave xyz
	variable ii, jj, npts
	string bondnam
	
	npts=dimsize(xyz,0)
	
	for (ii=0; ii<npts; ii+=1)
		for (jj=ii+1; jj<npts; jj+=1)
			bondnam="bond_"+num2str(ii)+"_"+num2str(jj)
			RemovefromGizmo/Z object=$bondnam
		endfor
	endfor
end

//---------------------------------------------------------------//
Function Chem3DmakeBonds(atoms, xyz) // make bonds as paths
	wave/T atoms
	wave xyz
	variable ii, jj, npts, r1, r2, dd
	string bondnam, atomstr, str
	nvar BondMult
	
	npts=dimsize(xyz,0)
	
	for (ii=0; ii<npts; ii+=1)
		for (jj=ii+1; jj<npts; jj+=1)
			str=atoms[ii]
			atomstr=str[0,0]
			if ((cmpstr(atomstr,"x")==0) && (cmpstr(str,"x")!=0)) // starts with "x", and more than "x"
				atomstr=str[1,inf]
			else
				atomstr=str
			endif
			r1=Chem3DcovRadius(atomstr)
			
			str=atoms[jj]
			atomstr=str[0,0]
			if ((cmpstr(atomstr,"x")==0) && (cmpstr(str,"x")!=0)) // starts with "x", and more than "x"
				atomstr=str[1,inf]
			else
				atomstr=str
			endif
			r2=Chem3DcovRadius(atomstr)
			dd=(xyz[ii][0]-xyz[jj][0])^2
			dd+=(xyz[ii][1]-xyz[jj][1])^2
			dd+=(xyz[ii][2]-xyz[jj][2])^2
			//print ii, jj, (r1+r2), sqrt(dd)
		 
			if (dd < bondtol*(r1+r2)^2) // note tolerance factor
				bondnam="bond_"+num2str(ii)+"_"+num2str(jj)
				make/O/N=(2,3) $bondnam
				wave bond=$bondnam
				bond[0][0]=xyz[ii][0]
				bond[0][1]=xyz[ii][1]
				bond[0][2]=xyz[ii][2]
				bond[1][0]=xyz[jj][0]
				bond[1][1]=xyz[jj][1]
				bond[1][2]=xyz[jj][2]
				AppendToGizmo path=bond,name=$bondnam
				ModifyGizmo ModifyObject=$bondnam,objectType=path,property={ pathColorType,1}
				ModifyGizmo ModifyObject=$bondnam,objectType=path,property={ lineWidthType,0}
				ModifyGizmo ModifyObject=$bondnam,objectType=path,property={ pathColor,0.2,0.2,0.2,1}
				ModifyGizmo ModifyObject=$bondnam,objectType=path,property={ drawTube,1}
				ModifyGizmo ModifyObject=$bondnam,objectType=path,property={ fixedRadius,bondthick*bondmult}
				ModifyGizmo DisplayLastObject, object=$bondnam
			endif
		endfor
	endfor
end

//---------------------------------------------------------------//
Function Chem3DvdWRadius(estr) // returns vsn der Waals radius in Angstrom
// vdW radii from
// M. Mantina; A.C. Chamberlin; R. Valero; C.J. Cramer; D.G. Truhlar (2009).
// "Consistent van der Waals Radii for the Whole Main Group". 
// J. Phys. Chem. A. 113 (19): 5806–12. doi:10.1021/jp8111556

// where not available from above ref, used vdW or calculated radius from 
// https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
// which are based on 
// E. Clementi; D.L.Raimondi; W.P. Reinhardt (1967). 
// "Atomic Screening Constants from SCF Functions. II. Atoms with 37 to 86 Electrons". 
// J. Chem. Phys. 47: 1300. doi:10.1063/1.1712084.

// no data for La, Ce => used Pr value
// no data for Ho, used mean of Dy and Er
// no data for actinides except U, used U value for all others
//
	string estr // symbol of element
	string teststr
	variable rr, ii
	variable nelem=104
	
	make/Free/N=(nelem) vdWradii
	vdWradii[0]= {1.5,1.1,1.4,1.81,1.53,1.92,1.7,1.55,1.52,1.47,1.54,2.27,1.73,1.84,2.1,1.8,1.8,1.75,1.88,2.75,2.31,2.11,1.76,1.71,1.66,1.61,1.56,1.52,1.63,1.4,1.39,1.87,2.11,1.85,1.9,1.83,2.02,3.03,2.49,2.12}
	vdWradii[40]= {2.06,1.98,1.9,1.83,1.78,1.73,1.63,1.72,1.58,1.93,2.17,2.06,2.06,1.98,2.16,3.43,2.68,2.47,2.47,2.47,2.06,2.05,2.38,2.31,2.33,2.25,2.28,2.27,2.26,2.22,2.22,2.17,2.08,2,2.93,1.88,1.85,1.8,1.75}
	vdWradii[79]= {1.66,1.55,1.96,2.02,2.07,1.97,2.02,2.2,3.48,2.83,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86,1.86}
	
	make/Free/T/N=(nelem) ElementSym
	Elementsym[0]= {"X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y"}
	Elementsym[40]= {"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir"}
	Elementsym[78]= {"Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"}

	teststr=estr[0,0]
	if ((cmpstr(teststr,"x")==0) && (cmpstr(estr,"x")!=0) && (cmpstr(estr,"Xe")!=0)) // starts with "x", but more than "x"
		teststr=estr[1,inf]
	else
		teststr=estr
	endif
	
	rr=1.1
	for (ii=0; ii<nelem; ii+=1)
		if (cmpstr(teststr,ElementSym[ii])==0)
			rr=vdWradii[ii]
			break
		endif
	endfor
	
	return rr
end

//---------------------------------------------------------------//
Function Chem3DcovRadius(estr) // returns single bond covalent radius in Angstrom
// from 
// Molecular Single-Bond Covalent Radii for Elements 1–118
// Pekka Pyykkö and Michiko Atsumi
// Chem. Eur. J. 2009, 15, 186 – 197
// DOI: 10.1002/chem.200800987
//
	string estr // symbol of element
	variable rr, ii
	variable nelem=104
	
	make/Free/N=(nelem) CovRadii  
	CovRadii[0]= {0.75,0.32,0.46,1.33,1.02,0.85,0.75,0.71,0.63,0.64,0.67,1.55,1.39,1.26,1.16,1.11,1.03,0.99,0.96,1.96,1.71,1.48,1.36,1.34,1.22,1.19,1.16,1.11,1.1,1.12,1.18,1.24,1.21,1.21,1.16,1.14,1.17,2.1,1.85}
	CovRadii[39]= {1.63,1.54,1.47,1.38,1.28,1.25,1.25,1.2,1.28,1.36,1.42,1.4,1.4,1.36,1.33,1.31,2.32,1.96,1.8,1.63,1.76,1.74,1.73,1.72,1.68,1.69,1.68,1.67,1.66,1.65,1.64,1.7,1.62,1.52,1.46,1.37,1.31,1.29,1.22}
	CovRadii[78]= {1.23,1.24,1.33,1.44,1.44,1.51,1.45,1.47,1.42,2.23,2.01,1.86,1.75,1.69,1.7,1.71,1.72,1.66,1.66,1.68,1.68,1.65,1.67,1.73,1.76,1.61}
	
	make/Free/T/N=(nelem) ElementSym
	Elementsym[0]= {"X","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y"}
	Elementsym[40]= {"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir"}
	Elementsym[78]= {"Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr"}

	for (ii=0; ii<nelem; ii+=1)
		if (cmpstr(estr,ElementSym[ii])==0)
			rr=CovRadii[ii]
			break
		endif
	endfor
	return rr
end

//------------------------------//
Function Chem3DatomScaleSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname
	variable ii, npts, r1
	
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				str=sa.win
				ii=strsearch(sa.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				
				Chem3DatomScale(gizname)

			endif
			break
	endswitch

	return 0
End

//------------------------------//
Function Chem3DatomScale(gname)
	string gname
	
	string str
	variable ii, r1, npts
	
	String savedDataFolder = GetDataFolder(1)
	str="root:"+gname
	setdatafolder str
	
	svar AtomList
	nvar AtomScale, AtomMode
	wave/T atoms=$atomlist
	wave AtomSizes
	nvar GCxmin, GCxmax, GCymin, GCymax, GCzmin, GCzmax
	variable dx, dy, dz

	if (AtomMode==0) // 0=atoms, 1=numbers
		dx=GCxmax-GCxmin
		dy=GCymax-GCymin
		dz=GCzmax-GCzmin
		npts=numpnts(atoms)
		for (ii=0; ii<npts; ii+=1)  
			r1=2*Chem3DvdwRadius(atoms[ii])
			AtomSizes[ii][0]=r1*atomscale /dx
			AtomSizes[ii][1]=r1*atomscale /dy
			AtomSizes[ii][2]=r1*atomscale /dz
		endfor
	else
		AtomSizes=AtomScale
	endif

	setdatafolder savedDataFolder
end

//------------------------------//
Function Chem3DbondThickSliderProc(sa) : SliderControl
	STRUCT WMSliderAction &sa
	string str, gizname, bondnam
	variable ii, jj, npts
	
	switch( sa.eventCode )
		case -1: // control being killed
			break
		default:
			if( sa.eventCode & 1 ) // value set
				Variable curval = sa.curval
				str=sa.win
				ii=strsearch(sa.win,"_GCpanel",0,2)
				gizname=str[0,ii-1]
				str="root:"+gizname+":xyzlist"
				svar xyzlist=$str
				wave xyz=$xyzlist 
				str="root:"+gizname+":bondmult"
				nvar bondmult=$str
	
				npts=dimsize(xyz,0)
				for (ii=0; ii<npts; ii+=1)
					for (jj=ii+1; jj<npts; jj+=1)
						bondnam="bond_"+num2str(ii)+"_"+num2str(jj)
						str="root:"+gizname+":"+bondnam
						if (exists(str))
							ModifyGizmo/N=$gizname ModifyObject=$bondnam,objectType=path,property={ fixedRadius,bondthick*bondmult}
						endif
					endfor
				endfor
			endif
			break
	endswitch

	return 0
End

//------------------------------//
Function Chem3DatomColors(atoms, cwv) // mostly similar to Jmol: en.wikipedia.org/wiki/CPK_coloring
	wave/T atoms
	wave cwv 
	variable ii, npts
	string teststr, str
	
	
	npts=numpnts(atoms)
	cwv=0  // black
	cwv[][3]=1 // fully opaque
	for (ii=0; ii<npts; ii+=1)  
		str=atoms[ii]
		teststr=str[0,0]
		if (cmpstr(teststr,"x")==0) // starts with "x"
			teststr="X"
		else
			teststr=str
		endif
	
		if (cmpstr(teststr,"H")==0)
			cwv[ii][0]=0.5 // medium gray
			cwv[ii][1]=0.5
			cwv[ii][2]=0.5
		elseif (cmpstr(teststr,"C")==0)
			cwv[ii][0]=0.25 // dark gray
			cwv[ii][1]=0.25
			cwv[ii][2]=0.25
		elseif (cmpstr(teststr,"N")==0)
			cwv[ii][2]=1 // blue
		elseif (cmpstr(teststr,"O")==0)
			cwv[ii][0]=1 // red
		elseif ((cmpstr(teststr,"F")==0)||(cmpstr(teststr,"Cl")==0))
			cwv[ii][1]=1 // green
		elseif ((cmpstr(teststr,"Ne")==0)||(cmpstr(teststr,"Ar")==0)||(cmpstr(teststr,"Kr")==0)||(cmpstr(teststr,"Xe")==0))
			cwv[ii][1]=1 // cyan
			cwv[ii][2]=1
		elseif (cmpstr(teststr,"S")==0)
			cwv[ii][0]=1 // yellow
			cwv[ii][1]=1
		elseif (cmpstr(teststr,"X")==0)
			cwv[ii][0]=1 // pink
			cwv[ii][2]=0.8
		endif
	endfor
end

//------------------------------//
Function Chem3Dnotes()
	String nb = "Chem3DNotesAndHelp"
	Dowindow $nb
	if (v_flag)
		Dowindow/F $nb
		return 0
	endif

	NewNotebook/N=$nb/F=1/V=1/K=1/ENCG={2,1}/W=(79.5,58.5,481.5,487.5)
	Notebook $nb defaultTab=36, magnification=125
	Notebook $nb showRuler=0, rulerUnits=2, updating={1, 1}
	Notebook $nb newRuler=Normal, justification=0, margins={0,0,468}, spacing={0,0,0}, tabs={}, rulerDefaults={"Arial",11,0,(0,0,0)}
	Notebook $nb ruler=Normal; Notebook $nb  justification=1, margins={0,0,288}, fStyle=1, text="\r"
	Notebook $nb ruler=Normal; Notebook $nb  tabs={106,302,389}, text="\t", fStyle=5, text="Chem3D\r"
	Notebook $nb ruler=Normal, fStyle=-1, text="\r"
	Notebook $nb ruler=Normal; Notebook $nb  tabs={27}, text="Chem3D is a package for 3D visualization of molecules. \r"
	Notebook $nb text="There are many excellent molecular graphics systems, but\r"
	Notebook $nb text="if you are using Igor for chemical purposes, it can be useful\r"
	Notebook $nb text="to have molecular visualization in Igor.\r"
	Notebook $nb text="\r"
	Notebook $nb text="Topics:\r"
	Notebook $nb text="1. General\r"
	Notebook $nb text="2. Menu and data requirements\r"
	Notebook $nb text="3. Chem3D minimal panel\r"
	Notebook $nb text="4. Chem3D panel\r"
	Notebook $nb text="5. Axes panel\r"
	Notebook $nb text="6. Range panel\r"
	Notebook $nb text="7. Rot/Rotation panel\r"
	Notebook $nb ruler=Normal, text="\r"
	Notebook $nb fStyle=1, text="1. General\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1
	Notebook $nb text="Chem3D uses a scatter plot to display the molecule in an Igor Gizmo window. The control panels are deriv"
	Notebook $nb text="ed from the GizmoControl package. This requires Igor 7.03 or higher, due to several Igor bugs that were "
	Notebook $nb text="found during development.\r"
	Notebook $nb text="\r"
	Notebook $nb text="The package is an Independent Module, to avoid conflict with  other namespaces. Normally it is not visib"
	Notebook $nb text="le, to see it, use:\r"
	Notebook $nb text="SetIgorOption IndependentModuleDev=1\r"
	Notebook $nb text="\r"
	Notebook $nb fStyle=4, text="Chemical bonds\r"
	Notebook $nb fStyle=-1
	Notebook $nb text="Bonds are automatically created, based on single bond covalent radii. See the function Chem3DcovRadius f"
	Notebook $nb text="or details and references. A 15% additional distance tolerance is used, this can be modified at the top "
	Notebook $nb text="of the procedure file: constant bondtol.\r"
	Notebook $nb text="\r"
	Notebook $nb fStyle=4, text="Atomic colors and radii\r"
	Notebook $nb fStyle=-1, text="Atom coloring is loosely based on Jmol. See function Chem3DatomColors.\r"
	Notebook $nb text="\r"
	Notebook $nb text="The atomic sphere sizes can be varied, the maximum size is currently the van der Waals radii. See functi"
	Notebook $nb text="on Chem3DvdWRadius for details and references.\r"
	Notebook $nb text="\r"
	Notebook $nb ruler=Normal, fStyle=4, text="Element X and prefix X\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1
	Notebook $nb text="In addition to all real elements, a fictitious element X is recognized. X atoms are colored bright pink,"
	Notebook $nb text=" with covalent radius 0.75 and van der Waals radius 1.5. This may be useful for emphasizing special atom"
	Notebook $nb text="s in the molecule.\r"
	Notebook $nb text="\r"
	Notebook $nb text="In addition, elements can be prefixed with x or X. For example xC or XAg. In this case the atom is color"
	Notebook $nb text="ed pink, but the radius is as for the unprefixed element. This is probably more useful than element X. \r"
	Notebook $nb text="\r"
	Notebook $nb fStyle=4, text="Changing the molecule\r"
	Notebook $nb fStyle=-1
	Notebook $nb text="Changes in the x,y,z and element waves are automatically recognized, and the Chem3D plot updated. \r"
	Notebook $nb text="\r"
	Notebook $nb ruler=Normal, fStyle=1, text="2. Menu and data requirements\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1, text="Use the ", fStyle=2, text="Make new Chem3D"
	Notebook $nb fStyle=-1
	Notebook $nb text=" menu item to create a new 3D plot. The panel shows waves which can be used. The molecular atomic coordi"
	Notebook $nb text="nates must be in an Nx3 (x,y,z) triplet wave, in Angstrom (0.1 nm). You must also specify an N-element t"
	Notebook $nb text="ext wave with the elemental symbols for each atom.\r"
	Notebook $nb text="\r"
	Notebook $nb text="From the menu you can also open windows allowing control of advanced aspects of the Gizmo plot, see the "
	Notebook $nb text="Igor helps, also via the Chem3D menu.\r"
	Notebook $nb text="\r"
	Notebook $nb ruler=Normal, fStyle=1, text="3. Chem3D Minimal Panel \r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1
	Notebook $nb text="By pressing the < button in the upper left corner, the panel is minimized, showing only a zoom slider an"
	Notebook $nb text="d a checkbox to enable or disable the axes. Pressing > opens the last active panel.\r"
	Notebook $nb text="\r"
	Notebook $nb ruler=Normal, fStyle=1, text="4. Chem3D Panel\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1
	Notebook $nb text="You probably seldom need more than this panel. Bond at atom size panels and buttons should be self-expla"
	Notebook $nb text="natory. \r"
	Notebook $nb text="\r"
	Notebook $nb text="The ", fStyle=2, text="Atom Numbers", fStyle=-1
	Notebook $nb text=" button shows the order of atoms in your lists, to help you identify and edit them.\r"
	Notebook $nb text="\r"
	Notebook $nb text="The ", fStyle=2, text="Scale Axes", fStyle=-1
	Notebook $nb text=" slider isotropically scales all axes, which thereby control clipping. For non-isotropic axes, see the "
	Notebook $nb fStyle=2, text="Ranges", fStyle=-1, text=" panel. See the ", fStyle=2, text="Ranges", fStyle=-1
	Notebook $nb text=" panel for discussion of the spring-loaded slider function. \r"
	Notebook $nb text="\r"
	Notebook $nb text="Using window zoom or zoomed axes may mean that the area of interest is out of the window. The ", fStyle=2
	Notebook $nb text="Move in Window ", fStyle=-1, text="sliders allow you to see your area of interest.\r"
	Notebook $nb text="\r"
	Notebook $nb text="The ", fStyle=2, text="Save 2D", fStyle=-1, text=" and ", fStyle=2, text="Graph 2D", fStyle=-1
	Notebook $nb text=" buttons make snapshot images of the current Chem3D window. The former saves to file, the latter makes a"
	Notebook $nb text=" new 2D image graph, which can be modified and saved. \r"
	Notebook $nb text="\r"
	Notebook $nb ruler=Normal, fStyle=1, text="5. Axes Panel\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1
	Notebook $nb text="Your preferred rotation of the molecule may cause the axes to be in the way, or you may want different l"
	Notebook $nb text="abels, ticks, etc. These can be easily controlled in this panel. See the Chem3D Advanced menu for more a"
	Notebook $nb text="xis control.\r"
	Notebook $nb text="\r"
	Notebook $nb text="The axis name (X0, Y3, etc) is bold if it is displayed (\"0/1\"), and colored if numerical tick labels (\"n"
	Notebook $nb text="um\")  are shown.\r"
	Notebook $nb text="\r"
	Notebook $nb text="If numerical labels are selected, ticks are automatically also selected. If ticks are deselected, numeri"
	Notebook $nb text="cal labels are also automatically turned off. \r"
	Notebook $nb text="\r"
	Notebook $nb ruler=Normal, fStyle=1, text="6. Range Panel\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}, fStyle=-1
	Notebook $nb text="In Chem3D you should rarely need the Range panel, because\r"
	Notebook $nb text="non-isotropic axes distort the molecule. Zooming and uniform scaling are available in the Chem3D panel. "
	Notebook $nb text="\r"
	Notebook $nb text="\r"
	Notebook $nb text="However, if you do set non-uniform axes, you can return to the standard view using the ", fStyle=2
	Notebook $nb text="AllData", fStyle=-1, text=" then ", fStyle=2, text="Cube", fStyle=-1
	Notebook $nb text=" buttons. The first ensures that the axes encompass the whole molecule, including atom spheres. The seco"
	Notebook $nb text="nd makes the axes isotropic.\r"
	Notebook $nb text="\r"
	Notebook $nb fStyle=4, text="Sliders\r"
	Notebook $nb fStyle=-1
	Notebook $nb text="The range sliders are proportional rather than absolute. A full-scale deflection changes the range by a "
	Notebook $nb text="factor of 2. Bigger to the right, smaller to the left.\r"
	Notebook $nb ruler=Normal, text="\r"
	Notebook $nb ruler=Normal; Notebook $nb  margins={0,0,287}
	Notebook $nb text="The sliders are also \"spring-loaded\", and return to the center position after each setting. \r"
	Notebook $nb text="\r"
	Notebook $nb text="In addition to \"pulling\" the slider from the center, you can also simply click on the slider scale.\r"
	Notebook $nb text="\r"
	Notebook $nb text="Using <shift> while moving or clicking the slider changes the axis offset rather than the (max-min) diff"
	Notebook $nb text="erence.\r"
	Notebook $nb text="\r"
	Notebook $nb fStyle=4, text="Appearance and clipping\r"
	Notebook $nb fStyle=-1
	Notebook $nb text="The effect of the range settings is very different depending on the gizmo aspect ratio setting. This is "
	Notebook $nb text="controlled by a button in the gizmo tools (4th from the top). See also the helps, by entering the follow"
	Notebook $nb text="ing on the command line:\r"
	Notebook $nb text="\r"
	Notebook $nb text="Displayhelptopic \"Gizmo Display Window Tool Palette\"\r"
	Notebook $nb text="\r"
	Notebook $nb text="In \"aspect ratio mode\", the displayed 3D objects are clipped depending on the ranges selected. If aspect"
	Notebook $nb text=" ratio mode is off, the objects alway fill a box with equal sides, only the size of the axis object is m"
	Notebook $nb text="odified. \r"
	Notebook $nb text="\r"
	Notebook $nb fStyle=1, text="7. Rot / Rotation Panel\r"
	Notebook $nb fStyle=-1, text="The 3D ", fStyle=4, text="view", fStyle=-1
	Notebook $nb text=" of the molecule can be rotated, as in any Gizmo window, but this does not change the molecular "
	Notebook $nb fStyle=4, text="coordinates", fStyle=-1
	Notebook $nb text=" in the XYZ wave. It can be useful to rotate the coordinates for further processing, this is possible in"
	Notebook $nb text=" the rotation panel.\r"
	Notebook $nb text="\r"
	Notebook $nb text="The sliders rotate about each axis, the limits are + and - pi. The X, Y and Z rotations are applied sequ"
	Notebook $nb text="entially to the current ", fStyle=4, text="reference coordinates", fStyle=-1
	Notebook $nb text=". The current rotation can be undone by using the Restore Reference button (or moving the sliders back t"
	Notebook $nb text="o the center position).\r"
	Notebook $nb text="\r"
	Notebook $nb text="The initial reference orientation is the X,Y,Z when the Chem3D was created. The reference can be set to "
	Notebook $nb text="the current orientation using the Save as Reference button. The rotation sliders are then reset to the n"
	Notebook $nb text="eutral position. The old reference orientation cannot be restored, there is no undo buffer."

	// move selection to the start of the notebook and display the selection
	Notebook $nb selection={startOfFile,startOfFile}, findText={"",1}
End

//---------------------------------------------------------------//
#if Exists("PanelResolution") != 3
Static Function PanelResolution(wName)	// For compatibility with Igor7
	String wName
	return 72
End
#endif

// I am happy if you use and modify Chem3D, but please do not sell my work.

// License: 
// Creative Commons Attribution-NonCommercial-ShareAlike 4.0  International License
// http://creativecommons.org/licenses/by-nc-sa/4.0/
//
// This software is provided by the copyright holder "as is", with no express or implied warranties. 
