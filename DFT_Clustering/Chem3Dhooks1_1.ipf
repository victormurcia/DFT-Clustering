#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//---------------------------------------------------------------//
Function AfterCompiledHook( ) // to reset Chem3D windows after expt reloaded
// this is always called after loading expt, but afteropenhook is not, depending on how expt was selected.
	Chem3Dmodule#GCcheckAllRefresh()
	return 0
end
//---------------------------------------------------------------//
Function BeforeExperimentSaveHook(rN,fileName,path,type,creator,kind)  // marks Chem3D windows as needing refresh
	Variable rN,kind
	String fileName,path,type,creator
	Chem3Dmodule#GCmarkAllRefresh() 
End
