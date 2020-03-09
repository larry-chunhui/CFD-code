# Microsoft Developer Studio Project File - Name="les" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=les - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "les.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "les.mak" CFG="les - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "les - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "les - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "les - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "les - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ  /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "les - Win32 Release"
# Name "les - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\bc__get.c
# End Source File
# Begin Source File

SOURCE=.\convective.c
# End Source File
# Begin Source File

SOURCE=.\convective_adv__2.c
# End Source File
# Begin Source File

SOURCE=.\convective_adv__4.c
# End Source File
# Begin Source File

SOURCE=.\convective_div__2.c
# End Source File
# Begin Source File

SOURCE=.\convective_div__4.c
# End Source File
# Begin Source File

SOURCE=.\correction.c
# End Source File
# Begin Source File

SOURCE=.\derivatives__2.c
# End Source File
# Begin Source File

SOURCE=.\derivatives__4.c
# End Source File
# Begin Source File

SOURCE=.\differentiate.c
# End Source File
# Begin Source File

SOURCE=.\error_check.c
# End Source File
# Begin Source File

SOURCE=.\fasttrig.c
# End Source File
# Begin Source File

SOURCE=.\fft.c
# End Source File
# Begin Source File

SOURCE=.\interpolation.c
# End Source File
# Begin Source File

SOURCE=.\interpolation__2.c
# End Source File
# Begin Source File

SOURCE=.\interpolation__4.c
# End Source File
# Begin Source File

SOURCE=".\interpolation__les copy.c"
# End Source File
# Begin Source File

SOURCE=.\interpolation__les.c
# End Source File
# Begin Source File

SOURCE=.\les.c
# End Source File
# Begin Source File

SOURCE=.\les_diagnostics.c
# End Source File
# Begin Source File

SOURCE=.\main.c
# End Source File
# Begin Source File

SOURCE=.\near_wall_ode.c
# End Source File
# Begin Source File

SOURCE=.\nrutil.c
# End Source File
# Begin Source File

SOURCE=.\parameters.c
# End Source File
# Begin Source File

SOURCE=.\pentasolver.c
# End Source File
# Begin Source File

SOURCE=.\poisson_solver_transpose.c
# End Source File
# Begin Source File

SOURCE=.\profiles.c
# End Source File
# Begin Source File

SOURCE=.\rhs.c
# End Source File
# Begin Source File

SOURCE=.\share.c
# End Source File
# Begin Source File

SOURCE=.\spiral.c
# End Source File
# Begin Source File

SOURCE=.\statistics.c
# End Source File
# Begin Source File

SOURCE=.\time_march.c
# End Source File
# Begin Source File

SOURCE=.\trisolver.c
# End Source File
# Begin Source File

SOURCE=.\velgrad_tensor.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\bc__get.h
# End Source File
# Begin Source File

SOURCE=.\convective.h
# End Source File
# Begin Source File

SOURCE=.\convective_adv__2.h
# End Source File
# Begin Source File

SOURCE=.\convective_adv__4.h
# End Source File
# Begin Source File

SOURCE=.\convective_div__2.h
# End Source File
# Begin Source File

SOURCE=.\convective_div__4.h
# End Source File
# Begin Source File

SOURCE=.\correction.h
# End Source File
# Begin Source File

SOURCE=.\definitions.h
# End Source File
# Begin Source File

SOURCE=.\derivatives__2.h
# End Source File
# Begin Source File

SOURCE=.\derivatives__4.h
# End Source File
# Begin Source File

SOURCE=.\differentiate.h
# End Source File
# Begin Source File

SOURCE=.\error_check.h
# End Source File
# Begin Source File

SOURCE=.\fasttrig.h
# End Source File
# Begin Source File

SOURCE=.\fft.h
# End Source File
# Begin Source File

SOURCE=.\float.h
# End Source File
# Begin Source File

SOURCE=.\interpolation.h
# End Source File
# Begin Source File

SOURCE=.\interpolation__2.h
# End Source File
# Begin Source File

SOURCE=.\interpolation__4.h
# End Source File
# Begin Source File

SOURCE=.\interpolation__les.h
# End Source File
# Begin Source File

SOURCE=.\les.h
# End Source File
# Begin Source File

SOURCE=.\les_diagnostics.h
# End Source File
# Begin Source File

SOURCE=.\near_wall_ode.h
# End Source File
# Begin Source File

SOURCE=.\nrutil.h
# End Source File
# Begin Source File

SOURCE=.\parameters.h
# End Source File
# Begin Source File

SOURCE=.\pentasolver.h
# End Source File
# Begin Source File

SOURCE=.\poisson_solver_transpose.h
# End Source File
# Begin Source File

SOURCE=.\profiles.h
# End Source File
# Begin Source File

SOURCE=.\rhs.h
# End Source File
# Begin Source File

SOURCE=.\share.h
# End Source File
# Begin Source File

SOURCE=.\spiral.h
# End Source File
# Begin Source File

SOURCE=.\statistics.h
# End Source File
# Begin Source File

SOURCE=.\time_march.h
# End Source File
# Begin Source File

SOURCE=.\trisolver.h
# End Source File
# Begin Source File

SOURCE=.\velgrad_tensor.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
