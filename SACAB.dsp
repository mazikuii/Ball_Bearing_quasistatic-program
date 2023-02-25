# Microsoft Developer Studio Project File - Name="SACAB" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=SACAB - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "SACAB.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "SACAB.mak" CFG="SACAB - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "SACAB - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "SACAB - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "SACAB - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
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
# ADD LINK32 dfor.lib kernel32.lib /nologo /subsystem:console /incremental:yes /machine:I386
# SUBTRACT LINK32 /debug /nodefaultlib

!ELSEIF  "$(CFG)" == "SACAB - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 dfor.lib kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /incremental:no

!ENDIF 

# Begin Target

# Name "SACAB - Win32 Release"
# Name "SACAB - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Group "SRC"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Group "HYBRD"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\HYBRD\dogleg.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\dpmpar.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\enorm.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\fdjac1.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\hybrd.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\hybrd1.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\qform.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\qrfac.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\r1mpyq.f
# End Source File
# Begin Source File

SOURCE=.\HYBRD\r1updt.f
# End Source File
# End Group
# Begin Group "BROYDN"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\BROYDN\BROYDN.F90
DEP_F90_BROYD=\
	".\Release\FMINLN.mod"\
	".\Release\fsbr.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\FDJAC.F90
DEP_F90_FDJAC=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\FMIN.F90
DEP_F90_FMIN_=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\LNSRCH.F90
DEP_F90_LNSRC=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\LUBKSB.F90
DEP_F90_LUBKS=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\LUDCMP.F90
DEP_F90_LUDCM=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\POLIN2.F90
DEP_F90_POLIN=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\POLINT.F90
DEP_F90_POLINT=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\PYTHAG.F90
DEP_F90_PYTHA=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\QRDCMP.F90
DEP_F90_QRDCM=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\QROMB.F90
DEP_F90_QROMB=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\QRUPDT.F90
DEP_F90_QRUPD=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\ROTATE.F90
DEP_F90_ROTAT=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\RSOLV.F90
DEP_F90_RSOLV=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\RTSAFE.F90
DEP_F90_RTSAF=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\BROYDN\TRAPZD.F90
DEP_F90_TRAPZ=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# End Group
# Begin Group "NR"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\nr.f90
DEP_F90_NR_F9=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\nrtype_alter.f90
# End Source File
# Begin Source File

SOURCE=.\nrutil_alter.f90
DEP_F90_NRUTI=\
	".\Release\nrtype.mod"\
	
# End Source File
# End Group
# Begin Group "NL2SOL"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\NL2SOL.FOR
DEP_F90_NL2SO=\
	".\Release\nrtype.mod"\
	
# End Source File
# End Group
# Begin Group "ALCON1"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\ALCON\alcon1.f
# End Source File
# Begin Source File

SOURCE=.\ALCON\linalg_alcon1.f
# End Source File
# Begin Source File

SOURCE=.\ALCON\zibconst.f
# End Source File
# End Group
# Begin Group "ALCON2"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\ALCON\alcon2.f
# End Source File
# Begin Source File

SOURCE=.\ALCON\linpack_alcon2.f
# End Source File
# End Group
# Begin Group "NITSOL"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\NITSOL\nitbd.f
DEP_F90_NITBD=\
	".\NITSOL\nitdflts.h"\
	".\NITSOL\nitparam.h"\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitbt.f
DEP_F90_NITBT=\
	".\NITSOL\nitparam.h"\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitdflts.h
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitdrv.f
DEP_F90_NITDR=\
	".\NITSOL\nitinfo.h"\
	".\NITSOL\nitparam.h"\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitgm.f
DEP_F90_NITGM=\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitinfo.h
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitjv.f
DEP_F90_NITJV=\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitparam.h
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitsol.f
DEP_F90_NITSO=\
	".\NITSOL\nitdflts.h"\
	".\NITSOL\nitinfo.h"\
	".\NITSOL\nitparam.h"\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nitstb.f
DEP_F90_NITST=\
	".\NITSOL\nitprint.h"\
	
# End Source File
# Begin Source File

SOURCE=.\NITSOL\nittfq.f
DEP_F90_NITTF=\
	".\NITSOL\nitprint.h"\
	
# End Source File
# End Group
# Begin Group "LAPACK"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\lapack\dasum.f
# End Source File
# Begin Source File

SOURCE=.\lapack\daxpy.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dcopy.f
# End Source File
# Begin Source File

SOURCE=.\lapack\ddot.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dgemv.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dgeqpf.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dgeqr2.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dgeqrf.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlabad.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlacn2.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlaic1.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlamch.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlantr.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlapy2.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlarf.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlarfb.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlarfp.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlarft.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlartg.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlassq.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dlatrs.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dnrm2.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dorg2r.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dorgqr.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dorm2r.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dormqr.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dormr2.f
# End Source File
# Begin Source File

SOURCE=.\lapack\drscl.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dscal.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dswap.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dtpmv.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dtpsv.f
# End Source File
# Begin Source File

SOURCE=.\lapack\dtrcon.f
# End Source File
# Begin Source File

SOURCE=.\lapack\idamax.f
# End Source File
# Begin Source File

SOURCE=.\lapack\ieeeck.f
# End Source File
# Begin Source File

SOURCE=.\lapack\iladlc.f
# End Source File
# Begin Source File

SOURCE=.\lapack\iladlr.f
# End Source File
# Begin Source File

SOURCE=.\lapack\ilaenv.f
# End Source File
# Begin Source File

SOURCE=.\lapack\iparmq.f
# End Source File
# Begin Source File

SOURCE=.\lapack\lsame.f
# End Source File
# End Group
# Begin Group "HOMOPACK"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\HOMPACK\HOMPACK90.f
# End Source File
# End Group
# Begin Group "LBFGS"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\Lbfgsb.f
# End Source File
# End Group
# Begin Group "BFGS"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\BFGS\dfpmin.f90
DEP_F90_DFPMI=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# End Group
# Begin Group "FRPRMN"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\FRPRMN\brent.f90
DEP_F90_BRENT=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FRPRMN\frprmn.f90
DEP_F90_FRPRM=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FRPRMN\linmin.f90
DEP_F90_LINMI=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FRPRMN\mnbrak.f90
DEP_F90_MNBRA=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# End Group
# Begin Group "UNCMIN"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\UNCMIN\bakslv_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\chlhsn_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\choldc_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\d1fcn_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\d2fcn_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\dfault_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\dogdrv_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\dogstp_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\explain_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\forslv_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\fstocd_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\fstofd_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\grdchk_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\heschk_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\hookdr_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\hookst_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\hsnint_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\lltslv_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\lnsrch_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\mvmltl_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\mvmlts_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\mvmltu_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\optchk_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\optdrv_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\optif0_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\optif9_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\optstp_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\qraux1_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\qraux2_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\qrupdt_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\r8_epsilon_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\r8vec_dot_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\r8vec_norm_l2_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\result_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\sclmul_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\secfac_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\secunf_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\sndofd_u.for
# End Source File
# Begin Source File

SOURCE=.\UNCMIN\tregup_u.for
# End Source File
# End Group
# Begin Group "NLEQ1"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\AFFINE\NLEQ1\linalg_nleq1.f
# End Source File
# Begin Source File

SOURCE=.\AFFINE\NLEQ1\nleq1.f
# End Source File
# Begin Source File

SOURCE=.\AFFINE\NLEQ1\wnorm.f
# End Source File
# Begin Source File

SOURCE=.\AFFINE\NLEQ1\zibmon.f
# End Source File
# Begin Source File

SOURCE=.\AFFINE\NLEQ1\zibsec.f
# End Source File
# End Group
# Begin Group "NLEQ2"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\AFFINE\NLEQ2\linalg_nleq2.f
# End Source File
# Begin Source File

SOURCE=.\AFFINE\NLEQ2\nleq2.f
# End Source File
# End Group
# Begin Group "GIANT_GBIT"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\GIANT_GBIT\gbit.c
# End Source File
# Begin Source File

SOURCE=.\GIANT_GBIT\giant.h
# End Source File
# Begin Source File

SOURCE=.\GIANT_GBIT\giant_gbit.c
# End Source File
# Begin Source File

SOURCE=.\GIANT_GBIT\itlin.h
# End Source File
# Begin Source File

SOURCE=.\GIANT_GBIT\test_giantgbit.c
# End Source File
# Begin Source File

SOURCE=.\GIANT_GBIT\utils_giant.c
# End Source File
# Begin Source File

SOURCE=.\GIANT_GBIT\utils_itlin.c
# End Source File
# End Group
# Begin Group "NESLV"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\NESLV\ddoglg.f
# End Source File
# Begin Source File

SOURCE=.\NESLV\DENORM.f
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEcndjac.F
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEcompmu.F
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEDDLG.F
DEP_F90_NEDDL=\
	".\Release\neintface.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEDFLT.F90
DEP_F90_NEDFL=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEdogleg.f90
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEDOGSTP.F
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEEXPLN.F90
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEFDJAC.F90
DEP_F90_NEFDJ=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEFDJAC4.F90
DEP_F90_NEFDJA=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEFMIN.F90
DEP_F90_NEFMI=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEFMQR.F90
DEP_F90_NEFMQ=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEINCHK.F90
DEP_F90_NEINC=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NELIQREV.F
# End Source File
# Begin Source File

SOURCE=.\NESLV\NELNSH.F90
DEP_F90_NELNS=\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEMOD.F90
DEP_F90_NEMOD=\
	".\Release\fsbr.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NEQRUDT.F90
DEP_F90_NEQRU=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NESLV.F90
DEP_F90_NESLV=\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NESLV4.F90
DEP_F90_NESLV4=\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NESLVSN.F90
DEP_F90_NESLVS=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NESTOP.F90
DEP_F90_NESTO=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NESTOP0.F90
DEP_F90_NESTOP=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NESLV\NETRUP.F
# End Source File
# Begin Source File

SOURCE=.\NESLV\TOOLS.F
# End Source File
# End Group
# Begin Group "SOLVERS"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=".\SLV01-BROYDN.F90"
DEP_F90_SLV01=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV02-LBFGSB.F90"
DEP_F90_SLV02=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV03-TENSOLVE.F90"
DEP_F90_SLV03=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV04-SA.FOR"
DEP_F90_SLV04=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV05-GCNEWTON.F90"
DEP_F90_SLV05=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV06-NEWTON-RAPHSON.F90"
DEP_F90_SLV06=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV07-HYBRD.F90"
DEP_F90_SLV07=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV08-NL2SOL.F90"
DEP_F90_SLV08=\
	".\Release\fsbr.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV09-SMSNO.F90"
DEP_F90_SLV09=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV10-DNEQNF.F90"
DEP_F90_SLV10=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	{$(INCLUDE)}"IMSLF90.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV11-DNEQNJ.F90"
DEP_F90_SLV11=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV12-DNEQBF.F90"
DEP_F90_SLV12=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV13-DNEQBJ.F90"
DEP_F90_SLV13=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV14-ALCON2.F90"
DEP_F90_SLV14=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV15-ALCON1.F90"
DEP_F90_SLV15=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV16-HOMPACK.F90"
DEP_F90_SLV16=\
	".\Release\fsbr.mod"\
	".\Release\HOMPACK90.mod"\
	".\Release\HOMPACK90_GLOBAL.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	".\Release\REAL_PRECISION.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV17-BFGS.F90"
DEP_F90_SLV17=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV18-FRPR.F90"
DEP_F90_SLV18=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV19-UNCMIN.FOR"
DEP_F90_SLV19=\
	".\Release\intface.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV20-NITSOL.F90"
DEP_F90_SLV20=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV21-NLEQ1.f"
DEP_F90_SLV21=\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV22-NLEQ2.f"
# End Source File
# Begin Source File

SOURCE=".\SLV23-GIANT_GBIT.F90"
# End Source File
# Begin Source File

SOURCE=".\SLV24-NESLV.F90"
DEP_F90_SLV24=\
	".\Release\fsbr.mod"\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\TT.F90
DEP_F90_TT_F9=\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# End Group
# Begin Source File

SOURCE=.\mnewt.f90
DEP_F90_MNEWT=\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\newt.f90
DEP_F90_NEWT_=\
	".\Release\FMINLN.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SMSNO.FOR
# End Source File
# Begin Source File

SOURCE=.\TENSOLVE.f
# End Source File
# End Group
# Begin Source File

SOURCE=.\FASTSIM.F90
DEP_F90_FASTS=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\FUNCS.F90
DEP_F90_FUNCS=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GETABK.F90
DEP_F90_GETAB=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GETCIJ.F90
DEP_F90_GETCI=\
	".\Release\fsbr.mod"\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GETCREEP.F90
DEP_F90_GETCR=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GETPCHAXIS.F90
DEP_F90_GETPC=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\GETSTATIC.F90
DEP_F90_GETST=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\global.f90
DEP_F90_GLOBA=\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\HYBRD_asym.F90
DEP_F90_HYBRD=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\INITBR.F90
DEP_F90_INITB=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\INITFS.F90
DEP_F90_INITF=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\INITS.F90
DEP_F90_INITS=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\INTFACE.F90
DEP_F90_INTFA=\
	".\Release\fsbr.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\NOTES.f90
# End Source File
# Begin Source File

SOURCE=.\OUTPUTS.f90
DEP_F90_OUTPU=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\PREPOST.F90
DEP_F90_PREPO=\
	".\Release\fsbr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\ratint.f90
DEP_F90_RATIN=\
	".\Release\nrtype.mod"\
	".\Release\nrutil.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SACAB.f90
DEP_F90_SACAB=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV2EQS-ASYM.F90"
DEP_F90_SLV2E=\
	".\Release\fsbr.mod"\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV2EQS-SYM.F90"
DEP_F90_SLV2EQ=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV3EQS-CREEP.F90"
DEP_F90_SLV3E=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV3EQS-INERT.F90"
DEP_F90_SLV3EQ=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV3EQS-ST.F90"
DEP_F90_SLV3EQS=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV6EQS-CREEP.F90"
DEP_F90_SLV6E=\
	".\Release\fsbr.mod"\
	".\Release\neintface.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=".\SLV6EQS-IDX.f90"
DEP_F90_SLV6EQ=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SLV6EQS.f90
DEP_F90_SLV6EQS=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SLV8EQS.f90
DEP_F90_SLV8E=\
	".\Release\fsbr.mod"\
	".\Release\nrtype.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\SLVCMBD.F90
DEP_F90_SLVCM=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\neintface.mod"\
	".\Release\nr.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\TEST.f90
DEP_F90_TEST_=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\WRTTM.F90
DEP_F90_WRTTM=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\XFVCT.f90
DEP_F90_XFVCT=\
	".\Release\fsbr.mod"\
	".\Release\intface.mod"\
	".\Release\nr.mod"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
