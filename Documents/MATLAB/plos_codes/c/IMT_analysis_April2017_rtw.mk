###########################################################################
## Makefile generated for MATLAB file/project 'IMT_analysis_April2017'. 
## 
## Makefile     : IMT_analysis_April2017_rtw.mk
## Generated on : Fri Jul 21 19:17:47 2017
## MATLAB Coder version: 3.3 (R2017a)
## 
## Build Info:
## 
## Final product: $(RELATIVE_PATH_TO_ANCHOR)\IMT_analysis_April2017.lib
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPUTER                Computer type. See the MATLAB "computer" command.
# PERL                    PERL Tool
# GEN_LNK_SCRIPT          Perl script to generate the command file
# CMD_FILE                Command file

PRODUCT_NAME              = IMT_analysis_April2017
MAKEFILE                  = IMT_analysis_April2017_rtw.mk
COMPUTER                  = PCWIN64
MATLAB_ROOT               = C:\PROGRA~1\MATLAB\R2017a
MATLAB_BIN                = C:\PROGRA~1\MATLAB\R2017a\bin
MATLAB_ARCH_BIN           = C:\PROGRA~1\MATLAB\R2017a\bin\win64
MASTER_ANCHOR_DIR         = 
START_DIR                 = C:\Users\remof\Documents\MATLAB\PlosCodes\Documents\MATLAB\plos_codes
ARCH                      = win64
RELATIVE_PATH_TO_ANCHOR   = .
PERL                      = $(MATLAB_ROOT)\sys\perl\win32\bin\perl.exe
GEN_LNK_SCRIPT            = $(MATLAB_ROOT)\rtw\c\tools\mkvc_lnk.pl
CMD_FILE                  = $(PRODUCT_NAME).lnk
C_STANDARD_OPTS           = 
CPP_STANDARD_OPTS         = 
NODEBUG                   = 1

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          Microsoft Visual C++ 2015 v14.0 | nmake (64-bit Windows)
# Supported Version(s):    14.0
# ToolchainInfo Version:   R2017a
# Specification Revision:  1.0
# 
#-------------------------------------------
# Macros assumed to be defined elsewhere
#-------------------------------------------

# NODEBUG
# cvarsdll
# cvarsmt
# conlibsmt
# ldebug
# conflags
# cflags

#-----------
# MACROS
#-----------

MEX_OPTS_FILE       = $(MATLAB_ROOT)\bin\$(ARCH)\mexopts\msvc2015.xml
MW_EXTERNLIB_DIR    = $(MATLAB_ROOT)\extern\lib\win64\microsoft
MW_LIB_DIR          = $(MATLAB_ROOT)\lib\win64
MEX_ARCH            = -win64
CPU                 = AMD64
APPVER              = 5.02
CVARSFLAG           = $(cvarsmt)
CFLAGS_ADDITIONAL   = -D_CRT_SECURE_NO_WARNINGS
CPPFLAGS_ADDITIONAL = -EHs -D_CRT_SECURE_NO_WARNINGS
LIBS_TOOLCHAIN      = $(conlibs)

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: Microsoft Visual C Compiler
CC = cl

# Linker: Microsoft Visual C Linker
LD = link

# C++ Compiler: Microsoft Visual C++ Compiler
CPP = cl

# C++ Linker: Microsoft Visual C++ Linker
CPP_LD = link

# Archiver: Microsoft Visual C/C++ Archiver
AR = lib

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_BIN)
MEX = "$(MEX_PATH)\mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: NMAKE Utility
MAKE = nmake


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -Zi
C_OUTPUT_FLAG       = -Fo
LDDEBUG             = /DEBUG
OUTPUT_FLAG         = -out:
CPPDEBUG            = -Zi
CPP_OUTPUT_FLAG     = -Fo
CPPLDDEBUG          = /DEBUG
OUTPUT_FLAG         = -out:
ARDEBUG             =
STATICLIB_OUTPUT_FLAG = -out:
MEX_DEBUG           = -g
RM                  = @del
ECHO                = @echo
MV                  = @ren
RUN                 = @cmd /C

#----------------------------------------
# "Faster Builds" Build Configuration
#----------------------------------------

ARFLAGS              = /nologo
CFLAGS               = $(cflags) $(CVARSFLAG) $(CFLAGS_ADDITIONAL) \
                       /Od /Oy-
CPPFLAGS             = /TP $(cflags) $(CVARSFLAG) $(CPPFLAGS_ADDITIONAL) \
                       /Od /Oy-
CPP_LDFLAGS          = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN)
CPP_SHAREDLIB_LDFLAGS  = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN) \
                         -dll -def:$(DEF_FILE)
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN)
MEX_CPPFLAGS         =
MEX_CPPLDFLAGS       =
MEX_CFLAGS           = $(MEX_ARCH) OPTIMFLAGS="/Od /Oy- $(MDFLAG) $(DEFINES)" $(MEX_OPTS_FLAG)
MEX_LDFLAGS          = LDFLAGS=='$$LDFLAGS'
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = $(ldebug) $(conflags) $(LIBS_TOOLCHAIN) \
                       -dll -def:$(DEF_FILE)

#--------------------
# File extensions
#--------------------

H_EXT               = .h
OBJ_EXT             = .obj
C_EXT               = .c
EXE_EXT             = .exe
SHAREDLIB_EXT       = .dll
HPP_EXT             = .hpp
OBJ_EXT             = .obj
CPP_EXT             = .cpp
EXE_EXT             = .exe
SHAREDLIB_EXT       = .dll
STATICLIB_EXT       = .lib
MEX_EXT             = .mexw64
MAKE_EXT            = .mk


###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = $(RELATIVE_PATH_TO_ANCHOR)\IMT_analysis_April2017.lib
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = $(START_DIR)\codegen\lib\IMT_analysis_April2017;$(START_DIR);$(MATLAB_ROOT)\extern\include;$(MATLAB_ROOT)\simulink\include;$(MATLAB_ROOT)\rtw\c\src;$(MATLAB_ROOT)\rtw\c\src\ext_mode\common;$(MATLAB_ROOT)\rtw\c\ert

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_STANDARD = -DMODEL=IMT_analysis_April2017 -DHAVESTDIO -DUSE_RTMODEL

DEFINES = $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)\codegen\lib\IMT_analysis_April2017\IMT_analysis_April2017_rtwutil.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\IMT_analysis_April2017_data.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\IMT_analysis_April2017_initialize.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\IMT_analysis_April2017_terminate.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\IMT_analysis_April2017.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\fprintf.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\fileManager.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\strcmp.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\rdivide.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\power.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\sum.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\sqrt.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\fminsearch.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\emgpdf.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\abs.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\exp.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\log.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\sortIdx.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\eps.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\onestagepdf2.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\onestagepdf_lag.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\convolv_2invG_adapt_nov.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\sort1.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\tailmass.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\gp_max.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\roots.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\anyNonFinite.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xgehrd.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xzlarfg.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xnrm2.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xdlapy3.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xscal.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\recip.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xzhseqr.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xzgeev.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xzlartg.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\xzhgeqz.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\onestagepdf_prime.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\approxconvolv.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\conv.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\convolv_3invG_nov.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\convolv_2invG_nov.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\IMT_analysis_April2017_emxutil.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\rt_nonfinite.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\rtGetNaN.c $(START_DIR)\codegen\lib\IMT_analysis_April2017\rtGetInf.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = IMT_analysis_April2017_rtwutil.obj IMT_analysis_April2017_data.obj IMT_analysis_April2017_initialize.obj IMT_analysis_April2017_terminate.obj IMT_analysis_April2017.obj fprintf.obj fileManager.obj strcmp.obj rdivide.obj power.obj sum.obj sqrt.obj fminsearch.obj emgpdf.obj abs.obj exp.obj log.obj sortIdx.obj eps.obj onestagepdf2.obj onestagepdf_lag.obj convolv_2invG_adapt_nov.obj sort1.obj tailmass.obj gp_max.obj roots.obj anyNonFinite.obj xgehrd.obj xzlarfg.obj xnrm2.obj xdlapy3.obj xscal.obj recip.obj xzhseqr.obj xzgeev.obj xzlartg.obj xzhgeqz.obj onestagepdf_prime.obj approxconvolv.obj conv.obj convolv_3invG_nov.obj convolv_2invG_nov.obj IMT_analysis_April2017_emxutil.obj rt_nonfinite.obj rtGetNaN.obj rtGetInf.obj

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = 

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = 

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_ = /openmp
CFLAGS_BASIC = $(DEFINES) 

CFLAGS = $(CFLAGS) $(CFLAGS_) $(CFLAGS_BASIC)

#-----------------
# C++ Compiler
#-----------------

CPPFLAGS_ = /openmp
CPPFLAGS_BASIC = $(DEFINES) 

CPPFLAGS = $(CPPFLAGS) $(CPPFLAGS_) $(CPPFLAGS_BASIC)

#---------------
# C++ Linker
#---------------

CPP_LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

CPP_LDFLAGS = $(CPP_LDFLAGS) $(CPP_LDFLAGS_)

#------------------------------
# C++ Shared Library Linker
#------------------------------

CPP_SHAREDLIB_LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS) $(CPP_SHAREDLIB_LDFLAGS_)

#-----------
# Linker
#-----------

LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

LDFLAGS = $(LDFLAGS) $(LDFLAGS_)

#--------------------------
# Shared Library Linker
#--------------------------

SHAREDLIB_LDFLAGS_ = /nodefaultlib:vcomp /LIBPATH:$(MATLAB_ROOT)\bin\$(ARCH) libiomp5md.lib

SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS) $(SHAREDLIB_LDFLAGS_)

###########################################################################
## INLINED COMMANDS
###########################################################################


!include $(MATLAB_ROOT)\rtw\c\tools\vcdefs.mak


###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute set_environment_variables


all : build
	@cmd /C "@echo ### Successfully generated all binary outputs."


build : set_environment_variables prebuild $(PRODUCT)


prebuild : 


download : build


execute : download


set_environment_variables : 
	@set INCLUDE=$(INCLUDES);$(INCLUDE)
	@set LIB=$(LIB)


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	$(PERL) $(GEN_LNK_SCRIPT) $(CMD_FILE) $(OBJS)
	@cmd /C "@echo ### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS) -out:$(PRODUCT) @$(CMD_FILE)
	@cmd /C "@echo ### Created: $(PRODUCT)"
	$(RM) $(CMD_FILE)


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(RELATIVE_PATH_TO_ANCHOR)}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(START_DIR)}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\IMT_analysis_April2017}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(START_DIR)\codegen\lib\IMT_analysis_April2017}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\rtw\c\src}.c.obj :
	$(CC) $(CFLAGS) -Fo"$@" "$<"


{$(MATLAB_ROOT)\rtw\c\src}.cpp.obj :
	$(CPP) $(CPPFLAGS) -Fo"$@" "$<"


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : $(MAKEFILE) rtw_proj.tmw


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@cmd /C "@echo ### PRODUCT = $(PRODUCT)"
	@cmd /C "@echo ### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@cmd /C "@echo ### BUILD_TYPE = $(BUILD_TYPE)"
	@cmd /C "@echo ### INCLUDES = $(INCLUDES)"
	@cmd /C "@echo ### DEFINES = $(DEFINES)"
	@cmd /C "@echo ### ALL_SRCS = $(ALL_SRCS)"
	@cmd /C "@echo ### ALL_OBJS = $(ALL_OBJS)"
	@cmd /C "@echo ### LIBS = $(LIBS)"
	@cmd /C "@echo ### MODELREF_LIBS = $(MODELREF_LIBS)"
	@cmd /C "@echo ### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@cmd /C "@echo ### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@cmd /C "@echo ### CFLAGS = $(CFLAGS)"
	@cmd /C "@echo ### LDFLAGS = $(LDFLAGS)"
	@cmd /C "@echo ### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@cmd /C "@echo ### CPPFLAGS = $(CPPFLAGS)"
	@cmd /C "@echo ### CPP_LDFLAGS = $(CPP_LDFLAGS)"
	@cmd /C "@echo ### CPP_SHAREDLIB_LDFLAGS = $(CPP_SHAREDLIB_LDFLAGS)"
	@cmd /C "@echo ### ARFLAGS = $(ARFLAGS)"
	@cmd /C "@echo ### MEX_CFLAGS = $(MEX_CFLAGS)"
	@cmd /C "@echo ### MEX_CPPFLAGS = $(MEX_CPPFLAGS)"
	@cmd /C "@echo ### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@cmd /C "@echo ### MEX_CPPLDFLAGS = $(MEX_CPPLDFLAGS)"
	@cmd /C "@echo ### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@cmd /C "@echo ### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@cmd /C "@echo ### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files..."
	@if exist $(PRODUCT) $(RM) $(PRODUCT)
	$(RM) $(ALL_OBJS)
	$(ECHO) "### Deleted all derived files."


