###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project 
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES         := no

CHEMICAL_REACTIONS  := no
CONTACT             := no
FLUID_PROPERTIES    := no
HEAT_CONDUCTION     := no
MISC                := no
NAVIER_STOKES       := no
PHASE_FIELD         := no
RDG                 := no
RICHARDS            := no
SOLID_MECHANICS     := no
STOCHASTIC_TOOLS    := no
TENSOR_MECHANICS    := no
WATER_STEAM_EOS     := no
XFEM                := no
POROUS_FLOW         := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := itk_testapp
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(FRAMEWORK_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
ADDITIONAL_LIBS := -L/opt/pozzis/itk/4.12/lib/ -lITKBiasCorrection-4.12.1 \
-lITKBioCell-4.12.1 \
-lITKCommon-4.12.1 \
-lITKDICOMParser-4.12.1 \
-litkdouble-conversion-4.12.1 \
-lITKEXPAT-4.12.1 \
-lITKFEM-4.12.1 \
-litkgdcmcharls-4.12.1 \
-litkgdcmCommon-4.12.1 \
-litkgdcmDICT-4.12.1 \
-litkgdcmDSED-4.12.1 \
-litkgdcmIOD-4.12.1 \
-litkgdcmjpeg12-4.12.1 \
-litkgdcmjpeg16-4.12.1 \
-litkgdcmjpeg8-4.12.1 \
-litkgdcmMEXD-4.12.1 \
-litkgdcmMSFF-4.12.1 \
-litkgdcmopenjpeg-4.12.1 \
-litkgdcmsocketxx-4.12.1 \
-litkgdcmuuid-4.12.1 \
-lITKgiftiio-4.12.1 \
-litkhdf5.1 \
-litkhdf5_cpp.1 \
-lITKIOBioRad-4.12.1 \
-lITKIOBMP-4.12.1 \
-lITKIOCSV-4.12.1 \
-lITKIOGDCM-4.12.1 \
-lITKIOGE-4.12.1 \
-lITKIOGIPL-4.12.1 \
-lITKIOHDF5-4.12.1 \
-lITKIOImageBase-4.12.1 \
-lITKIOIPL-4.12.1 \
-lITKIOJPEG-4.12.1 \
-lITKIOLSM-4.12.1 \
-lITKIOMesh-4.12.1 \
-lITKIOMeta-4.12.1 \
-lITKIOMRC-4.12.1 \
-lITKIONIFTI-4.12.1 \
-lITKIONRRD-4.12.1 \
-lITKIOPNG-4.12.1 \
-lITKIOSiemens-4.12.1 \
-lITKIOSpatialObjects-4.12.1 \
-lITKIOStimulate-4.12.1 \
-lITKIOTIFF-4.12.1 \
-lITKIOTransformBase-4.12.1 \
-lITKIOTransformHDF5-4.12.1 \
-lITKIOTransformInsightLegacy-4.12.1 \
-lITKIOTransformMatlab-4.12.1 \
-lITKIOVTK-4.12.1 \
-lITKIOXML-4.12.1 \
-litkjpeg-4.12.1 \
-lITKKLMRegionGrowing-4.12.1 \
-lITKLabelMap-4.12.1 \
-lITKMesh-4.12.1 \
-lITKMetaIO-4.12.1 \
-litknetlib-4.12.1 \
-litkNetlibSlatec-4.12.1 \
-lITKniftiio-4.12.1 \
-lITKNrrdIO-4.12.1 \
-lITKOptimizers-4.12.1 \
-lITKOptimizersv4-4.12.1 \
-lITKPath-4.12.1 \
-litkpng-4.12.1 \
-lITKPolynomials-4.12.1 \
-lITKQuadEdgeMesh-4.12.1 \
-lITKSpatialObjects-4.12.1 \
-lITKStatistics-4.12.1 \
-litksys-4.12.1 \
-litktestlib-4.12.1 \
-litktiff-4.12.1 \
-lITKTransform-4.12.1 \
-lITKTransformFactory-4.12.1 \
-litkv3p_netlib-4.12.1 \
-litkvcl-4.12.1 \
-lITKVideoCore-4.12.1 \
-lITKVideoIO-4.12.1 \
-litkvnl-4.12.1 \
-litkvnl_algo-4.12.1 \
-lITKVNLInstantiation-4.12.1 \
-lITKVTK-4.12.1 \
-lITKWatersheds-4.12.1 \
-litkzlib-4.12.1 \
-lITKznz-4.12.1 







app_INCLUDES += -I/opt/pozzis/itk/4.12/include/
