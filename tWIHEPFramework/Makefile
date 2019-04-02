####################################################################################
# Makefile for the IHEPAnalysisFrameWork package
# 
# to build a SingleTopRootAnalysis exectuable, simply type
# > gmake
#
# to clean up and delete all object and binary files and all temporary files, type
# > gmake clean
#
# Author: Huaqiao ZHANG
# date: 06/2015
#
####################################################################################

####################################################################################
# Options for Compilation
####################################################################################
# Compilers and linkers
CPP=g++
AR=ar ru
LD=g++
# Optimization or not?
#OPT= -g
OPT= -O3
startDir=$(shell pwd)
OSYS=$(shell uname -s)

# For MacOS (Darwin) include -bind_at_local
ifeq ($(OSYS),Darwin)
  SPECIALFLAGS= $(OPT) -bind_at_load #-fsanitize=address
else
  SPECIALFLAGS= $(OPT) #-fsanitize=address
endif

####################################################################################
# Flags that root demands:
####################################################################################

ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lTreePlayer -lMLP
#LIBS = externaltools/bunchcrossing/TrigBunchCrossingTool/standalone/libTrigBunchCrossingTool.so -lrt
# For gcc version 4 and greater use -iquote instead of -I.
CC-VERSION = $(shell gcc -v 2>/dev/stdout | grep "^gcc" | awk '{print $$3}' | sed -e "s/\.//g") 
VERSION-4 = $(shell if [ "$(call CC-VERSION)" -ge "400" ]; then echo "true"; else echo "false"; fi;)

####################################################################################
# Suffixes
####################################################################################
SUF_SRC= .cpp
SUF_INC= .hpp
SUF_OBJ= .o
SUF_DICT_INC= dict.hpp
SUF_DICT_OBJ= dict.o
# Shared library suffixes
ifeq ($(OSYS),Darwin)
  SUF_LIB= .dylib
else
  SUF_LIB= .so
endif
SUF_DEP= .d
SUF_DICT_DEP= dict.d
####################################################################################
# Directories
####################################################################################
DIR_SRC = src
DIR_INC = SingleTopRootAnalysis
DIR_OBJ = obj
DIR_LIB = lib
DIR_TMP = tmp
DIR_DEP= dep
####################################################################################
# Libraries
####################################################################################
LIB= SingleTopRootAnalysis

# Library flags are operating system dependent
ifeq ($(OSYS),Darwin)
  SO_OBJ_FLAG = -fno-common
  SO_LIB=  -dynamiclib -undefined dynamic_lookup -install_name lib$(LIB)$(SUF_LIB) -o lib/lib$(LIB)$(SUF_LIB)
else
  SO_OBJ_FLAG = -fPIC
  SO_LIB= -shared -Wl,-soname,$(LIB)$(SUF_LIB) -o lib/lib$(LIB)$(SUF_LIB)
endif

####################################################################################
# ldconfig 
####################################################################################
# Do not do ldconfig for MAC
ifneq ($(OSYS),Darwin)
  LDCONFIG= /sbin/ldconfig -n lib/
endif

####################################################################################
# Make compile flags from options
####################################################################################

ifeq ($(VERSION-4),true)
        CFLAGS = $(SPECIALFLAGS) -iquote- -iquote. -iquote${DIR_INC} -iquote${DIR_SRC}
else
        CFLAGS = $(SPECIALFLAGS) -I. -I${DIR_INC} -I${DIR_SRC}
endif

LFLAGS = $(SPECIALFLAGS)
RCXX=$(CFLAGS) $(ROOTCFLAGS) $(SO_OBJ_FLAG) $(CXXFLAGS) $(INCLUDES)
RLXX=$(LFLAGS) $(ROOTLIBS) $(LIBS) 


####################################################################################
# Source Directories
####################################################################################
DIR_BASE            = $(DIR_SRC)/Base/
DIR_BASE_DICTIONARY = $(DIR_SRC)/Base/Dictionary
DIR_PARTICLES       = $(DIR_SRC)/Particles/
DIR_TREES           = $(DIR_SRC)/Trees/
DIR_CUTS            = $(DIR_SRC)/Cuts/
DIR_VARS            = $(DIR_SRC)/Vars/
DIR_HIST            = $(DIR_SRC)/Histogramming/

####################################################################################
# Dictionary Files
####################################################################################
DICT_SRC =  $(subst $(DIR_BASE_DICTIONARY)   ,, $(shell find $(DIR_BASE_DICTIONARY) -name \*$(SUF_SRC) -type f ) ) \
            $(subst $(DIR_PARTICLES)         ,, $(shell find $(DIR_PARTICLES)       -name \*$(SUF_SRC) -type f ) ) \
            $(subst $(DIR_TREES)             ,, $(shell find $(DIR_TREES)           -name \*$(SUF_SRC) -type f ) ) 

DICT_OBJ_TMP =  $(subst $(SUF_SRC),$(SUF_DICT_OBJ),$(DICT_SRC))
DICT_INC_TMP =  $(subst $(SUF_SRC),$(SUF_DICT_INC),$(DICT_SRC))
DICT_DEP_TMP =  $(subst $(SUF_SRC),$(SUF_DICT_DEP),$(DICT_SRC))

DICT_OBJ =  $(subst $(DIR_SRC),$(DIR_OBJ),$(DICT_OBJ_TMP))
DICT_INC =  $(subst $(DIR_SRC),$(DIR_INC),$(DICT_INC_TMP))
DICT_DEP =  $(subst $(DIR_SRC),$(DIR_DEP),$(DICT_DEP_TMP))

####################################################################################
# Source Files
####################################################################################
SRC = $(subst $(DIR_BASE)      ,, $(shell find $(DIR_BASE)      -name \*$(SUF_SRC) -type f ) ) \
      $(subst $(DIR_PARTICLES) ,, $(shell find $(DIR_PARTICLES) -name \*$(SUF_SRC) -type f ) ) \
      $(subst $(DIR_TREES)     ,, $(shell find $(DIR_TREES)     -name \*$(SUF_SRC) -type f ) ) \
      $(subst $(DIR_CUTS)      ,, $(shell find $(DIR_CUTS)      -name \*$(SUF_SRC) -type f ) ) \
      $(subst $(DIR_VARS)      ,, $(shell find $(DIR_VARS)      -name \*$(SUF_SRC) -type f ) ) \
      $(subst $(DIR_HIST)      ,, $(shell find $(DIR_HIST)      -name \*$(SUF_SRC) -type f ) ) 

OBJ_TMP = $(subst $(SUF_SRC),$(SUF_OBJ),$(SRC))
INC_TMP = $(subst $(SUF_SRC),$(SUF_INC),$(SRC))
DEP_TMP = $(subst $(SUF_SRC),$(SUF_DEP),$(SRC))

OBJ = $(DICT_OBJ) $(subst $(DIR_SRC),$(DIR_OBJ),$(OBJ_TMP))
INC = $(subst $(DIR_SRC),$(DIR_INC),$(INC_TMP))
DEP = $(DICT_DEP) $(subst $(DIR_SRC),$(DIR_DEP),$(DEP_TMP))

####################################################################################
# Dependencies
####################################################################################
-include $(DEP)

####################################################################################
# Executables
####################################################################################
BINS_EventComparison= bin/EventComparison/EventComparisonMuon.x bin/EventComparison/EventComparisonElectron.x
BINS_Example= bin/example/test_analysis.x bin/example/test_analysis2.x
BINS_Matching= bin/Matching/MatchQuality.x
BINS_ParticleStudies= bin/ParticleStudies/electron_study.x bin/ParticleStudies/jet_study.x bin/ParticleStudies/MCParticle_study.x bin/ParticleStudies/muon_study.x
BINS_Wt= bin/Wt/Wt_generic.x
BINS_bstar = bin/bstar/bstar_generic.x
BINS_Vertex= bin/Wt/Wt_nVertOnly.x
BINS_Synch= bin/Wt/Wt_synch.x
#BINS_Wt= bin/Wt/Wt.x bin/Wt/Wt_mumu.x bin/Wt/Wteu.x bin/Wt/Wt_Tree.x bin/Wt/Wt_generic.x
BINS_Skimming= bin/Skimming/Skimming.x 
BINS_St= bin/St/St_subtractQCD.x
BINS_MCStudies= bin/MCstudies/MCttbarClassifier.x
BINS_BTag = bin/BTag/BTag_analysis.x bin/BTag/BTag_Electron_analysis.x bin/BTag/BTag_Muon_analysis.x

####################################################################################
# Target: Help
# Print message listing all options.  This is the default target.
####################################################################################
help:
	@echo
	@echo "------------------------------------"
	@echo "| Makefile for Single Top Analysis |"
	@echo "------------------------------------"
	@echo
	@echo "Compilation Options"
	@echo "-------------------"
	@echo "  make                  -->  Prints this message"
	@echo "  make  help            -->  Prints this message"
	@echo "  make  all             -->  Make all executables"
	@echo "  make  St              -->  St executables"
	@echo "  make  Wt              -->  Wt executables"
	@echo "  make  bstar              -->  bstar executables"
	@echo "  make  Skimming              -->  Skimming executables"
	@echo "  make  Example         -->  Example executable"
	@echo "  make  Matching        -->  Matching executables"
	@echo "  make  ParticleStudies -->  Particle study executables"
	@echo "  make  MCStudies       -->  MC-level studies executables"
	@echo "  make  EventComparison -->  Event Comparison executables"        
	@echo "  make  BTag            -->  Btagging study executables"       
	@echo
	@echo "Cleaning Options"
	@echo "----------------"
	@echo "  make clean                -->  Removes obj, lib, and executable files"
	@echo "  make cleanall             -->  Same as make clean"
	@echo "  make cleanExe             -->  Removes all executables"
	@echo "  make cleanSt              -->  Removes St executables"
	@echo "  make cleanWt              -->  Removes Wt executables"
	@echo "  make cleanbstar              -->  Removes bstar executables"
	@echo "  make cleanSkimming              -->  Removes Skimming executables"
	@echo "  make cleanExample         -->  Removes Example executables"
	@echo "  make cleanMatching        -->  Removes Matching executables"
	@echo "  make cleanParticleStudies -->  Removes Particle Studies executables"
	@echo "  make cleanMCStudies       -->  Removes MC Studies executables"
	@echo "  make cleanEventComparison -->  Removes Event Comparison executables"
	@echo "  make cleanBTag            -->  Removes BTag executables"	

####################################################################################
# Target: Debug
# Print debugging info
####################################################################################
dbg:
	@echo $(SRC)
	@echo ""
	@echo $(OBJ)
	@echo ""
	@echo $(INC) 
	@echo ""
	@echo $(DICT_INC) 
	@echo ""
	@echo $(DEP)

####################################################################################
# Target: Executable categories
# Make executables belonging to a certain category
####################################################################################
all: lib St Wt bstar Skimming FER2 Matching ParticleStudies MCStudies EventComparison BTag

St: lib $(BINS_St)

Wt: lib $(BINS_Wt)

bstar: lib $(BINS_bstar)

Synch: lib $(BINS_Synch)

nVert: lib $(BINS_Vertex)

Skimming: lib $(BINS_Skimming)

Example: lib $(BINS_Example)

Matching: lib $(BINS_Matching)

ParticleStudies: lib $(BINS_ParticleStudies)

MCStudies: lib $(BINS_MCStudies)

EventComparison: lib $(BINS_EventComparison)

BTag: lib $(BINS_BTag)


####################################################################################
# Target: Test
# For tests
####################################################################################
Test: $(OBJ)

####################################################################################
# Target: Dictionary 
# Make Dictionary files
####################################################################################
# Build rule for the dictionary files
$(DIR_OBJ)/%dict.o: $(DIR_SRC)/%$(SUF_SRC)
	@mkdir -p $(DIR_OBJ)/Base/CutFlow $(DIR_OBJ)/Base/Dictionary $(DIR_OBJ)/Base/Histograms 
	@mkdir -p $(DIR_OBJ)/Trees
	@mkdir -p $(DIR_OBJ)/Particles/Recon $(DIR_OBJ)/Particles/Truth $(DIR_OBJ)/Particles/TruthAll
	@mkdir -p $(DIR_OBJ)/Cuts/BoostedJet$(DIR_OBJ)/Cuts/Electron $(DIR_OBJ)/Cuts/IsolatedMuon $(DIR_OBJ)/Cuts/Jet $(DIR_OBJ)/Cuts/Lepton $(DIR_OBJ)/Cuts/Muon $(DIR_OBJ)/Cuts/Other $(DIR_OBJ)/Cuts/TaggedJet 
	@mkdir -p $(DIR_OBJ)/Histogramming/Matching  $(DIR_OBJ)/Histogramming/Recon  $(DIR_OBJ)/Histogramming/Topological  $(DIR_OBJ)/Histogramming/Truth  $(DIR_OBJ)/Histogramming/TruthAll $(DIR_OBJ)/Histogramming/Other
	@mkdir -p $(DIR_OBJ)/Vars

	@mkdir -p $(DIR_TMP)
	@mkdir -p $(DIR_TMP)/Base/CutFlow $(DIR_TMP)/Base/Dictionary $(DIR_TMP)/Base/Histograms
	@mkdir -p $(DIR_TMP)/Trees
	@mkdir -p $(DIR_TMP)/Particles/Recon $(DIR_TMP)/Particles/Truth $(DIR_TMP)/Particles/TruthAll
	@mkdir -p $(DIR_TMP)/Cuts/BoostedJet$(DIR_TMP)/Cuts/Electron $(DIR_TMP)/Cuts/IsolatedMuon $(DIR_TMP)/Cuts/Jet $(DIR_TMP)/Cuts/Lepton $(DIR_TMP)/Cuts/Muon $(DIR_TMP)/Cuts/Other $(DIR_TMP)/Cuts/TaggedJet $(DIR_TMP)/Cuts/Tau $(DIR_TMP)/Cuts/Weights
	@mkdir -p $(DIR_TMP)/Histogramming/Matching  $(DIR_TMP)/Histogramming/Recon  $(DIR_TMP)/Histogramming/Topological  $(DIR_TMP)/Histogramming/Truth  $(DIR_TMP)/Histogramming/TruthAll $(DIR_TMP)/Histogramming/Other
	@mkdir -p $(DIR_TMP)/Vars

	@mkdir -p $(DIR_DEP)
	@mkdir -p $(DIR_DEP)/Base/CutFlow $(DIR_DEP)/Base/Dictionary $(DIR_DEP)/Base/Histograms
	@mkdir -p $(DIR_DEP)/Trees
	@mkdir -p $(DIR_DEP)/Particles/Recon $(DIR_DEP)/Particles/Truth $(DIR_DEP)/Particles/TruthAll
	@mkdir -p $(DIR_DEP)/Cuts/BoostedJet$(DIR_DEP)/Cuts/Electron $(DIR_DEP)/Cuts/IsolatedMuon $(DIR_DEP)/Cuts/Jet $(DIR_DEP)/Cuts/Lepton $(DIR_DEP)/Cuts/Muon $(DIR_DEP)/Cuts/Other $(DIR_DEP)/Cuts/TaggedJet $(DIR_DEP)/Cuts/Tau $(DIR_DEP)/Cuts/Weights
	@mkdir -p $(DIR_DEP)/Histogramming/Matching  $(DIR_DEP)/Histogramming/Recon  $(DIR_DEP)/Histogramming/Topological  $(DIR_DEP)/Histogramming/Truth  $(DIR_DEP)/Histogramming/TruthAll $(DIR_DEP)/Histogramming/Other
	@mkdir -p $(DIR_DEP)/Vars

# Create dependency files                                                                                                                                                                       
	$(CPP) $(RCXX)  -MM  $(@:$(DIR_OBJ)/%dict$(SUF_OBJ)=$(DIR_SRC)/%$(SUF_SRC)) >  $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP))
# Move .d file to a temporary file                                                                                                                                                              
	@mv -f $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)) $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp
# Add path to .o file in target of dependency
	@sed -e 's|.*:|$@:|' < $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp > $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP))
# Treat missing .hpp files as changed so target is rebuilt                                                                                                                                      
	@sed -e 's/.*://' -e 's/\\$$//' < $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp  | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP))
# Remove temporary files                                                                                                                                                                        
	@rm -f $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp



	rootcint -f $(DIR_TMP)/$(@:$(DIR_OBJ)/%.o=%$(SUF_SRC)) -c -p -I${DIR_INC}  $(DIR_INC)/$(@:$(DIR_OBJ)/%dict$(SUF_OBJ)=%$(SUF_INC))
	$(CPP) $(RCXX) -I$(DIR_TMP)  -c $(DIR_TMP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_SRC)) -o $@

####################################################################################
# Target: Object
# Compile source code into object code
####################################################################################
# Build rule for object files
$(DIR_OBJ)/%.o: $(DIR_SRC)/%$(SUF_SRC)
	@mkdir -p $(DIR_OBJ)/Base/CutFlow $(DIR_OBJ)/Base/Dictionary $(DIR_OBJ)/Base/Histograms
	@mkdir -p $(DIR_OBJ)/Trees
	@mkdir -p $(DIR_OBJ)/Particles/Recon $(DIR_OBJ)/Particles/Truth $(DIR_OBJ)/Particles/TruthAll 
	@mkdir -p $(DIR_OBJ)/Cuts/BoostedJet$(DIR_OBJ)/Cuts/Electron $(DIR_OBJ)/Cuts/IsolatedMuon $(DIR_OBJ)/Cuts/Jet $(DIR_OBJ)/Cuts/Lepton $(DIR_OBJ)/Cuts/Muon $(DIR_OBJ)/Cuts/Other $(DIR_OBJ)/Cuts/TaggedJet $(DIR_OBJ)/Cuts/Tau $(DIR_OBJ)/Cuts/Weights
	@mkdir -p $(DIR_OBJ)/Histogramming/Matching  $(DIR_OBJ)/Histogramming/Recon  $(DIR_OBJ)/Histogramming/Topological  $(DIR_OBJ)/Histogramming/Truth  $(DIR_OBJ)/Histogramming/TruthAll $(DIR_OBJ)/Histogramming/Other

	@mkdir -p $(DIR_DEP)
	@mkdir -p $(DIR_DEP)/Base/CutFlow $(DIR_DEP)/Base/Dictionary $(DIR_DEP)/Base/Histograms
	@mkdir -p $(DIR_DEP)/Trees
	@mkdir -p $(DIR_DEP)/Particles/Recon $(DIR_DEP)/Particles/Truth $(DIR_DEP)/Particles/TruthAll
	@mkdir -p $(DIR_DEP)/Cuts/BoostedJet$(DIR_DEP)/Cuts/Electron $(DIR_DEP)/Cuts/IsolatedMuon $(DIR_DEP)/Cuts/Jet $(DIR_DEP)/Cuts/Lepton $(DIR_DEP)/Cuts/Muon $(DIR_DEP)/Cuts/Other $(DIR_DEP)/Cuts/TaggedJet $(DIR_DEP)/Cuts/Tau $(DIR_DEP)/Cuts/Weights
	@mkdir -p $(DIR_DEP)/Histogramming/Matching  $(DIR_DEP)/Histogramming/Recon  $(DIR_DEP)/Histogramming/Topological $(DIR_DEP)/Histogramming/Truth  $(DIR_DEP)/Histogramming/TruthAll $(DIR_DEP)/Histogramming/Other

        # Create dependency files
	$(CPP) $(RCXX) -MM  $(@:$(DIR_OBJ)/%$(SUF_OBJ)=$(DIR_SRC)/%$(SUF_SRC)) >  $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP))
        # Move .d file to a temporary file
	@mv -f $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)) $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp
        # Add path to .o file in target of dependency
	@sed -e 's|.*:|$@:|' < $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp > $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP))
        # Treat missing .hpp files as changed so target is rebuilt
	@sed -e 's/.*://' -e 's/\\$$//' < $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp  | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP))
        # Remove temporary files
	@rm -f $(DIR_DEP)/$(@:$(DIR_OBJ)/%$(SUF_OBJ)=%$(SUF_DEP)).tmp

        # Compile object files
	$(CPP) $(RCXX) -c $(@:$(DIR_OBJ)/%$(SUF_OBJ)=$(DIR_SRC)/%$(SUF_SRC)) -o $@

####################################################################################
# Target: Library
# Call target for making shared library
####################################################################################
lib: $(DIR_LIB)/lib$(LIB)$(SUF_LIB)
	@rm -f ${DIR_LIB}/*.pcm
#	@ln -s $(startDir)/${DIR_TMP}/*/*.pcm ${DIR_LIB}
#	@ln -s $(startDir)/${DIR_TMP}/*/*/*.pcm ${DIR_LIB}
####################################################################################
# Target: Library
# Make shared library
####################################################################################
lib/lib$(LIB)$(SUF_LIB): $(OBJ)
	@mkdir -p $(DIR_LIB)
	@mkdir -p $(DIR_OBJ)/Base/CutFlow $(DIR_OBJ)/Base/Dictionary $(DIR_OBJ)/Base/Histograms
	@mkdir -p $(DIR_OBJ)/Trees
	@mkdir -p $(DIR_OBJ)/Particles/Recon $(DIR_OBJ)/Particles/Truth $(DIR_OBJ)/Particles/TruthAll
	@mkdir -p $(DIR_OBJ)/Cuts/BoostedJet$(DIR_OBJ)/Cuts/Electron $(DIR_OBJ)/Cuts/IsolatedMuon $(DIR_OBJ)/Cuts/Jet $(DIR_OBJ)/Cuts/Lepton $(DIR_OBJ)/Cuts/Muon $(DIR_OBJ)/Cuts/Other $(DIR_OBJ)/Cuts/TaggedJet $(DIR_OBJ)/Cuts/Tau $(DIR_OBJ)/Cuts/Weights
	@mkdir -p $(DIR_OBJ)/Histogramming/Matching  $(DIR_OBJ)/Histogramming/Recon  $(DIR_OBJ)/Histogramming/Topological  $(DIR_OBJ)/Histogramming/Truth  $(DIR_OBJ)/Histogramming/TruthAll $(DIR_OBJ)/Histogramming/Other

#	$(AR) lib/lib$(LIB).a $(OBJ)
        # Use shared libraries          
	$(CPP) $(SO_LIB) $(OBJ) 

####################################################################################
# Target: Executable
# Link object files and libraries to form binary executable
####################################################################################
%.x: %$(SUF_SRC) lib/lib$(LIB)$(SUF_LIB) $(INC)
	$(CPP) $(RCXX) -c $<  -o $(DIR_OBJ)/bintmp$(SUF_OBJ) 
	$(LD) -o $@ $(RLXX) $(DIR_OBJ)/bintmp$(SUF_OBJ) -Llib -l$(LIB)  
	$(LDCONFIG)

####################################################################################
# Target: Cleaning
# Clean executables, object files, and libraries
####################################################################################
cleanEventComparison:    
	rm -f $(BINS_EventComparison)
cleanExample:
	rm -f $(BINS_Example)
cleanMatching:
	rm -f $(BINS_Matching)          
cleanParticleStudies:
	rm -f $(BINS_ParticleStudies)
cleanMCStudies:
	rm -f $(BINS_MCStudies)
cleanBTag:
	rm -f $(BINS_BTag)
cleanWt:
	rm -f $(BINS_Wt)
cleanbstar:
	rm -f $(BINS_bstar)
cleanSynch:
	rm -f $(BINS_Synch)
cleanSkimming:
	rm -f $(BINS_Skimming)
cleanSt:
	rm -f $(BINS_St)

cleanExe: cleanSt cleanWt cleanbstar cleanSynch cleanSkimming cleanEventComparison cleanExample cleanParticleStudies cleanMCStudies

cleanall: clean

clean:  cleanSt cleanWt cleanbstar cleanSynch cleanSkimming cleanEventComparison cleanExample cleanParticleStudies cleanMCStudies
	rm -rf ti_files
	rm -rf $(DIR_OBJ)/*
	rm -rf $(DIR_TMP)/*
	rm -f $(DIR_LIB)/*
	rm -rf $(DIR_DEP)/*
	rm -f core
	rm -f */*dict$(SUF_SRC)
	rm -f G__auto*inkDef.h
