###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Required Make variables
# APPLICATION_NAME  - the name of this application (all lower case)
# MOOSE_DIR	        - location of the MOOSE framework
#
# Optional Environment variables
# CURR_DIR	- current directory (DO NOT MODIFY THIS VARIABLE)
#
# Note: Make sure that there is no whitespace after the word 'yes' if enabling
# an application
###############################################################################
ROOT_DIR        ?= $(shell dirname `pwd`)

ifeq ($(MOOSE_DEV),true)
	MOOSE_DIR ?= $(ROOT_DIR)/devel/moose
else
	MOOSE_DIR ?= $(ROOT_DIR)/moose
endif

################################## ELK MODULES ################################
PHASE_FIELD       := yes
SOLID_MECHANICS   := yes
TENSOR_MECHANICS  := yes
HEAT_CONDUCTION   := yes
MISC              := yes
###############################################################################

# framework
include $(MOOSE_DIR)/build.mk
include $(MOOSE_DIR)/moose.mk

# modules
ELK_DIR ?= $(ROOT_DIR)/elk
include $(ELK_DIR)/elk.mk

# dep apps
APPLICATION_DIR    := $(ROOT_DIR)/hyrax
APPLICATION_NAME   := hyrax
BUILD_EXEC         := yes
DEP_APPS           := $(shell $(MOOSE_DIR)/scripts/find_dep_apps.py $(APPLICATION_NAME))
include            $(MOOSE_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
