PWD=$(shell pwd)
EXE_NAME=test_model
NAMELIST_SUFFIX=test
USE_FILTERS=true
override CPPFLAGS += -DCORE_TEST
FCINCLUDES += -I$(PWD)/src/filters

report_builds:
	@echo "CORE=test"
