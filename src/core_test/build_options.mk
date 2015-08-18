PWD=$(shell pwd)
EXE_NAME=test_model
NAMELIST_SUFFIX=test
override CPPFLAGS += -DCORE_TEST

ifneq ( , $(findstring filters, $(PRE_CORE_TARGETS) ) )
	override PRE_CORE_TARGETS += filters
	FCINCLUDES += -I$(PWD)/src/filters
	LIBS += -lfilters
endif

report_builds:
	@echo "CORE=test"
