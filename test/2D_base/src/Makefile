#shell definition
SHELL = /bin/bash

OPTIMIZATION_FLAG = -O

UFL_SOURCES = $(shell find . -name "*.ufl")
HEADER_FILES = $(UFL_SOURCES:.ufl=.h)
DEBUG_FILES = $(UFL_SOURCES:.ufl=_debug.py) 
DEBUG_FILES += $(UFL_SOURCES:.ufl=_debug.pyc) 

.DEFAULT_GOAL: all

.PHONY: all

all: $(HEADER_FILES)

%.h: %.ufl
	ffc -l dolfin $(OPTIMIZATION_FLAG) $<
	
cleanall: 
	$(RM) $(DEBUG_FILES) $(HEADER_FILES)
	
clean:
	$(RM) $(DEBUG_FILES)
