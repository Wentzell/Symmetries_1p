# MAKEFILE for fRG_O2 Project

# General
CC := g++ # This is the main compiler
SRCDIR := src
HEADDIR := include
BUILDDIR := build
TARGET := bin/run

# Sources and object files
SRCEXT := cpp
HEADEXT := h
SOURCES := $(shell find $(SRCDIR) -type f -name '*.$(SRCEXT)') 
HEADERS := $(shell find $(HEADDIR) -type f -name '*.$(HEADEXT)') 
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Compiler Settings
CFLAGS := -std=c++11 #-fopenmp # General compiler flags
DBFLAGS := -O0 -g # Compiler flags for debugging
PROFFLAGS := -O3 -g # Compiler flags for profiling
LIB := -lgsl -lgslcblas -lm #-fopenmp
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) -O3 $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) -O3 $(INC) -c -o $@ $<

debug: 	CFLAGS += $(DBFLAGS)
debug: 	$(TARGET)

prof: 	CFLAGS += $(PROFFLAGS)
prof: 	$(TARGET)

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean

