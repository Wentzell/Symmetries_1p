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
#HEADERS := $(shell find $(HEADDIR) -type f -name '*.$(HEADEXT)') 
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Compiler Settings
CFLAGS := -std=c++11 #-fopenmp # General compiler flags
CFLAGS += -D NO_MOMENTA
DBFLAGS := -O0 -g # Compiler flags for debugging
PROFFLAGS := -O2 -g # Compiler flags for profiling
OPTIMIZEFLAGS := -flto -march=native -Ofast -fopenmp # Compiler flags for optimal speed
#OPTIMIZEFLAGS := -march=native -O3 # Compiler flags for optimal speed
LIB := -lgsl -lgslcblas -lm #-fopenmp
INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

debug: 	CFLAGS += $(DBFLAGS)
debug: 	$(TARGET)

prof: 	CFLAGS += $(PROFFLAGS)
prof: 	$(TARGET)

optimize: CFLAGS += $(OPTIMIZEFLAGS)
optimize: $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(INC) $(SOURCES) -o $(TARGET) $(LIB)

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean

