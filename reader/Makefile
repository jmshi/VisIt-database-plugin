# General compiler specifications

CXX := icc
CPPFLAGS := 
CXXFLAGS := -O0 -g# -O3 
LDFLAGS  := 
LDLIBS   := 

# Preliminary definitions

# Preliminary definitions

EXE_DIR := 
EXECUTABLE := $(EXE_DIR)reader
SRC_FILES := globals.cpp \
	     io_wrapper.cpp \
	     parameter_input.cpp \
             meshblock_tree.cpp \
	     mesh.cpp \
	     main.cpp
	    
OBJ_FILES := globals.o \
	     io_wrapper.o \
	     parameter_input.o \
	     meshblock_tree.o \
	     mesh.o \
	     main.o


# Create objects from source files
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Link objects into executable
$(EXECUTABLE) : $(OBJ_FILES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $(OBJ_FILES) $(LDFLAGS) $(LDLIBS)


all : $(EXECUTABLE)

# Cleanup
clean :
	rm -rf *.o

#dependence
global.o:
io_wrapper.o:
parameter_input.o:
meshblock_tree.o:
mesh.o: global.o io_wrapper.o parameter_input.o meshblock_tree.o 
main.o: mesh.o parameter_input.o

