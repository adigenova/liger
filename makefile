ROOT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

#VPATH= ${ROOT_DIR}/obj 
CXX ?= g++
EXEC= liger
SRC_DIR := src
OBJ_DIR := obj

CXX ?= g++
EXEC= liger
LEMON= ${ROOT_DIR}/libs/lemon
EDLIB= ${ROOT_DIR}/libs/edlib/edlib/include
ISPOA= ${ROOT_DIR}/libs/spoa/include
#is necesary to install spoa lib
LSPOA= ${ROOT_DIR}/libs/spoa/build/lib

CFLAGS = -O3 -std=c++11 -lpthread -I ${LEMON}  -I ${ISPOA} -L ${LSPOA} -lspoa -I ${EDLIB}

SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))


ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

.PHONY: test clean deb


all: $(EXEC) 

clean:
	-rm -f  $(EXEC) $(OBJ_FILES)

$(EXEC): $(OBJ_FILES)
	$(CXX)  -o $@ $^ $(CFLAGS)

$(OBJ_DIR): 
	mkdir $(OBJ_DIR)
#creates on the flye the obj_dir if do not exist
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) -c -o $@ $< $(CFLAGS)


#init:
#	cd ${ROOT_DIR}/libs/spoa/
#	mkdir build
#	cmake ..
#	make
#	cd ${ROOT_DIR}
#	touch $@	


#spoalib:
#ifneq ($(wildcard ${ROOT_DIR}/libs/spoa/build/.),)
#	@echo "Found spoalib in current directory"
#else
#	@echo "Did not find spoalib. compiling"
#endif

