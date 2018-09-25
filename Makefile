TARGET = bcctools
BUILD_DIR = ./build
SRC_DIR = ./src

SRCS := $(shell find $(SRC_DIR) -type f -name *.cpp)
OBJS := $(patsubst $(SRC_DIR)/%,$(BUILD_DIR)/%,$(SRCS:.cpp=.o))

# Compiler
CXX = g++
CC = $(CXX)

# Set this to include SeqAn, SDSL and HTSlib, either system wide or download into
# current folder and set to .
SEQAN_INC = $(HOME)/local/include
SDSL_INC = $(HOME)/local/include
SDSL_LIB = $(HOME)/local/lib
HTSLIB_INC = $(HOME)/local/include

CXXFLAGS += -I$(SEQAN_INC) -DSEQAN_HAS_ZLIB=1 -DSEQAN_DISABLE_VERSION_CHECK -std=c++14
CXXFLAGS += -I$(SDSL_INC) -I$(HTSLIB_INC)

LDLIBS = -L$(SDSL_LIB) -lz -lpthread -lrt -lsdsl

# Date and version number from git
DATE := on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
VERSION := 0.0.1-$(shell git log --pretty=format:"%h" --date=iso | head -n 1)
CXXFLAGS += -DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"

# Enable warnings
CXXFLAGS += -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result
CXXFLAGS += -fno-stack-protector

# DEBUG build
#CXXFLAGS += -g -O0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=1

# RELEASE build
CXXFLAGS += -O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/%.h $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -f $(OBJS) $(TARGET)
