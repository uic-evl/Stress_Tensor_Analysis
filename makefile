TARGETS=Tensor

CXX=g++ -Wno-write-strings
CCFLAGS=-Wall -g -O3

# check the OS typing and set the flags depending
# on it
OSNAME := $(shell uname -s)
MACHNAME := $(shell uname -n | sed 's/[0-9].*//')
ifeq ($(OSNAME),Darwin)
	CCLIBS=-framework OpenGL -framework GLUT -lglew
else
	CCLIBS=-lGL -lGLU -lglut -lGLEW
endif

# get the source files 
SRCS = $(shell ls -R ./*cc ./lib/*/*cc)
# get all of the object names
OBJS = $(patsubst %.cc,%.o,$(patsubst %.cc,%.o,$(SRCS)))

# compile the assignment
Tensor: $(OBJS)
	$(CXX) $(CCFLAGS) -o $@ $^ $(CCLIBS)

# on make, compile Tensor
all: $(TARGETS)

clean:
	rm -f $(OBJS) $(TARGETS)
