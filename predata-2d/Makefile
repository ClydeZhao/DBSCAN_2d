CCOMPILE=mpic++

PLATFORM=Linux-amd64-64

CPPFLAGS= -I$(HADOOP_HOME)/src/c++/libhdfs -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux -I /home/clyde/pregelplus 
LIB = -L$(HADOOP_HOME)/c++/$(PLATFORM)/lib
LDFLAGS = -lhdfs -Wno-deprecated -O2

all: run

run: run.cpp
	$(CCOMPILE) run.cpp $(CPPFLAGS) $(LIB) $(LDFLAGS)  -o run

clean:
	-rm run
