CXX = g++
OBJS=NussRNAMain.o NussRNA.o
DEBUG=1
# -O0 -Wall
CFLAGS = -std=c++11
# LDFLAGS = 

NussRNA: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(CFLAGS) $(LDFLAGS)

%.o: %.cpp RNAFold.h
	$(CXX) -c $< $(CFLAGS) $(LDFLAGS)

clean:
	rm -f *.o