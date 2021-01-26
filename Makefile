CC=c++
UNAME_S := $(shell uname -s)
CFLAGS=-std=c++11 -Og -g
ifeq ($(UNAME_S),Linux)
    CFLAGS += -I $(TECINCLUDE)
	LIBS += libtecio.a libadkutil.a -lstdc++ -pthread
else
	LIBS += -lstdc++ -pthread
endif
DEPS = %.h
OBJ =  CAD2D/airfoil.o Dataprocessing.o PlungingMotion.o Tecplotwraper.o Body.o IncFlow.o Util.o StructuredData.o main.o

main.o: main.cpp
	${CC} $(CFLAGS) -c $<

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $<  -o $@

process: $(OBJ)
	$(CC)  -o $@ $^ $(CFLAGS) $(LIBS)