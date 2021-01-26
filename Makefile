CC=c++
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CCFLAGS += -I $TECINCLUDE -Og -g
	LIB += libtecio.a libadkutil.a -lstdc++ -pthread
else
    CCFLAGS += -Og -g
	LIB += -lstdc++ -pthread
endif
CFLAGS=-I.
DEPS = %.h
OBJ = main.o airfoil.o Dataprocessing.o PlungingMotion.o Tecplotwraper.o Body.o IncFlow.o Util.o StructuredData.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

process: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIB)