CC=c++
CFLAGS=-std=c++11 -Og -g
LIBS += -lstdc++
DEPS = %.h
OBJ =  airfoil.o Dataprocessing.o PlungingMotion.o FileIO.o Body.o IncFlow.o Util.o StructuredData.o main.o

main.o: main.cpp
	${CC} $(CFLAGS) -c $<

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $<  -o $@

process: $(OBJ)
	$(CC)  -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm *.o process