rm *.o process
c++ -std=c++11 -c -I $TECINCLUDE *cpp CAD2D/*.cpp -Og -g
c++ -std=c++11 -o process *.o libtecio.a libadkutil.a -lstdc++ -pthread -Og -g
