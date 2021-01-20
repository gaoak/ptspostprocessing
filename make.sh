rm *.o process
c++ -Og -c -I $TECINCLUDE *cpp -g
c++ -Og -o process *.o libtecio.a libadkutil.a -lstdc++ -pthread -g
./process