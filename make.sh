rm *.o process
c++ -Og -c -I $TECINCLUDE *cpp -g
c++ -Og -o process *.o $TECLIB/libtecio.a $TECLIB/libadkutil.a -lstdc++ -pthread -g
./process