rm *.o
c++ -c -I $TECINCLUDE *cpp
c++ -o process *.o $TECLIB/libtecio.a $TECLIB/libadkutil.a -lstdc++ -pthread