rm *.o process
if [[ $(uname -s) == Linux ]]
then
c++ -std=c++11 -c -I $TECINCLUDE *cpp CAD2D/*.cpp -Og -g
c++ -std=c++11 -o process *.o libtecio.a libadkutil.a -lstdc++ -pthread -Og -g
else
c++ -std=c++11 -c  *cpp CAD2D/*.cpp -Og -g
c++ -std=c++11 -o process *.o  -lstdc++ -pthread -Og -g
fi
