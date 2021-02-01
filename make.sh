rm *.o process
c++ -std=c++11 -DMAKEBASH -c *cpp -Og -g
c++ -std=c++11 -DMAKEBASH -o process *.o  -lstdc++ -Og -g