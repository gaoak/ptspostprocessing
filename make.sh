rm *.o process
c++ -std=c++11 -c *cpp -Og -g
c++ -std=c++11 -o process *.o  -lstdc++ -Og -g