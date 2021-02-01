rm *.o process
c++ -std=c++11 -c *cpp -O3
c++ -std=c++11 -o process *.o -lstdc++ -O3