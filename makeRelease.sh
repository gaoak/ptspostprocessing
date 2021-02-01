rm *.o process
c++ -std=c++11 -DMAKEBASH -c *cpp -O3
c++ -std=c++11 -DMAKEBASH -o process *.o -lstdc++ -O3