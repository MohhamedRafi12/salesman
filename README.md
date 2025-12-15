Starter code and data for traveling salesman problem


Files in this directory:

* datareader.cpp : example code to read in the data files (use Makefile)
* datareader.py  : example code to read in the data files
* cities23.dat : list of coordinates for 23 cities in North America
* cities150.dat : 150 cities in North America
* cities1k.dat : 1207 cities in North America
* cities2k.dat : 2063 cities around the world
* routeplot.py : code to plot the globe and salesman's path<br>
usage:<br>
python routeplot.py cities.dat [cities2.dat] -r [="NA"],"World"'<br>
NA = North America, World = Mercator projection of the whole earth
* earth.C : (just for fun) plotting the globe in ROOT



TO BUILD:
~~~
make 
~~~ 
This generates the results for 150,1k,2k datafiles and produces out puts

**Include a short description of you methods used to select trial configurations.**
Given the initial route, I produced a NN baseline configuration and then did a simulated annealing approach with multiple runs. 

cities150.data
Initial file-order tour length: 317298.645 km
FINAL Comparison (lower is better):
  Nearest Neighbor: 52617.304 km
  Best SA (multi):  48369.045 km
  Improvement:      4248.258 km  (8.074%)
Total execution time: 1.153 s

cities1k.data
Initial file-order tour length: 732177.737 km
FINAL Comparison (lower is better):
  Nearest Neighbor: 112854.217 km
  Best SA (multi):  96685.249 km
  Improvement:      16168.967 km  (14.327%)
Total execution time: 174.042 s

cities2k.dat
Initial file-order tour length: 10187617.637 km
FINAL Comparison (lower is better):
  Nearest Neighbor: 328476.234 km
  Best SA (multi):  281186.031 km
  Improvement:      47290.204 km  (14.397%)
Total execution time: 280.559 s