# dismat-graph-colouring-degeneracy-pbe
C++ implementation of graph colouring and degeneracy algorithms with implemented proof by exhaustion method as part of Discrete Mathematics course at Faculty of Electrical Engineering and Computing, University of Zagreb, academic year 2021./2022.

<br />

**About:**

- find shortest cycle length of given graph using proof by exhaustion

- find graph degeneracy of given graph

- find chromatic number of given graph 

<br />

**Input arguments:**

- **'--f'** - relative path to data file

- **'--m'** - data file format method (**'matrix'** or **'set'**)

- **'--a'** - optional variable, default 0

- **'--b'** - variable variable, default 0

<br />

**Optional variables:**

- in some cases, optional variables **a** and **b** can be set

- if set (**a** or **b** are **>0**), formula will be applied to each weight 

<br />

**Matrix format:**
- graphs are defined with matrix using weights between vertices
- first row defines number of vertices
- custom weights between vertices can be set
- **'matrix_format.txt'** - matrix format examples

<br />

**Set format:**
- let there be a set of numbers **S**, two vertices **i** and **j** are adjecent only if **|i - j|** is element of **S**
- 1st row defines number of vertices
- 2nd row - set elements (eg. **'1 3 2 5 7'**)
- **'set_format.txt'** - set format examples
