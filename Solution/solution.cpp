#include <bits/stdc++.h>
#include <math.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

// change formula here
#define formulaString "(a * i + b * j)^2 + 1"
#define formula(a, b, k, l) (pow(a * k + l * b, 2) + 1)
using namespace std;

// undirected graph
class Graph {
   public:
    int V = 0;  // number of vertices
    int degeneracyValue = 0;
    int cycle = 0;
    int colourNumber = 0;
    list<int> *adj;
    int *matrix;
    int *path;
    int *colour;

    bool isComplete = true;
    bool isNull = true;

    // constructor and destructor
    Graph(int V) {
        this->V = V;
        adj = new list<int>[V];
        matrix = new int[V * V];
        path = new int[V];
        colour = new int[V];
    }

    ~Graph() {
        delete[] adj;
        delete[] matrix;
        delete[] path;
        delete[] colour;
    }

    void addEdge(int v, int w);  // function to add an edge to graph
    void graphColoring();        // prints coloring of the vertices
    void degeneracy();           // graph degeneracy
    void proofByExhaustion();    // proof by exhaustion method
    void printData();            // prints calculation data
};

// add edge to graph
void Graph::addEdge(int v, int w) {
    adj[v].push_back(w);
    adj[w].push_back(v);  // the graph is undirected
}

void Graph::printData() {
    // print adjacency matrix
    cout << "Adjacency matrix after all calculations: " << endl;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            cout << matrix[i * V + j] << " ";
        }

        cout << endl;
    }

    cout << endl;
    // graph is null
    if (isNull) {
        cout << "Graph is null-graph" << endl;  // graph is null-graph
        cout << endl;
        cout << "Chromatic number of given graph is 1" << endl;
        cout << endl;
        cout << "Degeneracy of null-graph is 0" << endl;
        return;  // no need for further calculations
    }

    cout << "Proof by exhaustion cycle length found: " << cycle << endl;
    cout << endl;
    if (cycle != 0) {
        cout << "First best path found: ";
        for (int i = 0; i < V; i++) {
            cout << path[i] << " -> ";
        }

        cout << path[0] << endl;
    } else {
        cout << "No path was found" << endl;
    }

    cout << endl;

    if (isComplete) {
        cout << "Graph is a complete graph" << endl;  // graph is complete
    } else {
        cout << "Graph is not a complete graph" << endl;  // graph is not complete
    }

    cout << endl;
    cout << "Chromatic number of given graph is: " << colourNumber << endl;
    for (int i = 0; i < V; i++) {
        cout << "Vertex " << i << " -> Colour " << colour[i] << endl;
    }
    cout << endl;
    cout << "Degeneracy of complete graphs is: " << degeneracyValue << endl;
}

// find degeneracy of a given graph
void Graph::degeneracy() {
    if (isComplete) {
        degeneracyValue = V;
    }

    int ans = INT_MAX;

    // for all vertices
    for (int i = 0; i < V; i++) {
        // set distance to max
        vector<int> dist(V, INT_MAX);
        vector<int> par(V, -1);

        // distance to source is 0
        dist[i] = 0;
        queue<int> q;

        // push in line
        q.push(i);

        // Nastavi sve dok queue ne bude prazan
        while (!q.empty()) {
            // take first in line
            int x = q.front();
            q.pop();

            // traverse in ascending order
            for (int child : adj[x]) {
                // if not visited
                if (dist[child] == INT_MAX) {
                    dist[child] = 1 + dist[x];
                    par[child] = x;
                    q.push(child);
                } else if (par[x] != child and par[child] != x) {  // is visited?
                    ans = min(ans, dist[x] + dist[child] + 1);
                }
            }
        }
    }

    // what degeneracy is found?
    if (ans == INT_MAX)
        degeneracyValue = 0;
    else
        degeneracyValue = ans;
}

// color the given graph
void Graph::graphColoring() {
    int result[V];
    int resultFinal[V];
    int num;
    int minColors = INT_MAX;

    // color from each vertice
    for (int vertices = 0; vertices < V; vertices++) {
        result[vertices] = 0;

        // remaining vertices are yet to be colored
        for (int u = 0; u < V; u++) {
            if (u != vertices)
                result[u] = -1;  // u has no color yet
        }

        // adjacent vertices available to be colored
        bool available[V];
        for (int cr = 0; cr < V; cr++)
            available[cr] = false;

        // assign colors to V-1 vertices
        for (int u = 0; u < V; u++) {
            if (u == vertices)
                continue;
            // go through all adjacent vertices and record their colors as unavailable
            list<int>::iterator it;
            for (it = adj[u].begin(); it != adj[u].end(); it++) {
                if (result[*it] != -1)
                    available[result[*it]] = true;
            }

            // find the first available color
            int cr;
            for (cr = 0; cr < V; cr++) {
                if (available[cr] == false)
                    break;
            }

            // no need to search next coloring
            if (cr > minColors)
                break;

            result[u] = cr;  // available color is found
            // reset the values back to false for the next iteration
            for (it = adj[u].begin(); it != adj[u].end(); it++) {
                if (result[*it] != -1)
                    available[result[*it]] = false;
            }
        }

        // store result
        num = result[0];
        for (int i = 1; i < V; i++) {
            if (result[i] > num)
                num = result[i];
        }

        // is the new minimum found?
        if (num < minColors) {
            memcpy(&resultFinal, &result, sizeof(result));
            minColors = num;
        }
    }

    // print final result
    colourNumber = minColors + 1;
    memcpy(colour, &resultFinal, sizeof(resultFinal));

    return;
}

// check if the edge should be added
bool isValid(int i, int j, vector<int> *values) {
    vector<int>::iterator it;

    // as soon it meets the requirement return true
    for (it = values->begin(); it < values->end(); it++) {
        if (abs(i - j) == *it) {  // requirement
            return true;
        }
    }

    return false;
}

// add edge to the given graph
void assembleGraph(vector<int> *values, int *n, Graph *g) {
    for (int i = 0; i < *n; i++)        // for each vertice i pair
        for (int j = 0; j < *n; j++)    // with vertice j
            if (isValid(i, j, values))  // check if e(ij) should be added
                g->addEdge(i, j);
}

// calculate edge weight values
void calcMatrix(Graph *g, int a, int b) {
    for (int k = 1; k < g->V + 1; k++)
        for (int l = 1; l < g->V + 1; l++)
            if (k == l)
                continue;
            else if (*(g->matrix + (k - 1) * g->V + l - 1) == 0)
                continue;
            else if (k < l)
                *(g->matrix + (k - 1) * g->V + l - 1) = formula(a, b, k, l);
            else
                *(g->matrix + (k - 1) * g->V + l - 1) = formula(b, a, k, l);

    return;
}

// load data from given file
Graph *readDataM(string filename) noexcept {
    ifstream stream;
    int n = 0;      // number of vertices
    int count = 0;  // data counter
    string str;

    stream.open(filename);  // open file stream

    if (!stream)  // stream == nullptr?
        return nullptr;

    // read number of vertices
    if (stream >> n) {
        if (n < 1)
            return nullptr;
    }

    Graph *g = new Graph(n);

    // load matrix and add edge to graph
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int num = -1;
            stream >> num;
            if (i == j)
                if (num != 0)
                    return nullptr;

            if (num < 0)
                return nullptr;

            g->matrix[i * n + j] = num;
            count++;
        }
    }

    stream.close();  // close file stream

    if (count != n * n)
        return nullptr;

    // construct graph, avoids double edges
    for (int i = 0, k = 0; i < n; i++, k = 0) {
        for (int j = i + 1; k < j; k++) {
            if (g->matrix[i * n + k] > 0) {
                g->addEdge(i, k);
            }
        }
    }

    // check if matrix is valid, if graph is null or complete
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                continue;
            }

            if (g->matrix[i * n + j] > 0) {  // graph is not null

                g->isNull = false;

                if (g->matrix[j * n + i] != g->matrix[i * n + j])
                    return nullptr;  // make sure M[i,j] == M [j, i]
            } else {
                if (g->matrix[j * n + i] != g->matrix[i * n + j])
                    return nullptr;

                g->isComplete = false;
            }
        }
    }

    return g;
}

// make matrix form elements in S
void generateMatrixFromS(Graph *g, vector<int> values) {
    int absVal;

    for (int i = 0; i < g->V; i++) {      // for i
        for (int j = 0; j < g->V; j++) {  // for j

            if (i == j) {
                *(g->matrix + i * g->V + j) = 0;  // set to 0
            } else {
                absVal = abs(i - j);  // subtract
                vector<int>::iterator it;

                for (it = values.begin(); it != values.end(); it++) {
                    if (absVal == *it) {                  // is it matching
                        *(g->matrix + i * g->V + j) = 1;  // set 1 if yes
                        break;                            // go to next
                    }

                    *(g->matrix + i * g->V + j) = 0;  // else set to 0
                }
            }
        }
    }

    return;
}

// load data from given file
Graph *readDataS(string filename) noexcept {
    ifstream stream;
    vector<int> values;  // elements in S
    vector<int>::iterator it;
    int num;
    int n;
    int count = 0;

    stream.open(filename);  // open file stream

    if (!stream)  // stream == nullptr?
        return nullptr;

    bool isDouble = false;   // avoid adding double edges
    while (stream >> num) {  // read numbers

        if (count == 0) {
            n = num;
            if (n < 1)
                return nullptr;

            count++;
        } else {
            isDouble = false;
            for (it = values.begin(); it != values.end(); it++) {  // for elements already in S
                if (num == *it) {                                  // find if S already contains num
                    isDouble = true;                               // if yes, break loop
                    break;
                }
            }

            if (!isDouble) {            // if num is not in S
                values.push_back(num);  // add it in vector
                count++;
            } else {
                continue;
            }
        }
    }

    Graph *g = new Graph(n);

    stream.close();  // close file stream
    if (count < 2)
        return nullptr;  // return false if number of data is < 3

    // check for null-graph and if numbers are eligible for complete graph
    num = 0;
    for (it = values.begin(); it < values.end(); it++) {  // in vector
        if (*it > 0 && *it < n) {                         // is number eligible for complete graph?
            num += *it;                                   // if yes add the number to sum
            g->isNull = false;                            // if the number is eligible for complete graph, graph is not null-graph
        }
    }

    // generate adjacency matrix from S
    generateMatrixFromS(g, values);

    // if not null-graph, check for complete graph
    // complete graph has at least n - 1 elements in S an if the sum of eligible numbers makes the sum of first n numbers
    if (values.size() < n - 1 && num != n * (n - 1) / 2) {
        g->isComplete = false;
    }

    // construct graph, avoids double edges
    for (int i = 0, k = 0; i < n; i++, k = 0) {
        for (int j = i + 1; k < j; k++) {
            if (g->matrix[i * n + k] > 0) {
                g->addEdge(i, k);
            }
        }
    }

    return g;
}

void Graph::proofByExhaustion() {
    int xyz[V];  // 1 - n
    int best[V];

    int minSum = INT_MAX;
    int sum = 0;

    for (int i = 0; i < V; i++) {
        xyz[i] = i;
    }

    // permutations
    do {
        list<int>::iterator it;
        bool cycleFound = false;
        bool lastAndFirst = false;
        for (int i = 1; i < V; i++) {                                        // check if vertice i
            cycleFound = false;                                              // assume cycle does not continue
            for (it = adj[xyz[i]].begin(); it != adj[xyz[i]].end(); it++) {  // is adjacent to vertice j
                if (i == V - 1) {                                            // for last vertice
                    if (*it == xyz[0]) {
                        lastAndFirst = true;  // last and first are adjacent
                        continue;             // continue checking if i and i - 1 are adjacent
                    } else if (*it == xyz[i - 1]) {
                        cycleFound = true;
                        if (lastAndFirst) {
                            break;
                        } else {
                            continue;
                        }
                    }
                } else if (*it == xyz[i - 1]) {
                    cycleFound = true;  // cycle does continue
                    break;
                }
            }

            if (!cycleFound)
                break;
        }

        if (!(cycleFound && lastAndFirst))
            continue;

        sum = 0;
        for (int i = 0; i < V - 1; i++) {
            sum += matrix[xyz[i] * V + xyz[i + 1]];  // sum permutations of vertices
        }

        sum += matrix[xyz[0] * V + xyz[V - 1]];  // sum last and first

        if (minSum > sum) {
            minSum = sum;  // if new sum is lower, set new min sum
            memcpy(&best, &xyz, sizeof(xyz));
        }
    } while (next_permutation(xyz + 1, xyz + V));  // take next permutation

    if (minSum != INT_MAX) {
        cycle = minSum;
        memcpy(path, &best, sizeof(best));
        return;
    }

    return;
}

// create graph function
Graph *createGraph(string filename, string method) {
    if (method == "matrix") {
        return readDataM(filename);
    } else if (method == "set") {
        return readDataS(filename);
    }

    return nullptr;
}

// test program
int main(int argc, char *argv[]) {
    int a = 0;
    int b = 0;
    string filename;
    string method;

    int i = 1;
    do {
        if (string(argv[i]) == "--f")
            filename = argv[++i];
        else if (string(argv[i]) == "--m")
            method = argv[++i];
        else if (string(argv[i]) == "--a")
            a = stoi(string(argv[++i]));
        else if (string(argv[i]) == "--b")
            b = stoi(string(argv[++i]));

        i++;
    } while (i < argc);

    cout << endl;
    cout << "Filename: " << filename << endl;
    cout << "Input method: " << method << endl;
    cout << "a = " << a << endl;
    cout << "b = " << b << endl;
    cout << endl;
    cout << "Formula: " << formulaString << endl;
    cout << endl;

    Graph *g = createGraph(filename, method);  // read from given file

    if (!g) {
        cout << "Error, check input arguments or input file!" << endl;
    } else {
        if ((a || b) && !g->isNull) {
            calcMatrix(g, a, b);
        }

        g->proofByExhaustion();
        g->graphColoring();
        g->degeneracy();
        g->printData();
    }

    return 0;
}