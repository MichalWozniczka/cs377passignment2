#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <cmath>

using namespace std;

clock_t tbegin;
clock_t tend;

//Stores the from- and to-node as well as the edge weight
struct Pair {
	int x;
	int y;
	double w;

	Pair() {
	}

	Pair(int _x, int _y) {
		x = _x;
		y = _y;
	}

	Pair(int _x, int _y, double _w) {
		x = _x;
		y = _y;
		w = _w;
	}
};

//Used in sorting Pair objects
bool pairLessThan(Pair i, Pair j) {
	return i.x < j.x;
}

//Memory representation of CSR matrix
struct CSR {

	vector<int> rp;
	vector<int> ci;
	vector<double> va;
	vector<double> la;
	vector<double> laNxt;

	CSR() {
		rp;
		ci;
		va;
	}

	//Creates a CSR matrix from an uncompressed matrix
	//O(n^2)
	//depricated
	CSR(vector<vector<double>> matrix) {
		for(int i = 0; i < matrix.size(); i++) {
			bool containsEdge = false;
			for(int j = 0; j < matrix.at(0).size(); j++) {
				if(matrix.at(i).at(j) != -1) {
					
					if(containsEdge == false) {
						rp.push_back(ci.size());
						containsEdge = true;
					}
					ci.push_back(j);
					va.push_back(matrix.at(i).at(j));
				}
			}
			if(!containsEdge) {
				rp.push_back(ci.size());
			}
		}
		rp.push_back(ci.size());
		la = vector<double>(matrix.size());
		laNxt = la;
	}

	//Creates a transposed CSR matrix from an uncompressed matrix
	//O(n^2)
	//depricated
	CSR(vector<vector<double>> matrix, bool transpose) {
		for(int i = 0; i < matrix.size(); i++) {
			bool containsEdge = false;
			for(int j = 0; j < matrix.at(0).size(); j++) {
				if(matrix.at(j).at(i) != -1) {
					
					if(containsEdge == false) {
						rp.push_back(ci.size());
						containsEdge = true;
					}
					ci.push_back(j);
					va.push_back(matrix.at(j).at(i));
				}
			}
			if(!containsEdge) {
				rp.push_back(ci.size());
			}
		}
		rp.push_back(ci.size());
		la = vector<double>(matrix.at(0).size());
		laNxt = la;
	}

	//Creates a CSR matrix from a list of (x, y) pairs where values in the matrix exist,
	//along with a list of booleans representing whether the node represented by the index
	//of the list has any outgoing edges.
	//O(n logn)
	CSR(vector<Pair> pairs, vector<bool> outgoing) {
		//Sorts pairs by x-value
		sort(pairs.begin(), pairs.end(), pairLessThan);

		for(int i = 0; i < pairs.size(); i++) {
			//Add a row pointer value for the very first pair and each pair where the
			//x-value is different from the last
			if(i == 0 || pairs.at(i).x != pairs.at(i-1).x) {
				while(rp.size() < outgoing.size() && outgoing.at(rp.size()) == false) {
					rp.push_back(ci.size());
				}
				rp.push_back(ci.size());
			}

			//Add a column index value and edge weight value no matter what
			ci.push_back(pairs.at(i).y);
			va.push_back(pairs.at(i).w);
		}

		//Tacks on extra row pointer values for nodes that didn't point to anything else
		while(rp.size() < outgoing.size() && outgoing.at(rp.size()) == false) {
			rp.push_back(ci.size());
		}
		rp.push_back(ci.size());

		//Create the label vector
		la = vector<double>(outgoing.size());
		laNxt = la;
	}

	//Returns the transpose of the current matrix in CSR form
	//O(n logn)
	CSR transpose() {
		vector<int> new_rp;
		vector<int> new_ci;
		vector<double> new_va;
		vector<Pair> pairs;
		vector<bool> outgoing(rp.size()-1, false);

		//Builds the list of bools representing whether the node represented by the index
		//of the list has any outgoing edges
		for(int i = 0; i < rp.size()-1; i++) {
			for(int j = rp.at(i); j < rp.at(i+1); j++) {
				pairs.push_back(Pair(ci.at(j), i, va.at(j)));
				outgoing.at(ci.at(j)) = true;
			}
		}
		
		return CSR(pairs, outgoing);
	}

	//Returns the out degree of the given node in the current matrix
	//O(1)
	int out_degree(int n) {
		return rp.at(n+1) - rp.at(n);
	}

	//Returns a list of the neighbors of the given node in the current matrix
	//O(n)
	vector<int> neighbors(int n) {
		vector<int> neighbors;
		for(int i = rp.at(n); i < rp.at(n+1); i++) {
			neighbors.push_back(ci.at(i));
		}
		return neighbors;
	}

	//Pretty-prints the matrix in CSR form
	void print() {
		int rowctr = 0;
		int colctr = 0;
		while(colctr < ci.size()) {
			while(rowctr < rp.size() && rp.at(rowctr) == rp.at(rowctr+1)) {
				cout << rp.at(rowctr) << "\t\t\t" << la.at(rowctr) << "\n";
				rowctr++;
			}
			if(rowctr < rp.size() && rp.at(rowctr) == colctr) {
				cout << rp.at(rowctr) << "\t" << ci.at(colctr) << "\t" << va.at(colctr) << "\t" << la.at(rowctr) << "\n";
				rowctr++;
				colctr++;
				while(colctr < rp.at(rowctr)) {
					cout << "\t" << ci.at(colctr) << "\t" << va.at(colctr) << "\n";
					colctr++;
				}
			}
		}
		for(int i = rowctr; i < rp.size(); i++) {
			cout << rp.at(rowctr) << "\n";
		}

		cout << "\n";
	}
				
};

//Converts a DIMACS file into a CSR matrix and its transpose and returns both in a vector
vector<CSR> dimacs_to_csr_and_transpose(ifstream& file) {
	string line;
	int edgeCount;
	int nodeCount;

	vector<vector<double>> matrix;
	vector<Pair> pairs;
	vector<bool> outgoing;

	while(getline(file, line)) {
		stringstream ss(line);
		string token;
		vector<string> tokens;
		while(getline(ss, token, ' ')) {
			tokens.push_back(token);
		}

		if(line.size() > 0) {
			switch(line.at(0)) {
				case 'c':
					break;
				case 'p':
					nodeCount = stoi(tokens.at(2));
					edgeCount = stoi(tokens.at(3));
					outgoing = vector<bool>(nodeCount, false);
					matrix = vector<vector<double>>(nodeCount, vector<double>(nodeCount, -1));
					break;
				case 'a':
					int x = stoi(tokens.at(1));
					int y = stoi(tokens.at(2));
					int w = stoi(tokens.at(3));
					pairs.push_back(Pair(x-1, y-1, w));
					outgoing.at(x-1) = true;
					break;
			}
		}
	}

	cout << "File read. Edges: " << edgeCount << " " << " Vertices: " << nodeCount << "\n";
	cout << "Creating CSR model.\n";
	tbegin = clock();
	CSR csr(pairs, outgoing);
	tend = clock();
	cout << "Created CSR model. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";
	cout << "Creating transposed CSR.\n";
	tbegin = clock();
	CSR tcsr = csr.transpose();
	tend = clock();
	cout << "Created transposed CSR model. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";

	return {csr, tcsr};
}

//Returns either a transposed CSR or a normal one
CSR dimacs_to_csr(ifstream& file, bool transpose = true) {
	if(transpose) return dimacs_to_csr_and_transpose(file).at(1);
	return dimacs_to_csr_and_transpose(file).at(0);
}

//Initializes all values of nodes to 1/N
void page_rank_init(CSR& graph) {
	for(int i = 0; i < graph.la.size(); i++) {
		graph.la.at(i) = 1.0 / graph.la.size();
	}
}

//Pull-style pagerank
void page_rank_pull(CSR& graph) {
	double d = 0.85;
	page_rank_init(graph);
	double max_diff = 1;

	CSR transpose = graph.transpose();

	cout << "Beginning pull-style page rank.\n";
	tbegin = clock();

	//Loop until no node changes by more than 10^-4 between iterations
	while(max_diff > pow(10, -4)) {
		max_diff = 0;
		for(int j = 0; j < graph.rp.size()-1; j++) {
			vector<int> neighbors = graph.neighbors(j);
			//PR_i+1(v) = (1 - d) / N
			graph.laNxt.at(j) = (1.0 - d) / graph.la.size();
			//For each neighbor u of v,
			for(int k = 0; k < neighbors.size(); k++) {
				//PR_i+1(v) += d * (PR_i(u) / out-degree(u))
				graph.laNxt.at(j) += d * (graph.la.at(neighbors.at(k)) / transpose.out_degree(neighbors.at(k)));
			}
			max_diff = max(max_diff, abs(graph.laNxt.at(j) - graph.la.at(j)));
		}
		//PR_i = PR_i+1
		graph.la = graph.laNxt;
	}

	//Scale every pagerank value so their sum is one
	double sum = 0;
	for(int i = 0; i < graph.la.size(); i++) {
		sum += graph.la.at(i);
	}
	for(int i = 0; i < graph.la.size(); i++) {
		graph.la.at(i) /= sum;
	}
	tend = clock();
	cout << "Pull-style page rank completed. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";
}

//Push-style pagerank
void page_rank_push(CSR& graph) {
	double d = 0.85;
	page_rank_init(graph);
	double max_diff = 1;

	cout << "Beginning push-style page rank.\n";
	tbegin = clock();

	//Loop until no node changes by more than 10^-4 between iterations
	while(max_diff > pow(10, -4)) {
		max_diff = 0;
		//PR_i+1(v) = (1 - d) / N
		for(int j = 0; j < graph.laNxt.size(); j++) {
			graph.laNxt.at(j) = (1.0 - d) / graph.laNxt.size();
		}
		for(int j = 0; j < graph.rp.size()-1; j++) {
			vector<int> neighbors = graph.neighbors(j);
			//For each neighbor u of v,
			for(int k = 0; k < neighbors.size(); k++) {
				//PR_i+1(v) += d * (PR_i(u) / out-degree(u))
				graph.laNxt.at(neighbors.at(k)) += d * (graph.la.at(j) / graph.out_degree(j));
			}
		}
		for(int j = 0; j < graph.laNxt.size(); j++) {
			max_diff = max(max_diff, abs(graph.laNxt.at(j) - graph.la.at(j)));
		}
		//PR_i = PR_i+1
		graph.la = graph.laNxt;
	}

	//Scale every pagerank value so their sum is one
	double sum = 0;
	for(int i = 0; i < graph.la.size(); i++) {
		sum += graph.la.at(i);
	}
	for(int i = 0; i < graph.la.size(); i++) {
		graph.la.at(i) /= sum;
	}
	tend = clock();
	cout << "Push-style page rank completed. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";
}

int main(int argc, char** argv) {
	ifstream file;
	file.open((argc ? argv[1] : "input.txt"));
	cout << "File opened.\n";
	vector<CSR> csr = dimacs_to_csr_and_transpose(file);
	file.close();
	page_rank_pull(csr.at(1));
	page_rank_push(csr.at(0));

	//Checks to make sure the push- and pull-style pagerank solutions are roughly similar
	for(int i = 0; i < csr.at(0).la.size(); i++) {
		if(abs(csr.at(0).la.at(i) - csr.at(1).la.at(i)) > pow(10, -4)) {
			cout << "Node " << i << " does not match: " << csr.at(0).la.at(i) << " " << csr.at(1).la.at(i) << "\n";
		}
	}

	//Writes solutions to files
	ofstream pull;
	pull.open("pull.txt");

	ofstream push;
	push.open("push.txt");

	for(int i = 0; i < csr.at(1).la.size(); i++) {
		pull << i << " " << csr.at(1).la.at(i) << "\n";
	}

	for(int i = 0; i < csr.at(0).la.size(); i++) {
		push << i << " " << csr.at(0).la.at(i) << "\n";
	}

	pull.close();
	push.close();
}
