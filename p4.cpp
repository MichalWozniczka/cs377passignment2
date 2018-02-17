#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

using namespace std;

struct IntPoint {
	int r;
	vector<double>* p;
};

struct CSR {
	//vector<IntPoint> rp;
	//vector<vector<double>* > ci;

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

	double& weight(int x, int y) {
		for(int i = rp.at(x); i < rp.at(x+1); i++) {
			if(ci.at(i) == y) {
				return va.at(i);
			}
		}
		double d = -1;
		return d;
	}

	vector<int> neighbors(int n) {
		vector<int> neighbors;
		for(int i = rp.at(n); i < rp.at(n+1); i++) {
			neighbors.push_back(ci.at(i));
		}
		return neighbors;
	}

	void print() {
		int rowctr = 0;
		for(int i = 0; i < ci.size(); i++) {
			if(rowctr < rp.size() && rp.at(rowctr) < i) {
				cout << la.at(rowctr) << "\t" << rp.at(rowctr) << "\n";
				rowctr++;
			}
			if(rowctr < rp.size() && rp.at(rowctr) == i) {
				cout << la.at(rowctr) << "\t" << rp.at(rowctr) << "\t" << ci.at(i) << "\t" << va.at(i) << "\n";
				rowctr++;
			}
			else {
				cout << "\t\t" << ci.at(i) << "\t" << va.at(i) << "\n";
			}
		}
		for(int i = rowctr; i < rp.size(); i++) {
			cout << "\t" << rp.at(rowctr) << "\n";
		}

		cout << "\n";
	}
				
};

vector<CSR> dimacs_to_csr_and_transpose(ifstream& file) {
	string line;
	int edgeCount;
	int nodeCount;

	vector<vector<double>> matrix;

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
					matrix = vector<vector<double>>(nodeCount, vector<double>(nodeCount, -1));
					break;
				case 'a':
					int x = stoi(tokens.at(1));
					int y = stoi(tokens.at(2));
					int w = stoi(tokens.at(3));
					matrix.at(x-1).at(y-1) = w;
					break;
			}
		}
	}

	cout << edgeCount << " " << nodeCount << "\n";
	CSR csr(matrix);
	CSR tcsr(matrix, true);

	return {csr, tcsr};
}

CSR dimacs_to_csr(ifstream& file, bool transpose = true) {
	if(transpose) return dimacs_to_csr_and_transpose(file).at(1);
	return dimacs_to_csr_and_transpose(file).at(0);
}

void page_rank_init(CSR& graph) {
	for(int i = 0; i < graph.la.size(); i++) {
		graph.la.at(i) = 1.0 / graph.la.size();
	}
}

void page_rank_pull(CSR& graph, CSR& transpose) {
	double d = 0.85;
	page_rank_init(graph);

	for(int i = 0; i < 50; i++) {
		for(int j = 0; j < graph.rp.size()-1; j++) {
			vector<int> neighbors = graph.neighbors(j);
			graph.laNxt.at(j) = (1.0 - d) / graph.la.size();
			for(int k = 0; k < neighbors.size(); k++) {
				graph.laNxt.at(j) += d * (graph.la.at(neighbors.at(k)) / transpose.neighbors(neighbors.at(k)).size());
			}
		}
		graph.la = graph.laNxt;
	}
}

int main() {
	ifstream file;
	file.open("input.txt");
	vector<CSR> csr = dimacs_to_csr_and_transpose(file);
	csr.at(1).print();
	page_rank_pull(csr.at(1), csr.at(0));
	csr.at(1).print();

	double sum = 0;
	for(int i = 0; i < csr.at(1).la.size(); i++) {
		sum += csr.at(1).la.at(i);
	}
	for(int i = 0; i < csr.at(1).la.size(); i++) {
		csr.at(1).la.at(i) /= sum;
		cout << csr.at(1).la.at(i) << "\n";
	}
}
