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

bool pairLessThan(Pair i, Pair j) {
	return i.x < j.x;
}

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

	CSR(vector<int> _rp, vector<int> _ci, vector<double> _va) {
		rp = _rp;
		ci = _ci;
		va = _va;
		la = vector<double>(rp.size());
		laNxt = la;
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

	CSR(vector<Pair> pairs, vector<bool> outgoing) {
		sort(pairs.begin(), pairs.end(), pairLessThan);

		for(int i = 0; i < pairs.size(); i++) {
			while(outgoing.at(rp.size()) == false) {
				rp.push_back(ci.size());
			}
			if(i == 0) {
				rp.push_back(ci.size());
			}
			else if(pairs.at(i).x != pairs.at(i-1).x) {
				rp.push_back(ci.size());
			}
			ci.push_back(pairs.at(i).y);
			va.push_back(pairs.at(i).w);
		}
		while(rp.size() < outgoing.size() && outgoing.at(rp.size()) == false) {
			rp.push_back(ci.size());
		}
		rp.push_back(ci.size());
		la = vector<double>(outgoing.size());
		laNxt = la;
	}

	CSR transpose() {
		vector<int> new_rp;
		vector<int> new_ci;
		vector<double> new_va;
		vector<Pair> pairs;
		vector<bool> outgoing(rp.size(), false);

		for(int i = 0; i < rp.size()-1; i++) {
			for(int j = rp.at(i); j < rp.at(i+1); j++) {
				pairs.push_back(Pair(ci.at(j), i, va.at(j)));
				outgoing.at(ci.at(j)) = true;
			}
		}
		
		sort(pairs.begin(), pairs.end(), pairLessThan);

		for(int i = 0; i < pairs.size(); i++) {
			while(outgoing.at(new_rp.size()) == false) {
				new_rp.push_back(new_ci.size());
			}
			if(i == 0) {
				new_rp.push_back(new_ci.size());
			}
			else if(pairs.at(i).x != pairs.at(i-1).x) {
				new_rp.push_back(new_ci.size());
			}
			new_ci.push_back(pairs.at(i).y);
			new_va.push_back(pairs.at(i).w);
		}
		while(new_rp.size() < rp.size() && outgoing.at(new_rp.size()) == false) {
			new_rp.push_back(new_ci.size());
		}

		return CSR(new_rp, new_ci, new_va);
	}


	int out_degree(int n, bool transpose=false) {
		if(!transpose) {
			return rp.at(n+1) - rp.at(n);
		}
		int out_count = 0;
		for(int i = 0; i < ci.size(); i++) {
			if(ci.at(i) == n) {
				out_count++;
			}
		}
		return out_count;
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
				cout << rp.at(rowctr) << "\t\t\t" << la.at(rowctr) << "\n";
				rowctr++;
			}
			if(rowctr < rp.size() && rp.at(rowctr) == i) {
				cout << rp.at(rowctr) << "\t" << ci.at(i) << "\t" << va.at(i) << "\t" << la.at(rowctr) << "\n";
				rowctr++;
			}
			else {
				cout << "\t" << ci.at(i) << "\t" << va.at(i) << "\n";
			}
		}
		for(int i = rowctr; i < rp.size(); i++) {
			cout << rp.at(rowctr) << "\n";
		}

		cout << "\n";
	}
				
};

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
					//matrix.at(x-1).at(y-1) = w;
					pairs.push_back(Pair(x-1, y-1, w));
					outgoing.at(x-1) = true;
					break;
			}
		}
	}

	cout << "File read. Edges: " << edgeCount << " " << " Vertices: " << nodeCount << "\n";
	cout << "Creating CSR model.\n";
	tbegin = clock();
	//CSR csr(matrix);
	CSR csr(pairs, outgoing);
	//csr.print();
	tend = clock();
	cout << "Created CSR model. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";
	cout << "Creating transposed CSR.\n";
	tbegin = clock();
	//CSR tcsr2(matrix, true);
	CSR tcsr = csr.transpose();
	//tcsr2.print();
	//tcsr.print();
	//tcsr.transpose().print();
	tend = clock();
	cout << "Created transposed CSR model. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";

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

void page_rank_pull(CSR& graph) {
	double d = 0.85;
	page_rank_init(graph);
	double max_diff = 1;

	CSR transpose = graph.transpose();
	//transpose.transpose();

	cout << "Beginning pull-style page rank.\n";
	tbegin = clock();

	while(max_diff > pow(10, -4)) {
		max_diff = 0;
		for(int j = 0; j < graph.rp.size()-1; j++) {
			vector<int> neighbors = graph.neighbors(j);
			graph.laNxt.at(j) = (1.0 - d) / graph.la.size();
			for(int k = 0; k < neighbors.size(); k++) {
				graph.laNxt.at(j) += d * (graph.la.at(neighbors.at(k)) / transpose.out_degree(neighbors.at(k)));
			}
			max_diff = max(max_diff, abs(graph.laNxt.at(j) - graph.la.at(j)));
		}
		graph.la = graph.laNxt;
	}

	double sum = 0;
	for(int i = 0; i < graph.la.size(); i++) {
		sum += graph.la.at(i);
	}
	for(int i = 0; i < graph.la.size(); i++) {
		graph.la.at(i) /= sum;
		//cout << graph.la.at(i) << "\n";
	}
	tend = clock();
	cout << "Pull-style page rank completed. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";
}

void page_rank_push(CSR& graph) {
	double d = 0.85;
	page_rank_init(graph);
	double max_diff = 1;

	cout << "Beginning push-style page rank.\n";
	tbegin = clock();

	while(max_diff > pow(10, -4)) {
		max_diff = 0;
		for(int j = 0; j < graph.laNxt.size(); j++) {
			graph.laNxt.at(j) = (1.0 - d) / graph.laNxt.size();
		}
		for(int j = 0; j < graph.rp.size()-1; j++) {
			vector<int> neighbors = graph.neighbors(j);
			for(int k = 0; k < neighbors.size(); k++) {
				graph.laNxt.at(neighbors.at(k)) += d * (graph.la.at(j) / graph.out_degree(j));
			}
		}
		for(int j = 0; j < graph.laNxt.size(); j++) {
			max_diff = max(max_diff, abs(graph.laNxt.at(j) - graph.la.at(j)));
		}
		graph.la = graph.laNxt;
	}

	double sum = 0;
	for(int i = 0; i < graph.la.size(); i++) {
		sum += graph.la.at(i);
	}
	for(int i = 0; i < graph.la.size(); i++) {
		graph.la.at(i) /= sum;
		//cout << graph.la.at(i) << "\n";
	}
	tend = clock();
	cout << "Push-style page rank completed. Time elapsed: " << double(tend - tbegin) / CLOCKS_PER_SEC << "\n";
}

int main() {
	ifstream file;
	file.open("input.txt");
	cout << "File opened.\n";
	vector<CSR> csr = dimacs_to_csr_and_transpose(file);
	page_rank_pull(csr.at(1));
	page_rank_push(csr.at(0));
}
