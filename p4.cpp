#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

using namespace std;

struct IntPoint {
	int r;
	vector<int>* p;
};

struct CSR {
	vector<IntPoint> rp;
	vector<vector<int>* > ci;

	CSR() {
		rp;
		ci;
	}

	CSR(vector<vector<int>> matrix) {
		for(int i = 0; i < matrix.size(); i++) {
			bool containsEdge = false;
			IntPoint ip;
			ip.r = i;
			for(int j = 0; j < matrix.at(0).size(); j++) {
				if(matrix.at(i).at(j) != -1) {
					vector<int>* v = new vector<int>({j, matrix.at(i).at(j)});
					ci.push_back(v);
					
					if(containsEdge == false) {
						ip.p = v;
					}
					containsEdge = true;
				}
			}
			if(containsEdge) {
				rp.push_back(ip);
			}
		}
	}

	void print() {
		for(int i = 0; i < rp.size(); i++) {
			cout << rp.at(i).r << " " << rp.at(i).p << " " << rp.at(i).p->at(0) << "\n";
		}
		cout << "\n";
		for(int i = 0; i < ci.size(); i++) {
			cout << ci.at(i)->at(0) << " " << ci.at(i)->at(1) << "\n";
		}
	}
				
};

CSR dimacs_to_csr(ifstream& file) {
	string line;
	int edgeCount;
	int nodeCount;

	vector<vector<int>> matrix;

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
					matrix = vector<vector<int>>(nodeCount, vector<int>(nodeCount, -1));
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
	csr.print();

	return CSR();
}

int main() {
	ifstream file;
	file.open("input.txt");
	CSR csr = dimacs_to_csr(file);
}
