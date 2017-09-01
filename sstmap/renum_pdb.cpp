#include <fstream>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>

using namespace std;

int main (int argc, char** argv) {
	int i = 0; string infile;
	while (i < argc) {
		if (!strcmp(argv[i], "-i")) {
			infile = argv[++i];
		}
		i++;
	}

	int pos = 0; int watnum = 0;
	string temp;
	ifstream input(infile.c_str());
	ofstream output("fixedprob.pdb");
	while (!input.eof()) {
		getline(input, temp);
		if (!temp.empty()) {
			if (pos%3==0 && pos!=0) {
                                watnum++;
                        }
			if (pos < 10) {
				if (watnum < 10) {
					output << temp.substr(0,9) << " "  << pos << " " <<temp.substr(12,12) << " " << watnum << temp.substr(26,-1) << endl;
				}
				else {
					output << temp.substr(0,9) << " "  << pos << " " <<temp.substr(12,12) << watnum << temp.substr(26,-1) << endl;
				}
			}
			else if (pos < 100) {
				if (watnum <10) {
					output << temp.substr(0,9) << pos << " " << temp.substr(12,12) << " " <<watnum << temp.substr(26,-1) << endl;
				}
				else {
					output << temp.substr(0,9) << pos << " " << temp.substr(12,12) << watnum << temp.substr(26,-1) << endl;
				}
			}
			else {
				if (watnum <10) {
                                        output << temp.substr(0,9) << pos << " " << temp.substr(12,12) << " " <<watnum << temp.substr(26,-1) << endl;
                                }
                                else {
                                        output << temp.substr(0,9) << pos <<  temp.substr(12,12) << watnum << temp.substr(26,-1) << endl;
                                }
			}
			pos++;
		}
	}
	input.close();
	output.close();
}
