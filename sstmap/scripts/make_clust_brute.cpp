#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <math.h>

using namespace std;

/*
struct water {
    double oxygen[3];
    double hyd1[3];
    double hyd2[3];
    water() {
        for (int i = 0; i < 3; i++) {
            oxygen[i] = 0;
            hyd1[i] = 0;
            hyd2[i] = 0;
        }
    } //default constructor, will set all to 0
};
*/

int main (int argc, char** argv) {
    if (argc <= 1) {
        cerr << "\nUSAGE:\n\n"
        << "./bruteclust[-c clustercenterfile][-w within5Aofligand]\n"
        << "where:\n\n"
        << "clustercenterfile contains the coordinates of the chosen clusters\n"
        << "within5Aofligand contains all coordinates of waters within a certain distance of the ligand\n\n";
        exit(0);
    }

    int i = 0; string cfile; string wfile;
    while (i < argc) {
        if (!strcmp(argv[i], "-c")) {
            cfile = argv[++i];
        }
        if (!strcmp(argv[i], "-w")) {
            wfile = argv[++i];
        }
        i++;
    }

    if (cfile.empty() || wfile.empty()) {
        cerr << "Define the damn input files:\n"
        << "./bruteclust[-c clustercenterfile][-w within5Aofligand]\n"
        << "where:\n\n"
        << "clustercenterfile contains the coordinates of the chosen clusters\n"
        << "within5Aofligand contains all coordinates of waters within a certain distance of the ligand\n\n";
        exit(0);
    }


    vector<double > cens;
    vector<double > wats;


    string temp;
    ifstream cinput(cfile.c_str());
    getline(cinput, temp); //skip header
    while (!cinput.eof()) {
        getline(cinput, temp);
        if (!temp.empty()) {
            cens.push_back(atof(temp.substr(31, 7).c_str()));
            cens.push_back(atof(temp.substr(39, 7).c_str()));
            cens.push_back(atof(temp.substr(47, 7).c_str()));
        }
    }

    cinput.close();

    int numclust = cens.size()/3;

    ifstream winput(wfile.c_str());
    //getline(winput, temp); //skip header
    while (!winput.eof()) {
        getline(winput, temp);
        if (!temp.empty()) {
            wats.push_back(atof(temp.substr(31, 7).c_str()));
            wats.push_back(atof(temp.substr(39, 7).c_str()));
            wats.push_back(atof(temp.substr(47, 7).c_str()));
        }
    }

    winput.close();

    FILE* pFile;
    char fileName[80];
    int val;

    string name = "ATOM"; string atom = "H"; string resname = "T3P"; string chainid = "C"; int resseq = 1; double occupancy = 0.0; double T = 0.0;
    int pos = 0;
    //stringstream ss;
    //cout << "cencount" << cencount << endl;
    int j = 0;
    for (int i = 0; i < numclust; i++) {
        val = i+1;
        cout << val << endl;
        if (i < 9) {
            sprintf(fileName, "cluster.00000%i.pdb", val);
            //ss << "cluster.00000" << val << ".pdb";
            //ss >> fileName;
        }
        else if (i < 99){
            sprintf(fileName, "cluster.0000%i.pdb", val);
            //ss << "cluster.0000" << val << ".pdb";
            //ss >> fileName;
        }
	else {
	    sprintf(fileName, "cluster.000%i.pdb", val);
	}
        pFile = fopen(fileName, "w");
        pos = 0;
        /*for (int j = 0; j < watnum; j++) {
            if (wats[j].numclust != 0) {
                for (int k = 0; k < wats[j].numclust; k++) {
                    if (wats[j].clusters[k] == i) {
                        atom = "O";
                        fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, wats[j].oxygen[0], wats[j].oxygen[1], wats[j].oxygen[2], occupancy, T);
                        pos++;
                        atom = "H";
                        fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, wats[j].hyd1[0], wats[j].hyd1[1], wats[j].hyd1[2], occupancy, T);
                        pos++;
                        fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, wats[j].hyd2[0], wats[j].hyd2[1], wats[j].hyd2[2], occupancy, T);
                        pos++;
                    }
                }
            }
        }*/

        double dist = 10000;

        //for (int j = 0; j < cens.size(); j+=3) {
        j = i*3;
            for (int k = 0; k < wats.size(); k+=9) {
                dist = pow((cens[j] - wats[k]), 2) + pow((cens[j+1] - wats[k+1]), 2) + pow((cens[j+2] - wats[k+2]), 2);
                if (dist <= 4) {
                    atom = "O";
                    fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, wats[k], wats[k+1], wats[k+2], occupancy, T);
                    pos++;
                    atom = "H";
                    fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, wats[k+3], wats[k+4], wats[k+5], occupancy, T);
                    pos++;
                    fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, wats[k+6], wats[k+7], wats[k+8], occupancy, T);
                    pos++;
                }
            }
        //}

        fclose(pFile);
    }


}
