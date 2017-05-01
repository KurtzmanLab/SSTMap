#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include "6dimprobable.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc <= 1) {
        cerr << "\nUSAGE:\n\n"
        << "./6dimprob [-i inputfile]\n\n"
        << "where\n\n"
        << "inputfile is the file to read from (standard cluster)\n\n";
        //<< "x coordinate of center is from clustercenterfile\n"
        //<< "y coordinate of center is from clustercenterfile\n"
        //<< "z coordinate of center is from clustercenterfile\n\n";
        exit(0);
    }

    //double s = 0;

	ofstream probout("probcenters.pdb", ios::app);


    clock_t t;
    t = clock();
    int i = 0; string infile;
    while (i<argc) {
        if (!strcmp(argv[i], "-i")) {
            infile = argv[++i];
        }
        i++;
    }

    if (infile.empty()) {
        cerr << "infile needs to be defined\n"
        << "For full run instructions run executable with no arguments\n\n";
        exit(0);
    }

    double temp;
    string strtemp;
    
    vector<double > tmp4; //tmp4 contains all water atom coordinates
    ifstream stput(infile.c_str());
    getline(stput, strtemp); //skip header
    while (!stput.eof()) {
        getline(stput, strtemp);
        if (!strtemp.empty()) {
            temp = atof(strtemp.substr(31, 7).c_str());
            tmp4.push_back(temp);
            temp = atof(strtemp.substr(39, 7).c_str());
            tmp4.push_back(temp);
            temp = atof(strtemp.substr(47, 7).c_str());
            tmp4.push_back(temp);
        }
    }
	vector<double> tmp2; //storage for waters before adjusted for angles
    vector<double > tmp5; //tmp5 contains the oxygen x y z
    for (int i = 0; i < tmp4.size(); i++) {
        if (i%9 == 0 || i%9 == 1 || i%9 == 2) {
            tmp5.push_back(tmp4[i]);
        }
		tmp2.push_back(tmp4[i]);
    }

 
    /*
        Begin euler angle (quaternion) code
    */
    vector<double > tmp3; //will contain the quaternion description of the waters orientation
    double pi = 3.14159265359; double cenndist = 10000;
    int x_ref[3]; int y_ref[3]; int z_ref[3];
    x_ref[0] = 1; x_ref[1] = 0; x_ref[2] = 0;
    y_ref[0] = 0; y_ref[1] = 1; y_ref[2] = 0;
    z_ref[0] = 0; z_ref[1] = 0; z_ref[2] = 1;
    double ar[3]; double ar2; double h12; double h22; double h1length; double h2length; double arlength; double dotprohah1; double theta;
    double crossp_x_ref_h1[3]; double crossp_x_ref_h1_sign; double q[4]; double htemp[3]; double z_mol_vect[3];  double z_mol_vect2;
    double z_mol_vectlength; double dotproductz_mol_vectz_ref; double theta3p; double crossp_z_mol_vectz_ref[3]; double crossp_z_mol_vectz_ref_sign;
    double q2[4]; double e[4]; double singtest;
    for (int i = 0; i < tmp4.size(); i+=9) {
        tmp4[i+8] -= tmp4[i+2];
        tmp4[i+5] -= tmp4[i+2];
        tmp4[i+2] -= tmp4[i+2];
        tmp4[i+7] -= tmp4[i+1];
        tmp4[i+4] -= tmp4[i+1];
        tmp4[i+1] -= tmp4[i+1];
        tmp4[i+6] -= tmp4[i];
        tmp4[i+3] -= tmp4[i];
        tmp4[i] -= tmp4[i];
        h12 = pow(tmp4[i+3],2) + pow(tmp4[i+4],2) + pow(tmp4[i+5],2);
        h22 = pow(tmp4[i+6],2) + pow(tmp4[i+7],2) + pow(tmp4[i+8],2);
        h1length = pow(h12, 0.5);
        h2length = pow(h22, 0.5);
        if (tmp4[i+3] != 0) {
            tmp4[i+3] /= h1length;
        }
        if (tmp4[i+4] != 0) {
            tmp4[i+4] /= h1length;
        }
        if (tmp4[i+5] != 0) {
            tmp4[i+5] /= h1length;
        }
        if (tmp4[i+6] != 0) {
            tmp4[i+6] /= h1length;
        }
        if (tmp4[i+7] != 0) {
            tmp4[i+7] /= h1length;
        }
        if (tmp4[i+7] != 0) {
            tmp4[i+7] /= h1length;
        }
        ar[0]=tmp4[i+4]*x_ref[2] - tmp4[i+5]*x_ref[1];
        ar[1]=tmp4[i+5]*x_ref[0] - tmp4[i+3]*x_ref[2];
        ar[2]=tmp4[i+3]*x_ref[1] - tmp4[i+4]*x_ref[0];
        ar2 = pow(ar[0],2) + pow(ar[1],2) + pow(ar[2],2);
        arlength = pow(ar2, 0.5);
        if (ar[0] != 0) {
            ar[0] /= arlength;
        }
        if (ar[1] != 0) {
            ar[1] /= arlength;
        }
        if (ar[2] != 0) {
            ar[2] /= arlength;
        }
        dotprohah1 = 0;
        dotprohah1 += x_ref[0]*tmp4[i+3];
        dotprohah1 += x_ref[1]*tmp4[i+4];
        dotprohah1 += x_ref[2]*tmp4[i+5];
        theta = acos(dotprohah1);
        crossp_x_ref_h1[0]=tmp4[i+4]*x_ref[2] - tmp4[i+5]*x_ref[1];
        crossp_x_ref_h1[1]=tmp4[i+5]*x_ref[0] - tmp4[i+3]*x_ref[2];
        crossp_x_ref_h1[2]=tmp4[i+3]*x_ref[1] - tmp4[i+4]*x_ref[0];
        crossp_x_ref_h1_sign=crossp_x_ref_h1[0]*tmp4[i+3]+crossp_x_ref_h1[1]*tmp4[i+4]+crossp_x_ref_h1[2]*tmp4[i+5];
        if (crossp_x_ref_h1_sign > 0) {
            theta /=2;
        }
        else {
            theta /=-2;
        }
        q[0]=cos(theta);
        q[1]=ar[0]*sin(theta);
        q[2]=ar[1]*sin(theta);
        q[3]=ar[2]*sin(theta);

        htemp[0]= ((pow(q[0],2)+pow(q[1],2))-(pow(q[2],2)+pow(q[3],2)))* tmp4[i+3];
        htemp[0]= (2*(q[1]*q[2] + q[0]*q[3]) * tmp4[i+4] ) + htemp[0];
        htemp[0]= (2*(q[1]*q[3]-q[0]*q[2])*tmp4[i+5]) + htemp[0];
        htemp[1]= 2*( q[1]*q[2] - q[0]*q[3] ) * tmp4[i+3];
        htemp[1]= ( ( q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3] ) * tmp4[i+4] ) + htemp[1];
        htemp[1]= ( 2*( q[2]*q[3] + q[0]*q[1] ) * tmp4[i+5] ) + htemp[1];
        htemp[2]= 2*( q[1]*q[3] + q[0]*q[2]) * tmp4[i+3];
        htemp[2]= ( 2*( q[2]*q[3]-q[0]*q[1] ) * tmp4[i+4] ) + htemp[2];
        htemp[2]=  ( ( q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3] ) * tmp4[i+5] ) + htemp[2];
        tmp4[i+3]=htemp[0];
        tmp4[i+4]=htemp[1];
        tmp4[i+5]=htemp[2];
        htemp[0]= ((pow(q[0],2)+pow(q[1],2))-(pow(q[2],2)+pow(q[3],2)))* tmp4[i+6];
        htemp[0]= (2*(q[1]*q[2] + q[0]*q[3]) * tmp4[i+7] ) + htemp[0];
        htemp[0]= (2*(q[1]*q[3]-q[0]*q[2])*tmp4[i+8]) + htemp[0];
        htemp[1]= 2*( q[1]*q[2] - q[0]*q[3] ) * tmp4[i+6];
        htemp[1]= ( ( q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3] ) * tmp4[i+7] ) + htemp[1];
        htemp[1]= ( 2*( q[2]*q[3] + q[0]*q[1] ) * tmp4[i+8] ) + htemp[1];
        htemp[2]= 2*( q[1]*q[3] + q[0]*q[2]) * tmp4[i+6];
        htemp[2]= ( 2*( q[2]*q[3]-q[0]*q[1] ) * tmp4[i+7] ) + htemp[2];
        htemp[2]=  ( ( q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3] ) * tmp4[i+8] ) + htemp[2];
        tmp4[i+6]=htemp[0];
        tmp4[i+7]=htemp[1];
        tmp4[i+8]=htemp[2];
        z_mol_vect[0]=tmp4[i+4]*tmp4[i+8] - tmp4[i+5]*tmp4[i+7];
        z_mol_vect[1]=tmp4[i+5]*tmp4[i+6] - tmp4[i+3]*tmp4[i+8];
        z_mol_vect[2]=tmp4[i+3]*tmp4[i+7] - tmp4[i+4]*tmp4[i+6];
        z_mol_vect2= pow(z_mol_vect[0],2) + pow(z_mol_vect[1],2) + pow(z_mol_vect[2],2);
        z_mol_vectlength=pow(z_mol_vect2,0.5);
        if (z_mol_vect[0] !=0) {
            z_mol_vect[0] /= z_mol_vectlength;
        }
        if (z_mol_vect[1] !=0) {
            z_mol_vect[1] /= z_mol_vectlength;
        }
        if (z_mol_vect[2] !=0) {
            z_mol_vect[2] /= z_mol_vectlength;
        }
        dotproductz_mol_vectz_ref=0;
        for(int j=0;j<3;j++) {
            dotproductz_mol_vectz_ref+=z_mol_vect[j]*z_ref[j];
        }
        theta3p= acos(dotproductz_mol_vectz_ref);

        crossp_z_mol_vectz_ref[0]=z_mol_vect[1]*z_ref[2] - z_mol_vect[2]*z_ref[1];
        crossp_z_mol_vectz_ref[1]=z_mol_vect[2]*z_ref[0] - z_mol_vect[0]*z_ref[2];
        crossp_z_mol_vectz_ref[2]=z_mol_vect[0]*z_ref[1] - z_mol_vect[1]*z_ref[0];

        crossp_z_mol_vectz_ref_sign=crossp_z_mol_vectz_ref[0]*tmp4[i+3]+crossp_z_mol_vectz_ref[1]*tmp4[i+4]+crossp_z_mol_vectz_ref[2]*tmp4[i+5];

        if (crossp_z_mol_vectz_ref_sign < 0) {
            theta3p /=2;
        }
        else {
            theta3p /=-2;
        }

        q2[0]=cos(theta3p);
        q2[1]=x_ref[0]*sin(theta3p);
        q2[2]=x_ref[1]*sin(theta3p);
        q2[3]=x_ref[2]*sin(theta3p);

        e[0]= q[0]*q2[0]  -  q[1]*q2[1]  -  q[2]*q2[2]  -  q[3]*q2[3];
        e[1]= q[0]*q2[1]  +  q[1]*q2[0]  +  q[2]*q2[3]  -  q[3]*q2[2];
        e[2]= q[0]*q2[2]  -  q[1]*q2[3]  +  q[2]*q2[0]  +  q[3]*q2[1];
        e[3]= q[0]*q2[3]  +  q[1]*q2[2]  -  q[2]*q2[1]  +  q[3]*q2[0];

        tmp3.push_back(e[0]); tmp3.push_back(e[1]); tmp3.push_back(e[2]); tmp3.push_back(e[3]);

		
    }
	
	/*
		At this point we have tmp3 with quaternion values and tmp5 with oxygen coordinates
	*/ 

	vector<double> tmp;

	if (tmp5.size()/3 != tmp3.size()/4) {
		cout << "Somehow we have a different number of orientations than positions, back to the drawing board boys!\n\n";
	}
	
	for (int i = 0; i < tmp5.size(); i+=3) {
		int watpos = 0;
		tmp.push_back(tmp5[i]);
		tmp.push_back(tmp5[i+1]);
		tmp.push_back(tmp5[i+2]);
		watpos = i/3;
		watpos = watpos*4;
		tmp.push_back(tmp3[watpos]);
		tmp.push_back(tmp3[watpos+1]);
		tmp.push_back(tmp3[watpos+2]);
		tmp.push_back(tmp3[watpos+3]);

	}

	kdtree friendlytree(tmp);

	point pt; double* dists; int* winners;

	dists = new double[3];
	winners = new int[3];
	
	double dtemp;
	double windist = 10000;
	double winner[9];

	for (int i = 0; i < tmp.size(); i+=7) {
		pt.x[0] =  tmp[i];
		pt.x[1] =  tmp[i+1];
		pt.x[2] =  tmp[i+2];
		pt.x[3] =  tmp[i+3];
		pt.x[4] =  tmp[i+4];
		pt.x[5] =  tmp[i+5];
		pt.x[6] =  tmp[i+6];
		int watpos = i/7;
		friendlytree.nnearest(watpos, winners, dists, 3);
		watpos = watpos*9;
		for (int j = 0; j < 3; j++) {
			dtemp += dists[j];
		}
		dtemp = dtemp/3;
		if (dtemp < windist) {
			windist = dtemp;
			winner[0] = tmp2[watpos];
			winner[1] = tmp2[watpos+1];
			winner[2] = tmp2[watpos+2];
			winner[3] = tmp2[watpos+3];
			winner[4] = tmp2[watpos+4];
			winner[5] = tmp2[watpos+5];
			winner[6] = tmp2[watpos+6];
			winner[7] = tmp2[watpos+7];
			winner[8] = tmp2[watpos+8];

		}
		

	}

	char fileName[80];
	sprintf(fileName, "temp.dat");

	int pos = 0;
    FILE * pFile;
    pFile = fopen(fileName, "w");
    //Define the pdb file stuff
    string name = "ATOM"; string atom = "H"; string resname = "T3P"; string chainid = "C"; int resseq = 1; double occupancy = 0.0; double T = 0.0;
    fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos/1000, "O", resname.c_str(), chainid.c_str(), 0, winner[0], winner[1], winner[2], occupancy, T);
	fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos/1000, atom.c_str(), resname.c_str(), chainid.c_str(), 0, winner[3], winner[4], winner[5], occupancy, T);
	fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos/1000, atom.c_str(), resname.c_str(), chainid.c_str(), 0, winner[6], winner[7], winner[8], occupancy, T);

	fclose(pFile);

	string stemp;    
	ifstream input("temp.dat");
	getline(input, stemp); 
	probout << stemp << endl;
	getline(input, stemp); 
	probout << stemp << endl;
	getline(input, stemp); 
	probout << stemp << endl;

	delete dists;
	delete winners;
	input.close();
	probout.close();
}
    
