/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   main.cpp
 * Author: stevenramsey
 *
 * Created on March 16, 2016, 9:06 AM
 */


#include "probable.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc <= 1) {
        cerr << "\nUSAGE:\n\n"
        << "./probable [-i inputfile][-o outfile]\n\n"
        << "where\n\n"
        << "inputfile is the file to read from (a clusterfile with hydrogen atoms)\n"
        << "outfile is an outfile to append probable configs to\n"
        <<        "\t if not specified will be printed to probable.pdb \n\n";
        exit(0);
    }

    double s = 0;


    clock_t t;
    t = clock();
    int i = 0; string infile; string outfile;
    while (i<argc) {
        if (!strcmp(argv[i], "-i")) {
            infile = argv[++i];
        }
        if (!strcmp(argv[i], "-o")) {
            outfile = argv[++i];
        }
        i++;
    }

    if (infile.empty()) {
        cerr << "infile needs to be defined\n"
        << "For full run instructions run executable with no arguments\n\n";
        exit(0);
    }
    if (outfile.empty()) {
        cerr << "outfile needs to be defined\n"
        << "For full run instructions run executable with no arguments\n\n";
        exit(0);
    }


    FILE * pFile;
    pFile = fopen(outfile.c_str(), "a");
    //fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, vals[0], vals[1], vals[2], data[pos], T);



    vector<double > tmp;
    double temp;
    string strtemp;
    /*
    ifstream input(expfile.c_str());
    //getline(input, strtemp); //skip header
    while (!input.eof()) {
        getline(input, strtemp);
        if (!strtemp.empty()) {
            temp = atof(strtemp.substr(31, 7).c_str());
            tmp.push_back(temp);
            temp = atof(strtemp.substr(39, 7).c_str());
            tmp.push_back(temp);
            temp = atof(strtemp.substr(47, 7).c_str());
            tmp.push_back(temp);
        }
    }
    vector<double > tmp2;
    for (int i = 0; i < tmp.size(); i++) {
        if (i%9 == 0 || i%9==1 || i%9==2) {
            tmp2.push_back(tmp[i]);
        }
    }
     * /
    /*
    ofstream tout("test.dat"); tout.precision(16);
    tout << tmp2.size() << endl;
    for (int i = 0; i < tmp2.size(); i++) {
        tout << tmp2[i] << "\t";
        if (i%3== 2 && i!=0) {
            tout << endl;
        }
    }
    */
    //kdtree trans(tmp2);
    //cout << "made trans tree" << endl;

    vector<double > tmp4;
    ifstream stput(infile.c_str());
    getline(stput, strtemp); //skip header
    while (!stput.eof()) {
        getline(stput, strtemp);
        if (!strtemp.empty()) {
            temp = atof(strtemp.substr(31, 7).c_str());
            tmp4.push_back(temp);
            tmp.push_back(temp);
            temp = atof(strtemp.substr(39, 7).c_str());
            tmp4.push_back(temp);
            tmp.push_back(temp);
            temp = atof(strtemp.substr(47, 7).c_str());
            tmp4.push_back(temp);
            tmp.push_back(temp);
        }
    }
    vector<double > tmp5;
    for (int i = 0; i < tmp4.size(); i++) {
        if (i%9 == 0 || i%9 == 1 || i%9 == 2) {
            tmp5.push_back(tmp4[i]);
        }
    }

    kdtree trans(tmp5);
    int transi = 0; //index of closest trans
    int* indt;
    indt = new int[1];
    double* distt;
    distt = new double[1];
    double winner = 10000.00;
    for (i = 0; i < trans.npts; i++) {
        trans.nnearest(i, indt, distt, 1);
        if (distt[0] < winner) {
            winner = distt[0];
            transi = indt[0];
        }
    }

    delete distt;
    delete indt;
    //s = trans.run_tree_trans(tmp5);
    //transout << s << endl;
    //transout.close();
    /*
        Begin orientational code
    */
    vector<double > tmp3;
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

        singtest=((e[1]*e[2]) + (e[3]*e[0]));
        if (singtest > 0.4999) {
            tmp3.push_back(sin(pi/2));
            tmp3.push_back(0);
            tmp3.push_back(2*atan2(e[1],e[0]));
        }
        else if (singtest < -0.4999) {
            tmp3.push_back(sin(pi/-2));
            tmp3.push_back(0);
            tmp3.push_back(-2*atan2(e[1], e[0]));
        }
        else {
            tmp3.push_back(sin(asin(2*singtest)));
            tmp3.push_back(atan2(((2*e[1]*e[0])-(2*e[2]*e[3])) , (1 - (2*pow(e[1],2)) - (2*pow(e[3],2)))));
            tmp3.push_back(atan2(((2*e[2]*e[0])-(2*e[1]*e[3])) , (1 - (2*pow(e[2],2)) - (2*pow(e[3],2)))));
        }
    }
    kdtree orient(tmp3);
    //s = orient.run_tree_orient();
    //orientout << s << endl;
    //orientout.close();
    int orienti = 0; //index of closest orient
    int* indo;
    indo = new int[1];
    double* disto;
    disto = new double[1];
    winner = 10000.00;
    for (i = 0; i < orient.npts; i++) {
        orient.nnearest(i, indo, disto, 1);
        if (disto[0] < winner) {
            winner = disto[0];
            orienti = indo[0];
        }
    }

    delete disto;
    delete indo;

    /*
        Determined the best water orientation as orienti in array of pts
     *  Best oxygen position is the position of water at transi in pts
     *   need to find oxygen position of water with orienti in tmp array
     *      tmp array is for each water 9
     *        therefore position of oxygen is orienti * 9, orienti*9 + 1, and orienti*9 +2
     *   need to find distance of that water to transi
     *   translate
     */

    double orientO_x = tmp[orienti*9];
    double orientO_y = tmp[orienti*9 + 1];
    double orientO_z = tmp[orienti*9 + 2];
    double orientH1_x = tmp[orienti*9 + 3];
    double orientH1_y = tmp[orienti*9 + 4];
    double orientH1_z = tmp[orienti*9 + 5];
    double orientH2_x = tmp[orienti*9 + 6];
    double orientH2_y = tmp[orienti*9 + 7];
    double orientH2_z = tmp[orienti*9 + 8];

    string name = "ATOM";
    string atom = "O";
    string chainid = "X";
    int resseq = 1;
    string resname = "T3P";
    string testfile = "test.pdb";
    //FILE * tFile;
    //tFile = fopen(testfile.c_str(), "a");
    //fprintf (tFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), i, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, orientO_x, orientO_y, orientO_z, 0.0, 0.0);
    //fprintf (tFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), i, "H", resname.c_str(), chainid.c_str(), resseq, orientH1_x, orientH1_y, orientH1_z, 0.0, 0.0);
    //fprintf (tFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), i, "H", resname.c_str(), chainid.c_str(), resseq, orientH2_x, orientH2_y, orientH2_z, 0.0, 0.0);


    double transO_x = trans.pts[transi].x[0];
    double transO_y = trans.pts[transi].x[1];
    double transO_z = trans.pts[transi].x[2];

    double distx = (transO_x - orientO_x);
    double disty = (transO_y - orientO_y);
    double distz = (transO_z - orientO_z);

    orientO_x += distx;
    orientO_y += disty;
    orientO_z += distz;
    orientH1_x += distx;
    orientH1_y += disty;
    orientH1_z += distz;
    orientH2_x += distx;
    orientH2_y += disty;
    orientH2_z += distz;
    double ox = 0.0, oy = 0.0, oz = 0.0;
    for (i = 0; i < 3; i++) {
        if (i > 0) {atom = "H";}
        if (i == 0) {
            ox = orientO_x; oy = orientO_y; oz = orientO_z;
        }
        if (i == 1) {
            ox = orientH1_x; oy = orientH1_y; oz = orientH1_z;
        }
        if (i == 2){
            ox = orientH2_x; oy = orientH2_y; oz = orientH2_z;
        }
        fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), i, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, ox, oy, oz, 0.0, 0.0);
    }
}



