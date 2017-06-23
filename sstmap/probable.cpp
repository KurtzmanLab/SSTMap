#include "probable.h"

using namespace std;

point::point () {
    dim = D;
    x = new double[dim];
    for (int i = 0; i < dim; i++) {
        x[i] = 0;
    }

}

void point::set_point(double* vals) {
    //dim = d;
    //x = new double[dim];
    for (int i = 0; i < dim; i++) {
        x = &vals[i];
    }
}
void point::set_point(const point &p) {
    dim = p.dim;
    //x = new double[dim];
    for (int i = 0; i < dim; i++) {
        x[i] = p.x[i];
    }
}

/*point::~point () {
    delete x;
}*/

void point::zeros() {
    for (int i = 0; i < dim; i++) {
        x[i] = 0;
    }
}

void point::ones() {
    for (int i = 0; i < dim; i++) {
        x[i] = 1;
    }
}

void point::print_point() {
    for (int i = 0; i < dim; i++) {
        cout << x[i] << "\t";
    }
    cout << endl;
}


double dist(const point &p, const point &q) {
    if (p.dim != q.dim) {
        cerr << "Dimensions of points do not match in distance comparison!!\n";
        exit(EXIT_FAILURE);
    }
    double distance = 0.0;
    for (int i = 0; i < p.dim; i++) {
        distance += pow((p.x[i] - q.x[i]), 2);
    }
    double large = 10000;
    if (distance == 0) return large;
    return sqrt(distance);

}

bool same_points (const point &p, const point &q) {
    if (p.dim != q.dim) {
        return false;
    }
    else {
        for (int i = 0; i < p.dim; i++) {
            if (p.x[i] != q.x[i]) {
                return false;
            }
        }
    }
    return true;
}

boxnode::boxnode() {
    mom = 0; dau1 = 0; dau2 = 0; pthi = 0; ptlo = 0;
    //By default the points will be set to zero respectively
}

void boxnode::set_boxnode(point mylo, point myhi, int mymom, int myd1, int myd2, int myptlo, int mypthi) {
    //mybox.set_box(mylo, myhi);
    hi.set_point(myhi);
    lo.set_point(mylo);
    //cout << "set box\n";
    mom = mymom;
    dau1 = myd1;
    dau2 = myd2;
    ptlo = myptlo;
    pthi = mypthi;
    //cout << "done set box\n";
}

double dist(const boxnode &b, const point &p) {
    double distance = 0.0;
    if (p.dim != b.lo.dim || p.dim != b.hi.dim) {
        cerr << "Point and Box Points do not have the same dimensionality in distance calculation!!\n";
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < p.dim; i++) {
        if (p.x[i] < b.lo.x[i]) distance += pow((p.x[i]-b.lo.x[i]), 2);
        if (p.x[i] > b.hi.x[i]) distance += pow((p.x[i]-b.hi.x[i]), 2);
    }
    return sqrt(distance);
    //This will return 0 if the point is in the box
}

int selecti(const int k, int *indx, int n, double *arr) {
    int i, ia, ir, j, l, mid;
    double a;

    l = 0;

    ir = n-1;
    for (;;) {
        if (ir <= l+1) {
            if (ir == l+1 && arr[indx[ir]] < arr[indx[l]]) {
                swap(indx[l], indx[ir]);
            }
            return indx[k]; //final end point
        }
        else {
            mid = (l+ir) >> 1;
            swap(indx[mid], indx[l+1]);
            if (arr[indx[l]] > arr[indx[ir]]) swap(indx[l], indx[ir]);
            if (arr[indx[l+1]] > arr[indx[ir]]) swap(indx[l+1], indx[ir]);
            if (arr[indx[l]] > arr[indx[l+1]]) swap(indx[l], indx[l+1]);
            i = l+1;
            j = ir;
            ia = indx[l+1];
            a = arr[ia];
            for (;;) {
                do i++; while (arr[indx[i]] < a);
                do j--; while (arr[indx[j]] > a);
                if (j < i) break; //inner endpoint
                swap(indx[i], indx[j]);

            }

            indx[l+1] = indx[j];
            indx[j] = ia;
            if (j >= k) ir=j-1;
            if (j <= k) l = i;
        }
    }

}

const double kdtree::BIG(1.0e99);

kdtree::kdtree(std::vector< double > &vals) {
    /*
        This function assumes the doubles fed in through vals are only the pertinent ones. IE If this is 3d its the positions or orientations and nothing else.
    */
    /*for (int i = 0; i < 3; i++) {
        minn[i] = 10000;
        maxx[i] = 0;
    }*/
    //BIG = 1.0e99;
    dim = D;
    nd = new int[1];
    dn = new double[1];
    //within1 = new int[3000];
    npts = vals.size()/D;
    //cout << npts << endl;
    pts = new point[npts];
    /*
    cenn[0] = x;
    cenn[1] = y;
    cenn[2] = z;
    */
    int foo = 0;
    for (int i = 0; i < vals.size(); i++) {
        if (i%D==0) {
            pts[foo].x[0] = vals[i];
            //minn[0] = min(minn[0], vals[i]);
            //maxx[0] = max(maxx[0], vals[i]);
        }
        if (i%D==1) {
            pts[foo].x[1] = vals[i];
            //minn[1] = min(minn[1], vals[i]);
            //maxx[1] = max(maxx[1], vals[i]);
        }
        if (i%D==2) {
            pts[foo].x[2] = vals[i];
            //minn[2] = min(minn[2], vals[i]);
            //maxx[2] = max(maxx[2], vals[i]);
            foo++;
        }
    }

    /*for (int i = 0; i < 3; i++) {
        cenn[i] = (maxx[i]+minn[i])/2;
    }*/

    ptindx = new int[npts]; rptindx = new int[npts];
    int ntmp, m, kk, k, j, nowtask, jbox, np, tmom, tdim, ptlo, pthi;
    int *hp;
    double *cp;
    int taskmom[50], taskdim[50];
    for (k = 0; k < npts; k++) ptindx[k] = k;
    m = 1;
    for (ntmp = npts; ntmp; ntmp >>= 1) {
        m <<= 1;
    }
    //cout << "npts: " << npts << endl;
    numbox = 2*npts - (m>>1);
    if (m < numbox) numbox = m;
    //m = pow(2, (log(npts)/log(2)));
    //numbox = 2*npts - m/2;
    if (m < numbox) numbox = m;
    numbox--;
    //cout << "about to make boxes\n";
    //cout << numbox << endl;
    boxes = new boxnode[numbox];
    //cout << "made boxes\n";
    coord = new double[D*npts];
    //cout << "made coords\n";
    for (j = 0, kk = 0; j < D; j++, kk+= npts) {
        for (k = 0; k < npts; k++) {
            //cout << k << endl;
            coord[kk+k] = pts[k].x[j];
        }
    }
    //cout << "set coords\n";
    //double* nums;
    //double* nnums;
    point lo, hi;
    for (int i = 0; i < 3; i++) {
        hi.x[i] = BIG;
        lo.x[i] = -BIG;
    }

    //cout << hi.x[0] << "\t" << hi.x[1] << "\t" << hi.x[2] << endl;
    //cout << lo.x[0] << "\t" << lo.x[1] << "\t" << lo.x[2] << endl;
    boxes[0].set_boxnode(lo, hi, 0, 0, 0, 0, npts-1);
    /*if (D == 3) {
        //cout << "start D3 if\n";
        nums = new double[3];
        nnums = new double[3];
        for (int i =0; i <3; i++) {
            nums[i] = BIG;
            nnums[i] = -BIG;
        }
        cout << nums[0] << "\t" << nums[1] << "\t" << nums[2] << endl;
        cout << nnums[0] << "\t" << nnums[1] << "\t" << nnums[2] << endl;
        lo.set_point(nnums, D);
        hi.set_point(nums, D);
        //cout << "made BIG points\n";
        cout << hi.x[0] << "\t" << hi.x[1] << "\t" << hi.x[2] << endl;
        cout << lo.x[0] << "\t" << lo.x[1] << "\t" << lo.x[2] << endl;
        boxes[0].set_boxnode(lo, hi, 0, 0, 0, 0, npts-1);
        //cout << "made first box\n";
    }
    if (D == 6) {
        //cout << "start D6 if\n";
        nums = new double[6];
        nnums = new double[6];
        for (int i = 0; i < 6; i++) {
            nums[i] = BIG;
            nnums[i] = -BIG;
        }
        lo.set_point(nnums, D), hi.set_point(nums, D);
        boxes[0].set_boxnode(lo, hi, 0, 0, 0, 0, npts-1);
    }
    delete nums;
    delete nnums;
    */
    //cout << "Set initial box: \n\n";
    //cout << boxes[0].hi.x[0] << "\t" << boxes[0].hi.x[1] << "\t" << boxes[0].hi.x[2] << endl;
    //cout << boxes[0].lo.x[0] << "\t" << boxes[0].lo.x[1] << "\t" << boxes[0].lo.x[2] << endl;

    for (int i = 0; i < 3; i++) {
        boxes[0].hi.x[i] = BIG;
        boxes[0].lo.x[i] = -BIG;
    }

    //cout << "Fix initial box: \n\n";
    //cout << boxes[0].hi.x[0] << "\t" << boxes[0].hi.x[1] << "\t" << boxes[0].hi.x[2] << endl;
    //cout << boxes[0].lo.x[0] << "\t" << boxes[0].lo.x[1] << "\t" << boxes[0].lo.x[2] << endl;

    jbox = 0;
    taskmom[1] = 0;
    taskdim[1] = 0;
    nowtask = 1;
    //cout << "got to while loop\n";
    while (nowtask) {
        tmom = taskmom[nowtask];
        tdim = taskdim[nowtask--];
        ptlo = boxes[tmom].ptlo;
        pthi = boxes[tmom].pthi;
        hp =  &ptindx[ptlo];
        cp = &coord[tdim*npts];
        np = pthi - ptlo + 1;
        kk = (np-1)/2;
        selecti(kk, hp, np, cp);
        hi = boxes[tmom].hi;
        lo = boxes[tmom].lo;
        //hi.x[tdim] = lo.x[tdim] = coord[tdim*npts + hp[kk]];
        //cout << jbox << endl;
        boxes[++jbox].set_boxnode(boxes[tmom].lo, hi, tmom, 0, 0, ptlo, ptlo+kk);
        boxes[jbox].hi.x[tdim] = coord[tdim*npts + hp[kk]];
        //cout << jbox << endl;
        boxes[++jbox].set_boxnode(lo, boxes[tmom].hi, tmom, 0 , 0, ptlo+kk+1, pthi);
        boxes[jbox].lo.x[tdim] = coord[tdim*npts + hp[kk]];
        boxes[tmom].dau1 = jbox-1;
        boxes[tmom].dau2 = jbox;
        if (kk > 1) {
            taskmom[++nowtask] = jbox-1;
            taskdim[nowtask] = (tdim+1)%D;
        }
        if (np - kk > 3) {
            taskmom[++nowtask] = jbox;
            taskdim[nowtask] = (tdim+1)%D;
        }
    }
    for (j = 0; j < npts; j++) rptindx[ptindx[j]] = j;
    //cout << "made tree" << endl;

    //cout << "delete coord" << endl;
}


kdtree::~kdtree () {
    delete boxes;
    delete ptindx;
    delete rptindx;
    delete dn;
    delete nd;
    delete coord;
    //delete within1;
}

double kdtree::disti(int jpt, int kpt) {
    if (jpt == kpt) return BIG; //to avoid the closest neighbor is itself
    else return dist(pts[jpt], pts[kpt]);
}

int kdtree::locate(point pt) {
    int nb, d1, jdim;
    nb = jdim = 0;
    while (boxes[nb].dau1) { //basically keep going until bottom from root
        d1 = boxes[nb].dau1;
        if (pt.x[jdim] <= boxes[d1].hi.x[jdim]) nb = d1;
        else nb = boxes[nb].dau2;
        jdim = ++jdim%D;
    }
    return nb;
}

int kdtree::locate(int jpt) {
    int nb, d1, jh;
    jh = rptindx[jpt];
    nb = 0;
    while (boxes[nb].dau1) {
        d1 = boxes[nb].dau1;
        if (jh <= boxes[d1].pthi) nb = d1;
        else nb = boxes[nb].dau2;
    }
    return nb;
}

double kdtree::dnearest(point pt) {
    int i, k, nrst, ntask;
    int task[50];
    double dnrst = BIG, d;
    k = locate(pt);
    for (i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
        d = dist(pts[ptindx[i]], pt);
        if (d < dnrst && d != 0) {
            //this fix for != 0 may result in some uncertain behavior, will be necessary to check this out.
            nrst = ptindx[i];
            dnrst = d;
        }
    }
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (dist(boxes[k], pt) < dnrst) {
            if (boxes[k].dau1) {
                task[++ntask] = boxes[k].dau1;
                task[++ntask] = boxes[k].dau2;
            }
            else {
                for (i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
                    d = dist(pts[ptindx[i]], pt);
                    if (d < dnrst && d != 0) {
                        nrst = ptindx[i];
                        dnrst = d;
                    }
                }
            }
        }
    }
    return dnrst;
}

void kdtree::nnearest(int jpt, int* nn, double* dn, int n) {
    int i, k, ntask, kp;
    int task[50];
    double d;
    if (n > npts-1) throw("you're asking for too much buddy (nn > npts)");
    for (i = 0; i < n; i++) dn[i] = BIG;
    kp = boxes[locate(jpt)].mom;
    while (boxes[kp].pthi - boxes[kp].ptlo < n) kp = boxes[kp].mom;
    for (i = boxes[kp].ptlo; i <= boxes[kp].pthi; i++) {
        if (jpt == ptindx[i]) continue;
        d = disti(ptindx[i], jpt);
        if (d < dn[0]) {
            dn[0] = d;
            nn[0] = ptindx[i];
            if (n>1) sift_down(dn, nn, n);
        }
    }
    task[1] = 0;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (k == kp) continue;
        if (dist(boxes[k], pts[jpt]) < dn[0]) {
            if (boxes[k].dau1) {
                task[++ntask] = boxes[k].dau1;
                task[++ntask] = boxes[k].dau2;
            }
            else {
                for (i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
                    d = disti(ptindx[i], jpt);
                    if (d < dn[0]) {
                        dn[0] = d;
                        nn[0] = ptindx[i];
                        if (n > 1) sift_down(dn, nn, n);
                    }
                }
            }
        }
    }
    return;
}

void kdtree::sift_down(double* heap, int* ndx, int nn) {
    int n = nn - 1;
    int j, jold, ia;
    double a;
    a = heap[0];
    ia = ndx[0];
    jold = 0;
    j = 1;
    while (j <= n) {
        if (j < n && heap[j] < heap[j+1]) j++;
        if (a >= heap[j]) break;
        heap[jold] = heap[j];
        ndx[jold] = ndx[j];
        jold = j;
        j = 2*j+1;
    }
    heap[jold] = a;
    ndx[jold] = ia;
}


int kdtree::locatenear(point pt, double r, int *v, int nmax) {
    /*
        This fuction returns all the points within some distance of a target point. I dont think we will ever use it.
    */
    int k, i, nb, nbold, nret, ntask, jdim, d1, d2;
    int task[50];
    nb = jdim = nret = 0;
    if (r < 0.0) throw("radius must be nonnegative");
    while (boxes[nb].dau1) {
        nbold = nb;
        d1 = boxes[nb].dau1;
        d2 = boxes[nb].dau2;
        if (pt.x[jdim] + r <= boxes[d1].hi.x[jdim]) nb = d1;
        else if (pt.x[jdim] - r >= boxes[d2].lo.x[jdim]) nb = d2;
        jdim = ++jdim%D;
        if (nb == nbold) break;
    }
    //cout << nb << endl;
    task[1] = nb;
    ntask = 1;
    while (ntask) {
        k = task[ntask--];
        if (dist(boxes[k], pt) > r) {
            //cout << "box out of range: " << dist(boxes[k], pt) << endl;
            //cout << boxes[k].hi.x[0] << "\t" << boxes[k].hi.x[1] << "\t" << boxes[k].hi.x[2] << endl;
            continue;
        }
        else {
            //cout << "box in range\n";
        }
        if (boxes[k].dau1) {
            task[++ntask] = boxes[k].dau1;
            task[++ntask] = boxes[k].dau2;
        }
        else {
            for (i = boxes[k].ptlo; i <= boxes[k].pthi; i++) {
                if (dist(pts[ptindx[i]], pt) <= r && nret < nmax) {
                    v[nret++] = ptindx[i];
                }
                if (nret == nmax) return nmax;
            }
        }
    }
    return nret;
}
/*
double kdtree::run_tree() {
    //ofstream output;
    //output.open(OUT.c_str());
    //output.precision(16);
    //cout << "run tree start" << endl;
    double gd = 0;
    double s = 0.0;
    double T = 300.;
    double R = 8.314472;
    double pi = 3.14159265359;
    double cenndist = 10000;
    int npts2 = 0;
    int fcount = 10000;
    for (int i = 0; i < npts; i++) {
        //if (pts[i].x[0] > maxx[0]-5 || pts[i].x[0] < minn[0]+5 || pts[i].x[1] > maxx[1]-5 || pts[i].x[1] < minn[1]+5 || pts[i].x[2] > maxx[2]-5 || pts[i].x[2] < minn[2]+5) {
        //    continue;
        //}
        //else {
        cenndist = pow((pts[i].x[0] - cenn[0]), 2) + pow((pts[i].x[1] - cenn[1]), 2) + pow((pts[i].x[2] - cenn[2]), 2);
        //if (abs(pts[i].x[0] - cenn[0]) < 1 && abs(pts[i].x[1] - cenn[1]) < 1 && abs(pts[i].x[2] - cenn[2]) < 1) {
        if (cenndist <= 1) {
            nnearest(i, nd, dn, 1);
            //cout << dn[0] << endl;
            gd += log((0.0329223149*fcount*4*pi*pow(dn[0], 3))/3);
            npts2++;
        }
    }
    //cout << "\n\n";
    //cout << gd << endl;
    //cout << npts2 << endl;
    s = R*T*0.239*(gd/npts2 + 0.5772156649)/1000;
    //cout << s << endl;
    return s;
}
*/
double kdtree::run_tree_trans(std::vector<double > &cls) {
    point pt;
    double* dh;
    int numvals = cls.size()/3;
    dh = new double[numvals];
    int vecpos = 0;
    for (int i = 0; i < cls.size(); i+=3) {
        //run through the vector of the acknowledged standard cluster file, skipping hydrogens until the end
        pt.x[0] = cls[i]; pt.x[1] = cls[i+1]; pt.x[2] = cls[i+2];
        dh[vecpos] = dnearest(pt);
        vecpos++;
    }

    double gd = 0;
    double s = 0.0;
    double T = 300.;
    double R = 8.314472;
    double pi = 3.14159265359;

    int fcount = 10000;

    for (int i = 0; i < numvals; i++) {
        gd += log((0.0329223149*fcount*4*pi*pow(dh[i], 3))/3);
    }

    s = R*T*0.239*(gd/numvals + 0.5772156649)/1000;

    delete dh;
    return s;
}

double kdtree::run_tree_orient() {
    //ofstream ori("orientdists.txt"); ori.precision(16);
    double gd = 0;
    double s = 0.0;
    double T = 300.;
    double R = 8.314472;
    double pi = 3.14159265359;
    double de = 10000;
    point z;
    for (int i = 0; i < npts; i++) {
        nnearest(i, nd, dn, 1);
        //in order to implement this need to use a function which takes a point and returns a distance.
        //dnearest from before.
        if (pts[i].x[0] > pi/2) {
            z.x[0] = pts[i].x[0] - 2*pi;
            z.x[1] = pts[i].x[1];
            z.x[2] = pts[i].x[2];
            de = dnearest(z);
            if (de < dn[0] && de != 0) {
                dn[0] = de;
            }
        }
        else if (pts[i].x[0] < -pi/2) {
            z.x[0] = pts[i].x[0] + 2*pi;
            z.x[1] = pts[i].x[1];
            z.x[2] = pts[i].x[2];
            de = dnearest(z);
            if (de < dn[0] && de != 0) {
                dn[0] = de;
            }
        }
        else if (pts[i].x[1] > pi/2) {
            z.x[0] = pts[i].x[0];
            z.x[1] = pts[i].x[1] - 2*pi;
            z.x[2] = pts[i].x[2];
            de = dnearest(z);
            if (de < dn[0] && de != 0) {
                dn[0] = de;
            }
        }
        else if (pts[i].x[1] < -pi/2) {
            z.x[0] = pts[i].x[0];
            z.x[1] = pts[i].x[1] + 2*pi;
            z.x[2] = pts[i].x[2];
            de = dnearest(z);
            if (de < dn[0] && de != 0) {
                dn[0] = de;
            }
        }
        else if (pts[i].x[2] > pi/2) {
            z.x[0] = pts[i].x[0];
            z.x[1] = pts[i].x[1];
            z.x[2] = pts[i].x[2] - 2*pi;
            de = dnearest(z);
            if (de < dn[0] && de != 0) {
                dn[0] = de;
            }
        }
        else if (pts[i].x[2] < -pi/2) {
            z.x[0] = pts[i].x[0];
            z.x[1] = pts[i].x[1];
            z.x[2] = pts[i].x[2] + 2*pi;
            de = dnearest(z);
            if (de < dn[0] && de != 0) {
                dn[0] = de;
            }
        }
        //ori << dn[0] << endl;
        gd += log((pow(dn[0], 3)*npts)/(6*pi));
    }
    //cout << gd << endl;
    //cout << npts << endl;
    s = R*T*0.239*(gd/npts + 0.5772156649)/1000;
    //cout << s << endl;
    return s;
}

/*void kdtree::run_locate() {
    double dn = BIG;
    double d;
    int nr = 0;
    for (int i = 0; i < npts; i++) {
        nr = locatenear(pts[i], 0.25, within1, 3000);
        //cout << nr << endl;
        for (int j = 0; j < nr; j++) {
            d = disti(i, within1[j]);
            if (d < dn) {
                dn = d;
            }
        }
        //cout << dn << endl;
    }
}*/


