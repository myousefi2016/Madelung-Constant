#include <stdio.h>
#include <stdlib.h> // dynamic allocation
#include <math.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

const unsigned int nx=20;
const unsigned int ny=20;
const unsigned int nz=20;
const unsigned int n=8.1;
const double r0=2.82;
const double epsilon=1e-1;
const double beta=0.5;

void Madelung_constant1(vector<vector<vector<double>>> &M1) {

int h,k,l;

for (int xc=0; xc<=nx; xc += 4) {
for (int yc=0; yc<=ny; yc += 4) {
for (int zc=0; zc<=nz; zc += 4) {

 for (h=0; h<nx+1; h += 4) {
  for (k=0; k<ny+1; k += 4) {
   for (l=0; l<nz+1; l += 4) {
 
 if (h!=xc && k!=yc && l!=zc) {
 M1[xc][yc][zc] += pow(-1.0,h+k+l-xc-yc-zc)/(sqrt(pow(h-xc,2)+pow(k-yc,2)+pow(l-zc,2))*exp(beta*r0*sqrt(pow(h-xc,2)+pow(k-yc,2)+pow(l-zc,2))));
 }


 }
 }
 }

}
}
}

}

void Madelung_constant2(vector<vector<vector<double>>> &M2) {

int h,k,l;

for (int xc=0; xc<=nx; xc += 4) {
for (int yc=0; yc<=ny; yc += 4) {
for (int zc=0; zc<=nz; zc += 4) {

 for (h=0; h<nx+1; h += 2) {
  for (k=0; k<ny+1; k += 2) {
   for (l=0; l<nz+1; l += 2) {
 
 if (h!=xc && k!=yc && l!=zc) {
 M2[xc][yc][zc] += 1.0/sqrt(pow(pow(h-xc,2)+pow(k-yc,2)+pow(l-zc,2),n));
 }


 }
 }
 }

}
}
}

}

void Energy(vector<vector<vector<vector<vector<vector<double>>>>>> &U, double alpha, vector<vector<vector<double>>> &M1, vector<vector<vector<double>>> &M2) {

double r;

for (int xc=0; xc<=nx; xc += 4) {
for (int yc=0; yc<=ny; yc += 4) {
for (int zc=0; zc<=nz; zc += 4) {
for (int i=0; i<=nx; i++) {
for (int j=0; j<=ny; j++) {
for (int k=0; k<=nz; k++) {

if (i==xc && j==yc && k==zc) {
r = r0*epsilon;
} else {
r = r0*sqrt(pow(i-xc,2)+pow(j-yc,2)+pow(k-zc,2));
}

U[i][j][k][xc][yc][zc]=((-M1[xc][yc][zc])/r)+(alpha*M2[xc][yc][zc])/pow(r,n);

}
}
}
}
}
}

}


void TotalEnergy(vector<vector<vector<vector<vector<vector<double>>>>>> &U, vector<vector<vector<double>>> &Et) {

for (int i=0; i<=nx; i++) {
for (int j=0; j<=ny; j++) {
for (int k=0; k<=nz; k++) {
for (int xc=0; xc<=nx; xc += 4) {
for (int yc=0; yc<=ny; yc += 4) {
for (int zc=0; zc<=nz; zc += 4) {

Et[i][j][k] += U[i][j][k][xc][yc][zc];

}
}
}
}
}
}

}

void write_output_vtk(vector<vector<vector<double>>>& U, int nx, int ny, int nz)
{
    string name = "./outputTotalEnergymodified.vtk";
    ofstream ofile (name);

    // vtk preamble
    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << "OUTPUT by LIBM\n";
    ofile << "ASCII" << endl;

    // write grid
    ofile << "DATASET RECTILINEAR_GRID" << endl;
    ofile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    ofile << "X_COORDINATES " << nx << " float" << endl;
    for(size_t i = 0; i < nx; i++)
        ofile << i << "\t";
    ofile << endl;
    ofile << "Y_COORDINATES " << ny << " float" << endl;
    for(size_t i = 0; i < ny; i++)
        ofile << i << "\t";
    ofile << endl;
    ofile << "Z_COORDINATES " << nz << " float" << endl;
    for(size_t i = 0; i < nz; i++)
        ofile << i << "\t";
    ofile << endl;

    // point data
    ofile << "POINT_DATA " << nx*ny*nz << endl;

    // write rho
    ofile << "SCALARS " << "U" << " double" << endl;
    ofile << "LOOKUP_TABLE default" << endl;
  for (int k = 0; k < nz; k++) 
    for(int j = 0; j < ny; j++)
        for(int i = 0; i < nx; i++)
            ofile << U[i][j][k] << endl;


}


int main(int argc, char** argv) {

double alpha = 684;
vector<vector<vector<double>>> Et       (nx+1,vector<vector<double>>(ny+1,vector<double>(nz+1,0)));
vector<vector<vector<double>>> M1       (nx+1,vector<vector<double>>(ny+1,vector<double>(nz+1,0)));
vector<vector<vector<double>>> M2       (nx+1,vector<vector<double>>(ny+1,vector<double>(nz+1,0)));
vector<vector<vector<vector<vector<vector<double>>>>>> U       (nx+1,vector<vector<vector<vector<vector<double>>>>>(ny+1,vector<vector<vector<vector<double>>>>(nz+1,vector<vector<vector<double>>>(nx+1,vector<vector<double>>(ny+1,vector<double>(nz+1,0))))));

Madelung_constant1(M1);

Madelung_constant2(M2);

Energy(U,alpha,M1,M2);

TotalEnergy(U,Et);

write_output_vtk(Et,nx+1,ny+1,nz+1); 

}
