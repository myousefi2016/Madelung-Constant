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

const unsigned int nx=100;
const unsigned int ny=100;
const unsigned int nz=100;
const unsigned int n=8.1;
const double r0=2.82;
const double epsilon=1e-1;

double Madelung_constant1() {

int h,k,l;
double sum;
sum=0.0;

 for (h=0; h<nx+1; h++) {
  for (k=0; k<ny+1; k++) {
   for (l=0; l<nz+1; l++) {
 
 if (h!=0 && k!=0 && l!=0) {
 sum += pow(-1.0,h+k+l)/sqrt(pow(h,2)+pow(k,2)+pow(l,2));
 }


 }
 }
 }

 return sum;

}

void Madelung_constant1_convergence(FILE* stream) {

int h,k,l;
double sum;
sum=0.0;

for (int x=3; x<201; x++) {

 for (h=0; h<x+1; h++) {
  for (k=0; k<x+1; k++) {
   for (l=0; l<x+1; l++) {
 
 if (h!=0 && k!=0 && l!=0) {
 sum += pow(-1.0,h+k+l)/sqrt(pow(h,2)+pow(k,2)+pow(l,2));
 }


 }
 }
 }

 fprintf(stream, "%d\t%f\n", x, sum);

 sum=0.0;
 
 }

}

double Madelung_constant2() {

int h,k,l;
double sum;
sum=0.0;

 for (h=0; h<nx+1; h++) {
  for (k=0; k<ny+1; k++) {
   for (l=0; l<nz+1; l++) {
 
 if (h!=0 && k!=0 && l!=0) {
 sum += 1.0/sqrt(pow(pow(h,2)+pow(k,2)+pow(l,2),n));
 }


 }
 }
 }

 return sum;

}

double Madelung_constant2_convergence(FILE* stream) {

int h,k,l;
double sum;
sum=0.0;

for (int x=3; x<201; x++) {

 for (h=0; h<x+1; h++) {
  for (k=0; k<x+1; k++) {
   for (l=0; l<x+1; l++) {
 
 if (h!=0 && k!=0 && l!=0) {
 sum += 1.0/sqrt(pow(pow(h,2)+pow(k,2)+pow(l,2),n));
 }


 }
 }
 }

 fprintf(stream, "%d\t%f\n", x, sum);

 sum=0.0; 

 }

}

void Energy(vector<vector<vector<double>>> &U, double alpha, double M1, double M2) {

double r;

for (int i=0; i<=nx; i++) {
for (int j=0; j<=ny; j++) {
for (int k=0; k<=nz; k++) {

if (i==0 && j==0 && k==0) {
r = r0*epsilon;
} else {
r = r0*sqrt(pow(i,2)+pow(j,2)+pow(k,2));
}

U[i][j][k]=((-M1)/r)+(alpha*M2)/pow(r,n);

}
}
}

}

void write_output_vtk(vector<vector<vector<double>>>& U, int nx, int ny, int nz)
{
    string name = "./output.vtk";
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

FILE* f_convergence;
char filename[40];
double M1, M2;
double alpha = 684;
vector<vector<vector<double>>> U       (nx+1,vector<vector<double>>(ny+1,vector<double>(nz+1,0)));

M1 = Madelung_constant1();

M2 = Madelung_constant2();

sprintf(filename, "./Madelung1.txt");

f_convergence = fopen(filename, "w");

Madelung_constant1_convergence(f_convergence);

fclose(f_convergence);

sprintf(filename, "./Madelung2.txt");

f_convergence = fopen(filename, "w");

Madelung_constant2_convergence(f_convergence);

fclose(f_convergence);

Energy(U,alpha,M1,M2);

write_output_vtk(U,nx+1,ny+1,nz+1);

printf("The first Madelung constant is = %f and the second one is = %f\n",M1,M2); 

}
