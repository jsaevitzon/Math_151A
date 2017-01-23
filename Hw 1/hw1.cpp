#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <vector>

#define pi 3.14159265359

using namespace std;

double absolute( double x ) {

  if( x < 0 )
    return -x;
  return x;
}

double f(float x);

double f_prime(float x) {

  float h = .000001;

  return (f(x+h) - f(x))/h;

}

// f(x) = sin(x) - x = 0 on [-π/4, π/2]
double f(float x){

  return sin(x) - x;
  
}


void build_graph_file(vector<double> seq, string func){
  
  ofstream file(func);
  file << "#x y" << endl;
  
  for( int i = 0; i < seq.size(); i++ ){
    file << i+1 << ' ' << seq[i] << endl;
    cout << i+1 << ' ' << seq[i] << endl;
  }

  file.close();

}

double Newton( double p_0, double TOL, int N) {

  //vector x stores the sequence
  //for graphing purposes
  vector<double> seq;
  seq.push_back(p_0);

  //f_p is f(p) and f_pp is f'(p)  
  int i = 1;
  double  f_p = f(p_0), f_pp =  f_prime(p_0);
  double  p_1;

  while ( i <= N ) {
    
    
    //get current p value
    p_1 = p_0 - f_p/f_pp;
    seq.push_back(p_1);
    
    //test error
    if ( absolute(p_1 - p_0) < TOL ) {
      build_graph_file(seq, "graph_newt.dat");
      return p_1;
    }
    //increment iteration, update p_0, f_p, f_prime
    i++;
    p_0  = p_1;
    f_p  = f(p_0);
    f_pp = f_prime(p_0); 

  }

  cout << "Too many iterations have passed, no solution found";
  return -1;


}

double Secant( double p_0, double p_1, double TOL, int N ) {
  
  //vector x stores the sequence
  //for graphing purposesk
  vector<double> seq;
  seq.push_back(p_0);
  seq.push_back(p_1);
  //f_p is f(p) and f_pp is f'(p)                                                                                                                                                                                                           
  int i = 1;
  double  f_p0 = f(p_0), f_p1 =  f(p_1);
  double  p_2;

  while ( i <= N ) {
    
    p_2 = p_1 - (f_p1*(p_1-p_0))/(f_p1-f_p0);
    seq.push_back(p_2);

    if ( absolute(p_2 - p_0) < TOL ) {
      build_graph_file(seq, "graph_secant.dat");
      return p_2;
    }
    
    //increment iteration, update p_0, f_p, f_prime
    i++;
    p_0   = p_1;
    p_1   = p_2;
    f_p0  = f_p1;
    f_p1  = f(p_1);

  }
  
  cout << "Too many iterations have passed, no solution found";
  return -1;
}

int main()
{

  double TOL = pow(10, -7), p_0 = pi/4, p_1 = pi*3.0/8.0;
  int N = 100;

  cout << Newton(p_0, TOL, N);
  

  cout << Secant(p_0, p_1, TOL, N);


  
  

}
