#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

/*
Matthew Golden- May 2018

This code simulates an orbit around a Schwarzschild black hole.
We use Schwarschild coordinates so that spheres have surface area 4*pi*r^2

ds^2 = -f(r) dt^2 + f(r)^-1 dr^2 + r^2 dOmega^2

f(r) = 1 - 2M/r

We take G=c=M=1 so f(r) = 1 - 2/r
*/

/*
The characteristic function of a black hole
We assume WLOG M=1
double r - the coordinate radius of the spacetime
*/
double f( double r)
{
    return 1. - 2./r;
}

/*
The derivative of f(r) with respect to r.
This function is needed for many Christoffel symbols
double r - the coordinate radius of the spacetime 
*/
double df( double r )
{
    return 2./(r*r);
}

/*
This calculates the time derivative of all system variables
double* td     - a pointer to store these values at.
                 allocate this before using. Need 6 doubles
double* system - a pointer to doubles (t,r,theta,dt,dr,dtheta) 
                 dvar is time deriv of var
*/
void calc_td( double* td, double* system)
{
    //First three are trivial
    td[0] = system[3];
    td[1] = system[4];
    td[2] = system[5];
    //Use christoffel symbols
    td[3] = -df(system[1])/f(system[1])*system[3]*system[4];
    td[4] = -0.5*df(system[1])*f(system[1])*system[3]*system[3]
            +0.5*df(system[1])/f(system[1])*system[4]*system[4]
            +f(system[1])*system[1]*system[5]*system[5];
    td[5] = -(2./system[1])*system[4]*system[5];
    
    return;
}

/*
Writes the system to a file for processing.
ofstream* data - The file to pass the data to. The function needs a pointer
                 to this object, it cannot call by value
double* system - a pointer to the system variables
double   tau   - the proper time of the particle
*/
void write_system( ofstream* data, double* system, double tau )
{
    *data << scientific << tau << "\t" << system[1]*cos(system[2]) 
         << "\t" << system[1]*sin(system[2]) << "\n";
}

/*
This reads in parameters for the simulation so that you don't have to recompile every time you change initial conditions.
double* system  - pointer to the physical variables.
double* max_tau - a pointer to the end proper time of the simulation 
double* dtau    - the timestep of the simulation
int*    every   - the script will write to file every "every" cycles. 
                  e.g. if *every = 100, every 100 timesteps the system will be recorded
*/
void read_parameters( double* system, double* max_tau, double* dtau, int* every)
{
    ifstream parameters;
    parameters.open("parameters.dat");
    string junk;     //place to dump text from file
    int use_circ; //If 1, fixes l
    double l;
    system[0] = 0.; //time  always starts at 0
    system[2] = 0.; //theta always starts at 0
    parameters >> junk >> system[1]
               >> junk >> system[4]
               >> junk >> use_circ
               >> junk >> system[5]
               >> junk >> *max_tau
               >> junk >> *dtau
               >> junk >> *every;
    if(use_circ == 1)
    {
        l         = sqrt(system[1]*system[1]/(system[1]-3));  
        system[5] = l/(system[1]*system[1]);
    }

    //Normalize the four-vector
    system[3] = sqrt( (  system[4]*system[4]/f(system[1]) 
                       + system[1]*system[1]*system[5]*system[5]+1
                       )/f(system[1]) );


    parameters.close();

    return;
}

int main( int argc, char* argv[])
{
    /*
    system is a pointer to the following variables:
    system[0] = t
    system[1] = r
    system[2] = theta
    system[3] = dt/dtau
    system[4] = dr/dtau
    system[5] = dtheta/dtau
    I am using third order Runge-Kutta to evolve the system, 
    so we need two extra copies. 
    td is the time derivaitve of the system
    */
    double system[6];
    double system1[6];
    double system2[6];
    double td[6];
    double td1[6];
    double td2[6];


    //Additional parameters to be read
    double max_tau;
    double dtau;
    int    every;
    read_parameters( system, &max_tau, &dtau, &every ); 
 
    ofstream data;
    data.open("orbit.dat");

    //Evolution loop
    double tau     = 0;
    int    counter = 0;
    while( tau < max_tau )
    {
        if(system[1] < 2.)
        {
            cout << "Program terminated, r < 2.\n";
            return 1;
        }
 
        //Change this to record timesteps as needed
        if(counter%every==0) write_system( &data, system, tau);

        /*
        Using the third order runge-kutta scheme
        b1 = 2/9  b2 = 3/9  b3=4/9
                  c2 = 1/2  c3=3/4
        There is nothing special about this one, 
        I just derived it once and like it
        */
        calc_td( td, system);
        for(int i=0; i<6; i++)
            system1[i] = system[i] + 0.5*dtau*td[i];

        calc_td( td1, system1);
        for(int i=0; i<6; i++)
            system2[i] = system[i] + 0.75*dtau*td1[i];

        calc_td( td2, system2 );
        for(int i=0; i<6; i++)
            system[i] += dtau*(2.*td[i] + 3.*td1[i] + 4*td2[i])/9.;
 
        tau += dtau;
        counter++;
    }
    //Write data a final time for compatibility with gnuplot script
    write_system(&data, system, tau);

    data.close();

    return 0;
}
