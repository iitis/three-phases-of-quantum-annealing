//////////////////////////////////////////////////////////////////////////////////////////////////
// Time-dependent dynamics of the transverse field Ising chain // Source code by: Artur Soriani //
//////////////////////////////////////////////////////////////////////////////////////////////////

// Hamiltonian (in LaTeX): H = - \Gamma/2 \sum_{j=1}^{N} \sigma_j^x - J/2 \sum_{j=1}^{N} \sigma_j^z \sigma_{j+1}^z

// load libraries
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <complex>

// define pi
#ifndef M_PI
#define M_PI	acos(-1)
#endif

// namespace for complex numbers
using namespace std;
using namespace std::complex_literals;

// some strings to output to terminal
string			program		= "TI_optimal_driving: varying tau while crossing the critical point";
string			protocol	= "LIN";

// define parameters to be used
int				N_p			;						// number of spins; to be entered when running the program
long double		h			= pow(10.,-3.);			// timestep
long double		tau			;						// process duration; will be varied throughout the run of the program
long double		J			= 1.;					// coupling constant
long double		a			= 1.;					// lattice constant
int				points		= 100+1;				// number of points to plot
const int		gs			= 1;					// gs = 1 for starting in the ground state; changing may break the code

// x is the variable in the horizontal axis
long double		xmin		;						// minimum value of x
long double		xmax		;						// maximum value of x
bool			debug		= false;				// boolean for debugging
long double		t_i			;						// initial time
long double		t_f			;						// final time

// Gamma is the external field
long double		Gamma_i		= .5*J;				// initial value of Gamma
long double		Delta		= 1.*J;				// total variation of Gamma
long double		Gamma_f		= Gamma_i + Delta;		// final value of Gamma

// lambda is a redefinition of the field; not sure why I did this, probably because I was working exclusively with lambda (and not Gamma) initially
long double		lambda_i	= (Gamma_i - J)/Delta;	// initial value of lambda
long double		lambda_f	= (Gamma_f - J)/Delta;	// final value of lambda

// define the time dependence of Gamma
long double Gamma(long double t){
	long double s = (t - t_i)/tau;
	return Gamma_i + Delta * (
		1.* s						// in this case, we have a linear protocol
	);
}

// define lambda
long double lambda(long double t){
	return ( Gamma(t) - J )/Delta;
}

// define the possible values of momentum
long double k(int n){
	return ( n + 1./2. ) * 2.*M_PI/(N_p*a);
}

// define the energy eigenvalues in each momentum subspace; essentially a Landau-Zener system for each positive k
long double e(long double k, long double t){
	return pow( pow(J + Delta*lambda(t) - J*cos(k*a),2.) + pow(J*sin(k*a),2.) , 1./2.);
}

// redefine the arctangent function
long double aTan(long double x){
	if(x >= 0){ return atan(x);  }
	else{ return atan(x) + M_PI; }
}

// define the theta function, necessary to calculate eigenstates
long double theta(long double k, long double t){
	return 1./2. * aTan( J*sin(k*a)/( J + Delta*lambda(t) - J*cos(k*a) ) );
}

// define the time evolution of single timestep
void evolution(std::complex<long double> u[][gs], std::complex<long double> v[][gs], long double t[][gs]){
	// define necessary variables for the algorithm
	std::complex<long double> U1, U2, U3, U4, V1, V2, V3, V4;
	// loop over every eigenstate
	int i = 0, n;
	while(i < gs){
		
		// loop over every momentum subspace; again, a Landau-Zener system for each positive k
		n = 0;
		while(n < N_p/2){
		
			// fourth-order Runge-Kutta algorithm
			U1   =  1il*( ( J + Delta*lambda(t[n][i]) - J*cos(k(n)*a) )*u[n][i] + J*sin(k(n)*a)*v[n][i] );
			V1   =  1il*( J*sin(k(n)*a)*u[n][i] - ( J + Delta*lambda(t[n][i]) - J*cos(k(n)*a) )*v[n][i] );
		
			U2   =  1il*( ( Delta*lambda(t[n][i]+h/2.) + J - J*cos(k(n)*a) )*( u[n][i] + h/2. * U1 ) + J*sin(k(n)*a)*( v[n][i] + h/2. * V1 ) );
			V2   =  1il*( J*sin(k(n)*a)*( u[n][i] + h/2. * U1 ) - ( Delta*lambda(t[n][i]+h/2.) + J - J*cos(k(n)*a) )*( v[n][i] + h/2. * V1 ) );
		
			U3   =  1il*( ( Delta*lambda(t[n][i]+h/2.) + J - J*cos(k(n)*a) )*( u[n][i] + h/2. * U2 ) + J*sin(k(n)*a)*( v[n][i] + h/2. * V2 ) );
			V3   =  1il*( J*sin(k(n)*a)*( u[n][i] + h/2. * U2 ) - ( Delta*lambda(t[n][i]+h/2.) + J - J*cos(k(n)*a) )*( v[n][i] + h/2. * V2 ) );
		
			U4   =  1il*( ( Delta*lambda(t[n][i]+h) + J - J*cos(k(n)*a) )*( u[n][i] + h * U3 ) + J*sin(k(n)*a)*( v[n][i] + h * V3 ) );
			V4   =  1il*( J*sin(k(n)*a)*( u[n][i] + h * U3 ) - ( Delta*lambda(t[n][i]+h) + J - J*cos(k(n)*a) )*( v[n][i] + h * V3 ) );
		
			u[n][i] += h/6. * ( U1 + 2.L*U2 + 2.L*U3 + U4 );
			v[n][i] += h/6. * ( V1 + 2.L*V2 + 2.L*V3 + V4 );
			t[n][i] += h;
	
			n++;
		}
		
		i++;
	}
	
}

// start the program
int main(int narg, char **arg){
	
	if(narg - 1 != 4){
		cout<<"Four arguments required."<<endl;
		return 1;
	}
	
	// four inputs for runing the program: 1 - number of spins (must be even); 2 - 'lin' or 'log' scale; 3 - minimum value of x; 4 - maximum value of x;
	N_p = atoi(arg[1]);
	string scale_test = arg[2];
	if(scale_test == "lin"){ // when in linear scale, the inputs must be exactly the initial and final values of x
		xmin = atof(arg[3]);
		xmax = atof(arg[4]);
	}
	else if(scale_test == "log"){ // when in logarithmic scale, the inputs must be the log(base 10) of the initial and final values of x
		xmin = pow(10,atof(arg[3]));
		xmax = pow(10,atof(arg[4]));
	}
	else{
		cout<<"Second argument must be 'lin' or 'log'."<<endl;
		return 1;
	}
	
	if( N_p % 2 == 1){
		cout<<"Number of spins is not even."<<endl;
		return 1;
	}
	
	long double prex = M_PI*M_PI/(N_p*N_p) * J*J/Delta; // define the prefactor of tau in the horizonal axis, i.e., x = prex tau
	
	long double base = pow(xmax/xmin,1./(points-1)); // calculate the horizontal distance between two points in the plot
	
	// a few outputs to the terminal
	cout<<program<<endl;
	cout<<"Protocol: "<<protocol<<";"<<endl;
	cout<<"Number of spins: "<<arg[1]<<";"<<endl;
	if(scale_test == "lin"){
	cout<<"Minimum re-scaled duration: "<<arg[3]<<";"<<endl;
	cout<<"Maximum re-scaled duration: "<<arg[4]<<";"<<endl;
	}
	if(scale_test == "log"){
	cout<<"Minimum re-scaled duration: 10^"<<arg[3]<<";"<<endl;
	cout<<"Maximum re-scaled duration: 10^"<<arg[4]<<";"<<endl;
	}
    
    long double time = clock(); // define variable used for calculating the program duration
	
	// define the arrays containing the state vectors and time variables
	auto u = new std::complex<long double>[N_p/2][gs];
    auto v = new std::complex<long double>[N_p/2][gs];
	auto t = new long double[N_p/2][gs];
	
	int n; // define iterator for later use
	
	long double Wex; // define the variable containing the excess work for each tau
	long double taumax = xmax/prex; // calculate the maximum value of tau
	tau = xmin/prex;				// calculate the initial value of tau, which is also the minimum value of tau
	
	// define and open the files for outputing the data to be plotted; the name of the files contain the number of spins and the initial and final values of x
	ofstream fileWorkEx( ("TI_("+ string(arg[1]) + "_" + string(arg[3]) + "_" + string(arg[4]) + ")workEx.dat" ).c_str() );
	ofstream fileDebug;
	if(debug){fileDebug.open( ("TI_(" + string(arg[1]) + "_" + string(arg[3]) + "_" + string(arg[4]) + ")debug.dat" ).c_str() );}
	
	// loop over every value of tau
	while(tau < taumax+taumax/(2000.*points)){
		
		// calculate the initial and final times
		t_i = -tau/2;
		t_f = t_i + tau;
		
		// loop over every momentum subspace
		n = 0;
		while(n < N_p/2){
		
			// calculate the initial states, which are eigenstates of the initial Hamiltonian
			t[n][0] = t_i;
			u[n][0] = cos( theta(k(n),t[n][0]) );
			v[n][0] = sin( theta(k(n),t[n][0]) );
			
			/* 
			t[n][1] = t_i;
			u[n][1] = - sin( theta(k(n),t[n][1]) );
			v[n][1] =   cos( theta(k(n),t[n][1]) );
			*/
			
			n++;
		}
		
		while( t[0][0] < t_f ){evolution(u, v, t);} // loop over every timestep, numerically solving the dynamics
		
		Wex = 0; // resets the excess work from previous valus of tau
		
		// loop over every momentum subspace
		n = 0;
		while(n < N_p/2){
			
			// calculate the excess work
			t[n][0] = t_f;
			//t[n][1] = t_f;
			Wex += 2. * e(k(n),t[n][0]) * norm( sin(theta(k(n),t[n][0])) * u[n][0] - cos(theta(k(n),t[n][0])) * v[n][0] );
			
			n++;
		}
		
		fileWorkEx<<prex*tau<<" "<<Wex/(N_p * J)<<endl; // outputs x and the excess work per spin per coupling constant to file
		
		if(debug){fileDebug<<prex*tau<<" "<<norm(u[0][0])+norm(v[0][0])<<endl;} // when debugging, outputs x and the norm of the evolved lowest-momentum state to file (norm should be 1)
		
		// terminal text to follow the progress of the program
		cout<<"\33[2k\r"<<flush;
		cout<<"Progress: "<<round( 100* log( tau / (xmin/prex) ) / log( taumax / (xmin/prex) ) )<<"%"<<flush;
		
		// update tau for the next iteration
		tau = base * tau;
		
	}
	
	// terminal text to follow the progress of the program
	cout<<"\33[2k\r"<<flush;
	//cout<<"Progress: 100%"<<endl;
	
	// close files
	fileWorkEx.close();
	if(debug){fileDebug.close();}
	
	// delete arrays to save memory
	delete[] u;
	delete[] v;
	delete[] t;
	
	// calculate the program duration
	time = (clock()-time)/CLOCKS_PER_SEC;
	
	// outputs program duration to terminal
	if(time<60) cout<<time<<" second";
	else if(time/60<60) cout<<time/60<<" minute";
	else cout<<time/3600<<" hour";
	cout<<"s to compute.";
	return 0;
}