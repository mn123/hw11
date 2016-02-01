#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);

void step(	const double dt, const double dx, const int Nx, const double xmin, 
				const double k, cmplx* const f0, cmplx* const f1);

void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40;
  	const double xmax = 40;
	const double Tend = 10*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = dx;
  	double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10;
  	const double omega = 0.2;
	const double k = omega*omega;
	const double alpha = pow(k,0.25);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];
	cmplx* psi1 = new cmplx[Nx];
	cmplx* h = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {

         t+=dt;
			step(dt,dx,Nx,xmin,k,psi0,psi1);
			h = psi0;
			psi0 = psi1;
			psi1 = h;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;

	delete[] psi0;
  	delete[] psi1;
	return 0;
}
//-----------------------------------

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
//----------------------------------

void step(const double dt, const double dx, const int Nx, const double xmin, const double k, cmplx* const f0, cmplx* const f1)
{
	cmplx a = cmplx(0.0,-dt/(4*dx*dx));

	double x;
	double* V = new double[Nx];
	cmplx* d = new cmplx[Nx];
	cmplx* r = new cmplx[Nx]; // rechte Seite der Gleichung

	//V berechnen
	for(int j=0; j<Nx; j++){
		x = xmin + j*dx;
		V[j] = 0.5*k*x*x;
		d[j] = cmplx(1,dt/(2*dx*dx)+dt*V[j]/2);
	}
	
	//rechte Seite berechnen
	r[0] = f0[0]*(cmplx(1,0)+cmplx(2.0,0.0)*a-cmplx(0.0,dt*V[0]/2))-a*f0[1];
	
	for(int j=1; j<Nx-1; j++){
		r[j] = a*(-f0[j-1]+cmplx(2,0)*f0[j]-f0[j+1])+f0[j]-cmplx(0,dt*V[j]/2)*f0[j];
	}
	
	r[Nx-1] = -a*f0[Nx-2]+f0[Nx-1]*(cmplx(1,0)+cmplx(2,0)*a-cmplx(0,dt*V[Nx-1]/2));
	

	//forward substitution
	for(int j=1; j<Nx; j++){
		d[j] -= a*a/d[j-1];
		r[j] -= r[j-1]*a/d[j-1];
	}

	//backward substitution
	f1[Nx-1] = r[Nx-1]/d[Nx-1];
	
	for(int j=Nx-2; j>=0; j--){
		f1[j] = (r[j]-a*f1[j+1])/d[j];
	}

	delete[] V;
	delete[] r;
	delete[] d;
	
}
