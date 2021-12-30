#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <time.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <iomanip>
#include <limits>

using namespace std;

//this is the function that solves the forward  problem (BTE) using adjoint (backward) Monte Carlo (MC) method:
double BTE_bkwrd(double x[]);

int main ()
{ 
    ofstream result;   //the file that the final optimized parameters will be written to
    int n={12};        //size of the optimized parameters (12=2*3*2: 2 functions (2 phonon branches), each a piecewise linear functions with 3 lines, where each line is parameterized with 2 parameters
    int maxfun=10000;  //max number of function calculations allowed in the optimization
    int maxiter=10000; //max number of iterations allowed in the optimization
    double tolf=1e-30; //function tolerance
    double tolx=1e-30; //parameter tolerance
    double x[12]={6,11.3,3.5e+012,8,6.5e+013,11,6,11.3,3.4e+012,8.5,2.2e+013,10};       //initial condition
    
    double IC[n]; 
    for (int i=0;i<n;++i) {x[i]=x[i]*1.0; IC[i]=x[i];}                                  //initial simplex of NM
    
    //NM parameters:
    double rho={1}; 
    double chi={2}; 
    double psi={0.5}; 
    double sigma={0.95};

    //initialization of the simplex:
    double onesn[n];
    for (int i=0;i<n;++i) {onesn[i] = 1.0;}
    double two2np1[n];
    for (int i=0;i<n;++i) {two2np1[i] = i+2.0;}
    double one2n[n];
    for (int i=0;i<n;++i) {one2n[i] = i+1.0;}
    double xin[n];
    for (int i=0;i<n;++i) {xin[i]=x[i];}
    double v[n][n+1];
    for (int i=0;i<n;++i) {for (int j=0;j<n+1;++j) {v[i][j]=0.0;}}
    double fv[n+1];
    for (int i=0;i<n+1;++i) {fv[i] = 0.0;}
    for (int i=0;i<n;++i) {v[i][0]=xin[i];}
    fv[0]=BTE_bkwrd(x);
    int func_evals={1};
    int itercount={0};
    cout<<"Nelder-Mead algorithm with Lagarias correction. Created by Mojtaba Forghani \n\n"; 
    cout<<setw(5)<<"iteration"<<setw(15)<<"function eval";
    cout<<setw(15)<<"function val"<<"\n";
    cout<<setw(5)<<itercount<<setw(15)<<func_evals<<setw(15)<<fv[0]<<"\n";
    double usual_delta={0.02};
    double zero_term_delta ={0.00025};

    //construction of the initial simplex:
    for (int j=0;j<n;++j)
    {
        double y[n];
        for (int i=0;i<n;++i) {y[i]= xin[i];} 
        if (y[j] != 0) { y[j]=(1.0+usual_delta)*y[j];}
        else { y[j]=zero_term_delta;}
        for (int i=0;i<n;++i) { v[i][j+1]=y[i]; }
        for (int i=0;i<n;++i) {x[i]=y[i];}
        double f;
        f=BTE_bkwrd(x);
        fv[j+1]=f;
    }

    //sorting the simplex:
    int J[n+1]; for (int i=0;i<(n+1);++i) { J[i]=i; }
    for (int i=(n);i>=0;i--) { for (int j=1;j<=i;j++) { if (fv[j-1]>fv[j])
    {
        double a=fv[j-1]; int b=J[j-1];
        fv[j-1]=fv[j]; J[j-1]=J[j];
        fv[j]=a; J[j]=b; } } 
    }
    double v_temp[n][n+1];
    for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) {v_temp[j][i]=v[j][J[i]];}}
    for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) {v[j][i]=v_temp[j][i];}}
    string how = "initial simplex";
    itercount=itercount+1;
    func_evals=n+1;
    cout<<setw(5)<<itercount<<setw(15)<<func_evals;
    cout<<setw(15)<<fv[0]<<setw(25)<<how<<"\n";

    //iterations of the NM algorithm:
    while (func_evals<maxfun && itercount<maxiter) 
    {
        //sorting of simplex values at each iteration:
        double max_fv=abs(fv[0]-fv[1]);
        for (int i=1;i<(n+1);++i) 
        { 
            if (abs(fv[0]-fv[i])<max_fv) {max_fv=abs(fv[0]-fv[i]);}
        }
        double max_v=abs(v[0][0]-v[0][1]);
        for (int i=1;i<(n+1);++i) {for (int j=0;j<n;++j) 
        { 
            if (abs(v[j][0]-v[j][i])<max_v) {max_v=abs(v[j][0]-v[j][i]);}}
        }

        //stopping criteria:
        if (max_fv <= tolf && max_v<=tolx) {break;}
        double xbar[n];
        for (int i=0;i<n;++i) {xbar[i]=0.0;}
        for (int i=0;i<n;++i) {for (int j=0;j<n;++j) 
        {
            xbar[i]=xbar[i]+v[i][j]/n;}
        }

        //check if "expand", "reflect", "shrink" or "contract" step should take place:
        double xr[n];
        for (int i=0;i<n;++i) {xr[i]=(1.0+rho)*xbar[i]-rho*v[i][n];}
        for (int i=0;i<n;++i) {x[i]=xr[i];}
        double fxr;
        fxr=BTE_bkwrd(x);
        func_evals=func_evals+1;
        if (fxr < fv[0])
        {
            double xe[n];
            for (int i=0;i<n;++i) {xe[i]=(1.0+rho*chi)*xbar[i]-rho*chi*v[i][n];}
            for (int i=0;i<n;++i) {x[i]=xe[i];}
            double fxe;
            fxe=BTE_bkwrd(x);
            func_evals=func_evals+1;
            if (fxe<fxr) 
            { 
                for (int i=0;i<n;++i) {v[i][n]=xe[i];}
                fv[n]=fxe;
                how="expand";
            }
            else
            {
                for (int i=0;i<n;++i) {v[i][n]=xr[i];}
                fv[n]=fxr;
                how="reflect";
            } 
        }
        else
        { 
            if (fxr<fv[n-1]) 
            { 
                for (int i=0;i<n;++i) {v[i][n]=xr[i];}
                fv[n]=fxr;
                how="reflect";
            }
            else
            {
                if (fxr<fv[n])
                {
                    double xc[n];
                    for (int i=0;i<n;++i) 
                    {
                        xc[i]=(1.0+psi*rho)*xbar[i]-psi*rho*v[i][n];
                    }
                    for (int i=0;i<n;++i) {x[i]=xc[i];}
                    double fxc;
                    fxc=BTE_bkwrd(x);
                    func_evals=func_evals+1;
                    if (fxc<=fxr)
                    {
                        for (int i=0;i<n;++i) {v[i][n]=xc[i];}
                        fv[n]=fxc;
                        how="contract outside";
                    }
                    else {how="shrink";}
                }
                else
                {
                    double xcc[n];
                    for (int i=0;i<n;++i) 
                    {
                        xcc[i]=(1.0-psi)*xbar[i]+psi*v[i][n];
                    }
                    for (int i=0;i<n;++i) {x[i]=xcc[i];}
                    double fxcc;
                    fxcc=BTE_bkwrd(x);
                    func_evals=func_evals+1;
                    if (fxcc<fv[n])
                    {
                        for (int i=0;i<n;++i) {v[i][n]=xcc[i];}
                        fv[n]=fxcc;
                        how="contract inside";
                    }
                    else {how="shrink";}
                }
                if (how=="shrink")
                { 
                    for (int j=1;j<(n+1);++j)
                    {
                        for (int i=0;i<n;++i) 
                        {
                            v[i][j]=v[i][0]+sigma*(v[i][j]-v[i][0]);
                        }
                        for (int i=0;i<n;++i) {x[i]=v[i][j];}
                        fv[j]=BTE_bkwrd(x);
                    }
                    func_evals=func_evals+n;
                }
            }
        }

        //sorting of simplex values:
        int J[n+1]; for (int i=0;i<(n+1);++i) { J[i]=i; }
        for (int i=(n);i>=0;i--) { for (int j=1;j<=i;j++) { if (fv[j-1]>fv[j])
        {
            double a=fv[j-1]; int b=J[j-1];
            fv[j-1]=fv[j]; J[j-1]=J[j];
            fv[j]=a; J[j]=b; } } 
        }
        double v_temp[n][n+1];
        for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) 
        {
            v_temp[j][i]=v[j][J[i]];}
        }
        for (int i=0;i<(n+1);++i) {for (int j=0;j<n;++j) 
        {
            v[j][i]=v_temp[j][i];}
        }

        //updating and outputing each iteration result:
        itercount=itercount+1;
        cout<<setw(5)<<itercount<<setw(15)<<func_evals;
        cout<<setw(15)<<fv[0]<<setw(25)<<how<<"\n";
        result.open("result_5n_ord_3_2_1e4_loglog.txt");
        for (int i=0;i<n;++i) {result<<IC[i]<<" ";} result<<"\n";
        result<<fv[0]<<"\n"; for (int i=0;i<n;++i) { result<<v[i][0]<<" "; }
        result<<"\n";
        result.close();
    }
    for (int i=0;i<n;++i) {x[i]=v[i][0];}
    cout<<"final x is: \n";
    for (int i=0;i<n;++i) {cout<<"   "<<setw(15)<<x[i]<<"\n";}
    //result.close();
    cin>>n;
}

//**************************************************************************
//**************************************************************************

double BTE_bkwrd ( double x[12] ) 
{
    double y={0.0};     //objective function value

    //random generation setup:
    srand (time(NULL));
    typedef boost::mt19937 RNGType;
    RNGType rng( time(0) );
    boost::uniform_real<> uni_dist(0,1);
    boost::variate_generator< RNGType, boost::uniform_real<> > uni(rng, uni_dist);

    //physical (model) parameters:
    double PI = 3.14159265;
    double hbar=1.054517e-34;                            
    double boltz=1.38065e-23;                              
    double Teq = 300;   //base temperature
    int Nmode_Si_first=1;
    int Nmode_Si_LA_TA=1000;
    int Nmode_Si_last=1399;
    int Nmode_Al_first=1;
    int Nmode_Al_last=1498;
    int Ntt=101;        //number of measurements

    //load and setup Al and Si material data:
    ifstream file_Si("BVK_Si.txt");
    double BVK_Si[4][1399];
    if(file_Si.is_open())
    {
        for(int i = 0; i < 1399; ++i)
        {
        for(int j = 0; j < 4; ++j)
        {
            file_Si >> BVK_Si[j][i];
        }
        }
        file_Si.close();
    }  
    ifstream file_Al("BVK_Al.txt");
    double BVK_Al[3][1498];
    if(file_Al.is_open())
    {
        for(int i = 0; i < 1498; ++i)
        {
        for(int j = 0; j < 3; ++j)
        {
                file_Al >> BVK_Al[j][i];
        }
        }
        file_Al.close();
    }
    double SD_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        SD_Si[i] = BVK_Si[1][i+Nmode_Si_first-1];
    }
    double SD_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        SD_Al[i] = BVK_Al[1][i+Nmode_Al_first-1];
    }
    double Dom_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        Dom_Si[i] = BVK_Si[3][i+Nmode_Si_first-1];
    }
    double Dom_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        Dom_Al[i]=Dom_Si[0]; 
    }
    double V_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        V_Si[i] = BVK_Si[2][i+Nmode_Si_first-1];
    }
    double V_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        V_Al[i] = BVK_Al[2][i+Nmode_Al_first-1];
    }
    double F_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0; i < Nmode_Si_last-Nmode_Si_first+1; ++i)
    {
        F_Si[i] = BVK_Si[0][i+Nmode_Si_first-1];
    }
    double F_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i < Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        F_Al[i] = BVK_Al[0][i+Nmode_Al_first-1];
    }
    double BL=1/5.32e18/0.9/2/PI/2/PI;
    double BT=1/5.07e18/1.6/4/pow(PI,2);
    double wb=1.2e6;
    double A=3e-45;
    double tau_inv_Si[Nmode_Si_last-Nmode_Si_first+1];
    double tau_Si[Nmode_Si_last-Nmode_Si_first+1];
    
    //**************************************************************************
    //**************************************************************************
    //relaxation times function parameterization:
    
    //ORDER 1 TAU:
    double Del=0.25e13;
	double tau_inv_Si_temp[Nmode_Si_last-Nmode_Si_first+1];
    /*for (int i=0; i < Nmode_Si_LA_TA-Nmode_Si_first+1; ++i)
    {
        tau_inv_Si_temp[i]=pow(10.0,(x[1]-x[0])/(log10(F_Si[Nmode_Si_LA_TA-Nmode_Si_first])-
            log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
        tau_Si[i]=1.0/tau_inv_Si[i];
    }
    for (int i=(Nmode_Si_LA_TA-Nmode_Si_first+1);
	    i < (Nmode_Si_last-Nmode_Si_first+1); ++i)
    {
    	tau_inv_Si_temp[i]=pow(10.0,(x[1]-x[0])/(log10(F_Si[Nmode_Si_LA_TA-Nmode_Si_first])-
		    log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
        tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
        tau_Si[i]=1.0/tau_inv_Si[i];
    }*/
    
    //************************************************************************** 
    //************************************************************************** 
    
    //ORDER 2 TAU:    
    /*double a_1=(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]));
    double b_1=-log10(F_Si[0])*a_1+x[0];
    double a_2=(x[1]-x[3])/(log10(F_Si[999])-log10(x[2]));
    double b_2=-log10(x[2])*a_2+x[3];
    double X_1=log10(x[2]-Del);
    double Y_1=x[0]+x[3]-x[0]*log10((x[2]-Del)/F_Si[0])/log10(x[2]/F_Si[0]);
    double X_2=log10(x[2]+Del);
    double Y_2=x[3]+(x[1]-x[3])*log10(1.0+Del/x[2])/log10(F_Si[999]/x[2]);
    double a=-((x[1]-x[3])/log10(F_Si[999]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
	         pow(log10((x[2]+Del)/(x[2]-Del)),3.0)*log10(1.0-(Del/x[2])*(Del/x[2]));
    double b=0.5*(((x[1]-x[3])/log10(F_Si[999]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
             log10((x[2]+Del)/(x[2]-Del))-3*a*log10(x[2]*x[2]-Del*Del));
    double c=(x[3]-x[0])/log10(x[2]/F_Si[0])-3*a*log10(x[2]-Del)*log10(x[2]-Del)-2*b*log10(x[2]-Del);
    double d=x[0]-log10(F_Si[0])*(x[3]-x[0])/log10(x[2]/F_Si[0])+0.5*((x[1]-x[3])/
	         log10(F_Si[999]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/log10((x[2]+Del)/
			 (x[2]-Del))*log10(x[2]-Del)*log10(x[2]-Del)-0.5*a*log10(x[2]-Del)*log10(x[2]-
			 Del)*log10((x[2]+Del)*(x[2]+Del)*(x[2]+Del)/(x[2]-Del));
    for (int i=0;i<1000;++i)
    {
	    if (log10(F_Si[i])<X_1)
	    {
            tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else if (log10(F_Si[i])<X_2)
        {
            tau_inv_Si_temp[i]=pow(10,a*pow(log10(F_Si[i]),3.0)+b*pow(log10(F_Si[i]),2.0)+c*log10(F_Si[i])+d);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else
        {
            tau_inv_Si_temp[i]=pow(10,(x[1]-x[3])/(log10(F_Si[999])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
    }
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || log10(x[2])>log10(F_Si[999]))
    {
    	y=y+1e5;
    }
    //**************************************************************************
    //**************************************************************************
    /*a_1=(x[7]-x[4])/(log10(x[6])-log10(F_Si[0]));
    b_1=-log10(F_Si[0])*a_1+x[4];
    a_2=(x[5]-x[7])/(log10(F_Si[999])-log10(x[6]));
    b_2=-log10(x[6])*a_2+x[7];
    X_1=log10(x[6]-Del);
    Y_1=x[4]+x[7]-x[4]*log10((x[6]-Del)/F_Si[0])/log10(x[6]/F_Si[0]);
    X_2=log10(x[6]+Del);
    Y_2=x[7]+(x[5]-x[7])*log10(1.0+Del/x[6])/log10(F_Si[999]/x[6]);
    a=-((x[5]-x[7])/log10(F_Si[999]/x[6])-(x[7]-x[4])/log10(x[6]/F_Si[0]))/
	         pow(log10((x[6]+Del)/(x[6]-Del)),3.0)*log10(1.0-(Del/x[6])*(Del/x[6]));
    b=0.5*(((x[5]-x[7])/log10(F_Si[999]/x[6])-(x[7]-x[4])/log10(x[6]/F_Si[0]))/
             log10((x[6]+Del)/(x[6]-Del))-3*a*log10(x[6]*x[6]-Del*Del));
    c=(x[7]-x[4])/log10(x[6]/F_Si[0])-3*a*log10(x[6]-Del)*log10(x[6]-Del)-2*b*log10(x[6]-Del);
    d=x[4]-log10(F_Si[0])*(x[7]-x[4])/log10(x[6]/F_Si[0])+0.5*((x[5]-x[7])/
	         log10(F_Si[999]/x[6])-(x[7]-x[4])/log10(x[6]/F_Si[0]))/log10((x[6]+Del)/
			 (x[6]-Del))*log10(x[6]-Del)*log10(x[6]-Del)-0.5*a*log10(x[6]-Del)*log10(x[6]-
			 Del)*log10((x[6]+Del)*(x[6]+Del)*(x[6]+Del)/(x[6]-Del));*/
    /*for (int i=1000;i<(1399);++i)
    {
	    if (log10(F_Si[i])<X_1)
	    {
            tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else if (log10(F_Si[i])<X_2)
        {
            tau_inv_Si_temp[i]=pow(10,a*pow(log10(F_Si[i]),3.0)+b*pow(log10(F_Si[i]),2.0)+c*log10(F_Si[i])+d);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else
        {
            tau_inv_Si_temp[i]=pow(10,(x[1]-x[3])/(log10(F_Si[999])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
    }
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || log10(x[2])>log10(F_Si[999]))
    {
    	y=y+1e5;
	}
	/*for (int i=1000;i<(1399);++i)
    {
	    if (log10(F_Si[i])<X_1)
	    {
            tau_inv_Si_temp[i]=pow(10,(x[7]-x[4])/(log10(x[6])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[4]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else if (log10(F_Si[i])<X_2)
        {
            tau_inv_Si_temp[i]=pow(10,a*pow(log10(F_Si[i]),3.0)+b*pow(log10(F_Si[i]),2.0)+c*log10(F_Si[i])+d);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
        else
        {
            tau_inv_Si_temp[i]=pow(10,(x[5]-x[7])/(log10(F_Si[999])-log10(x[6]))*(log10(F_Si[i])-log10(x[6]))+x[7]);
            tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
        }
    }
    if (log10(F_Si[0])>x[6] || X_1>log10(x[6]) || log10(x[6])>X_2 || log10(x[6])>log10(F_Si[1398]))
    {
    	y=y+1e5;
    }*/
	
    //**************************************************************************
    //**************************************************************************
	
    //ORDRER 3 TAU:
    double a_1=(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]));
    double b_1=-log10(F_Si[0])*a_1+x[0];
    double a_2=(x[5]-x[3])/(log10(x[4])-log10(x[2]));
    double b_2=-log10(x[2])*a_2+x[3];
    double a_3=(x[1]-x[5])/(log10(F_Si[999])-log10(x[4]));
    double b_3=-log10(x[4])*a_3+x[5];
    double X_1=log10(x[2]-Del);
    double Y_1=x[0]+x[3]-x[0]*log10((x[2]-Del)/F_Si[0])/log10(x[2]/F_Si[0]);
    double X_2=log10(x[2]+Del);
    double Y_2=x[3]+(x[5]-x[3])*log10(1.0+Del/x[2])/log10(x[4]/x[2]);
    double X_3=log10(x[4]-Del);
    double Y_3=x[2]+x[5]-x[3]*log10((x[4]-Del)/x[2])/log10(x[4]/x[2]);
    double X_4=log10(x[4]+Del);
    double Y_4=x[5]+(x[1]-x[5])*log10(1.0+Del/x[4])/log10(F_Si[999]/x[4]);
    double a1=-((x[5]-x[3])/log10(x[4]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
             pow(log10((x[2]+Del)/(x[2]-Del)),3.0)*log10(1.0-(Del/x[2])*(Del/x[2]));
    double b1=0.5*(((x[5]-x[3])/log10(x[4]/x[2])-(x[3]-x[0])/log10(x[2]/F_Si[0]))/
             log10((x[2]+Del)/(x[2]-Del))-3*a1*log10(x[2]*x[2]-Del*Del));
    double c1=(x[3]-x[0])/log10(x[2]/F_Si[0])-3*a1*log10(x[2]-Del)*log10(x[2]-Del)-2*b1*log10(x[2]-Del);
    double d1=x[0]-log10(F_Si[0])*(x[3]-x[0])/log10(x[2]/F_Si[0])+0.5*((x[5]-x[3])/log10(x[4]/x[2])-
	         (x[3]-x[0])/log10(x[2]/F_Si[0]))/log10((x[2]+Del)/(x[2]-Del))*log10(x[2]-Del)*log10(x[2]-Del)
			  -0.5*a1*log10(x[2]-Del)*log10(x[2]-Del)*log10((x[2]+Del)*(x[2]+Del)*(x[2]+Del)/(x[2]-Del));
    double a2=-((x[1]-x[5])/log10(F_Si[999]/x[4])-(x[5]-x[3])/log10(x[4]/x[2]))/
              pow(log10((x[4]+Del)/(x[4]-Del)),3.0)*log10(1.0-(Del/x[4])*(Del/x[4]));
    double b2=0.5*(((x[1]-x[5])/log10(F_Si[999]/x[4])-(x[5]-x[3])/log10(x[4]/x[2]))/
              log10((x[4]+Del)/(x[4]-Del))-3*a2*log10(x[4]*x[4]-Del*Del));
    double c2=(x[5]-x[3])/log10(x[4]/x[2])-3*a2*log10(x[4]-Del)*log10(x[4]-Del)-2*b2*log10(x[4]-Del);
    double d2=x[3]-log10(x[2])*(x[5]-x[3])/log10(x[4]/x[2])+0.5*((x[1]-x[5])/log10(F_Si[999]/x[4])-
	          (x[5]-x[3])/log10(x[4]/x[2]))/log10((x[4]+Del)/(x[4]-Del))*log10(x[4]-Del)*log10(x[4]-Del)
	          -0.5*a2*log10(x[4]-Del)*log10(x[4]-Del)*log10((x[4]+Del)*(x[4]+Del)*(x[4]+Del)/(x[4]-Del));
    for (int i=0;i<1000;++i)
    {
    	if (log10(F_Si[i])<X_1)
    	{
    		tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
    		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_2)
        {
        	tau_inv_Si_temp[i]=pow(10,a1*pow(log10(F_Si[i]),3.0)+b1*pow(log10(F_Si[i]),2.0)+c1*log10(F_Si[i])+d1);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_3)
        {
        	tau_inv_Si_temp[i]=pow(10,(x[5]-x[3])/(log10(x[4])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_4)
        {
        	tau_inv_Si_temp[i]=pow(10,a2*pow(log10(F_Si[i]),3.0)+b2*pow(log10(F_Si[i]),2.0)+c2*log10(F_Si[i])+d2);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else
        {
        	tau_inv_Si_temp[i]=pow(10,(x[1]-x[5])/(log10(F_Si[999])-log10(x[4]))*(log10(F_Si[i])-log10(x[4]))+x[5]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
	}
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || X_2>X_3 || X_3>log10(x[4]) || log10(x[4])>X_4 || log10(x[4])>log10(F_Si[999]))
    {
    	y=y+1e5;
    }
    //**************************************************************************
    a_1=(x[9]-x[6])/(log10(x[8])-log10(F_Si[0]));
    b_1=-log10(F_Si[0])*a_1+x[6];
    a_2=(x[11]-x[9])/(log10(x[10])-log10(x[8]));
    b_2=-log10(x[8])*a_2+x[9];
    a_3=(x[7]-x[11])/(log10(F_Si[999])-log10(x[10]));
    b_3=-log10(x[10])*a_3+x[11];
    X_1=log10(x[8]-Del);
    Y_1=x[6]+x[9]-x[6]*log10((x[8]-Del)/F_Si[0])/log10(x[8]/F_Si[0]);
    X_2=log10(x[8]+Del);
    Y_2=x[9]+(x[11]-x[9])*log10(1.0+Del/x[8])/log10(x[10]/x[8]);
    X_3=log10(x[10]-Del);
    Y_3=x[8]+x[11]-x[9]*log10((x[10]-Del)/x[8])/log10(x[10]/x[8]);
    X_4=log10(x[10]+Del);
    Y_4=x[11]+(x[7]-x[11])*log10(1.0+Del/x[10])/log10(F_Si[999]/x[10]);
    a1=-((x[11]-x[9])/log10(x[10]/x[8])-(x[9]-x[6])/log10(x[8]/F_Si[0]))/
        pow(log10((x[8]+Del)/(x[8]-Del)),3.0)*log10(1.0-(Del/x[8])*(Del/x[8]));
    b1=0.5*(((x[11]-x[9])/log10(x[10]/x[8])-(x[9]-x[6])/log10(x[8]/F_Si[0]))/
        log10((x[8]+Del)/(x[8]-Del))-3*a1*log10(x[8]*x[8]-Del*Del));
    c1=(x[9]-x[6])/log10(x[8]/F_Si[0])-3*a1*log10(x[8]-Del)*log10(x[8]-Del)-2*b1*log10(x[8]-Del);
    d1=x[6]-log10(F_Si[0])*(x[9]-x[6])/log10(x[8]/F_Si[0])+0.5*((x[11]-x[9])/log10(x[10]/x[8])-
	    (x[9]-x[6])/log10(x[8]/F_Si[0]))/log10((x[8]+Del)/(x[8]-Del))*log10(x[8]-Del)*log10(x[8]-Del)
		 -0.5*a1*log10(x[8]-Del)*log10(x[8]-Del)*log10((x[8]+Del)*(x[8]+Del)*(x[8]+Del)/(x[8]-Del));
    a2=-((x[7]-x[11])/log10(F_Si[999]/x[10])-(x[11]-x[9])/log10(x[10]/x[8]))/
        pow(log10((x[10]+Del)/(x[10]-Del)),3.0)*log10(1.0-(Del/x[10])*(Del/x[10]));
    b2=0.5*(((x[7]-x[11])/log10(F_Si[999]/x[10])-(x[11]-x[9])/log10(x[10]/x[8]))/
        log10((x[10]+Del)/(x[10]-Del))-3*a2*log10(x[10]*x[10]-Del*Del));
    c2=(x[11]-x[9])/log10(x[10]/x[8])-3*a2*log10(x[10]-Del)*log10(x[10]-Del)-2*b2*log10(x[10]-Del);
    d2=x[9]-log10(x[8])*(x[11]-x[9])/log10(x[10]/x[8])+0.5*((x[7]-x[11])/log10(F_Si[999]/x[10])-
        (x[11]-x[9])/log10(x[10]/x[8]))/log10((x[10]+Del)/(x[10]-Del))*log10(x[10]-Del)*log10(x[10]-Del)
        -0.5*a2*log10(x[10]-Del)*log10(x[10]-Del)*log10((x[10]+Del)*(x[10]+Del)*(x[10]+Del)/(x[10]-Del));
	/*for (int i=1000;i<(1399);++i)
    {
    	if (log10(F_Si[i])<X_1)
    	{
    		tau_inv_Si_temp[i]=pow(10,(x[3]-x[0])/(log10(x[2])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[0]);
    		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_2)
        {
        	tau_inv_Si_temp[i]=pow(10,a1*pow(log10(F_Si[i]),3.0)+b1*pow(log10(F_Si[i]),2.0)+c1*log10(F_Si[i])+d1);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_3)
        {
        	tau_inv_Si_temp[i]=pow(10,(x[5]-x[3])/(log10(x[4])-log10(x[2]))*(log10(F_Si[i])-log10(x[2]))+x[3]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_4)
        {
        	tau_inv_Si_temp[i]=pow(10,a2*pow(log10(F_Si[i]),3.0)+b2*pow(log10(F_Si[i]),2.0)+c2*log10(F_Si[i])+d2);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else
        {
        	tau_inv_Si_temp[i]=pow(10,(x[1]-x[5])/(log10(F_Si[999])-log10(x[4]))*(log10(F_Si[i])-log10(x[4]))+x[5]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
	}
    if (log10(F_Si[0])>x[2] || X_1>log10(x[2]) || log10(x[2])>X_2 || X_2>X_3 || X_3>log10(x[4]) || log10(x[4])>X_4 || log10(x[4])>log10(F_Si[999]))
    {
    	y=y+1e5;
	}*/
	for (int i=1000;i<(1399);++i)
    {
    	if (log10(F_Si[i])<X_1)
    	{
    		tau_inv_Si_temp[i]=pow(10,(x[9]-x[0])/(log10(x[8])-log10(F_Si[0]))*(log10(F_Si[i])-log10(F_Si[0]))+x[6]);
    		tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_2)
        {
        	tau_inv_Si_temp[i]=pow(10,a1*pow(log10(F_Si[i]),3.0)+b1*pow(log10(F_Si[i]),2.0)+c1*log10(F_Si[i])+d1);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_3)
        {
        	tau_inv_Si_temp[i]=pow(10,(x[11]-x[9])/(log10(x[10])-log10(x[8]))*(log10(F_Si[i])-log10(x[8]))+x[9]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else if (log10(F_Si[i])<X_4)
        {
        	tau_inv_Si_temp[i]=pow(10,a2*pow(log10(F_Si[i]),3.0)+b2*pow(log10(F_Si[i]),2.0)+c2*log10(F_Si[i])+d2);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
        else
        {
        	tau_inv_Si_temp[i]=pow(10,(x[7]-x[11])/(log10(F_Si[999])-log10(x[10]))*(log10(F_Si[i])-log10(x[10]))+x[11]);
        	tau_inv_Si[i]=abs(tau_inv_Si_temp[i]);
            tau_Si[i]=1.0/tau_inv_Si[i];
		}
	}
    if (log10(F_Si[0])>x[8] || X_1>log10(x[8]) || log10(x[8])>X_2 || X_2>X_3 || X_3>log10(x[10]) || log10(x[10])>X_4 || log10(x[10])>log10(F_Si[1398]))
    {
    	y=y+1e5;
    }

    //**************************************************************************
    //**************************************************************************
	
    double tau_inv_Al[Nmode_Al_last-Nmode_Al_first+1];
    double tau_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0; i< Nmode_Al_last-Nmode_Al_first+1; ++i)
    {
        tau_inv_Al[i]=1e11;
        tau_Al[i]=1.0/tau_inv_Al[i];
    }
    double de_dT_Si[Nmode_Si_last-Nmode_Si_first+1];
    for (int i=0;i <Nmode_Si_last-Nmode_Si_first+1;++i)
    {
        de_dT_Si[i] = pow(hbar*F_Si[i]/Teq,2)/boltz*exp(hbar*F_Si[i]/boltz/Teq)/pow(exp(hbar*F_Si[i]/boltz/Teq)-1,2);
    }
    double de_dT_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0;i <Nmode_Al_last-Nmode_Al_first+1;++i)
    {
        de_dT_Al[i] = pow(hbar*F_Al[i]/Teq,2)/boltz*exp(hbar*F_Al[i]/boltz/Teq)/pow(exp(hbar*F_Al[i]/boltz/Teq)-1,2);
    } 
    double Lambda_Si[Nmode_Si_last-Nmode_Si_first+1];
    double Lambda_Si_max=0.0;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    { 
        Lambda_Si[i]=V_Si[i]/tau_inv_Si[i];
        if (Lambda_Si[i]>Lambda_Si_max)
        {
        Lambda_Si_max=Lambda_Si[i];
        }
    }
    double Lambda_Al[Nmode_Al_last-Nmode_Al_first+1];
    for (int i=0;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    { 
        Lambda_Al[i]=V_Al[i]/tau_inv_Al[i];
    }
    double C_Si[Nmode_Si_last-Nmode_Si_first+1];
    double C_Si_ave=0.0;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    { 
        C_Si[i]=SD_Si[i]*de_dT_Si[i]*Dom_Si[i];
        C_Si_ave=C_Si_ave+C_Si[i];
    }
    double C_Al[Nmode_Al_last-Nmode_Al_first+1];
    double C_Al_ave=0.0;
    for (int i=0;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    { 
        C_Al[i]=SD_Al[i]*de_dT_Al[i]*Dom_Al[i];
        C_Al_ave=C_Al_ave+C_Al[i];
    }
    double C_Al_Si_ratio;
    C_Al_Si_ratio=C_Al_ave/C_Si_ave;

    //------------------------------------------------------------------------------  
    //calculating the objective function:  

    int L_num={8};    //number of length-scales in the TTG experiments:
    double **tt;
    tt=new double*[L_num];
    for(int i=0;i<L_num;i++) {tt[i] = new double[Ntt];}

    //assigning time-steps to each length-scale:
    for(int j=0;j<Ntt;j++)
    { 
        tt[0][j]=1.0*j*1e-12;
        tt[1][j]=2.0*j*1e-12;  //tt[1][j]=0.5*j*1e-12;
        tt[2][j]=3.0*j*1e-12;  //tt[2][j]=0.5*j*1e-12;
        tt[3][j]=12.0*j*1e-12; //tt[3][j]=0.5*j*1e-12;
        tt[4][j]=40.0*j*1e-12; //tt[4][j]=0.5*j*1e-12;
        tt[5][j]=50.0*j*1e-12; //tt[5][j]=0.5*j*1e-12;
        tt[6][j]=50.0*j*1e-12; //tt[6][j]=0.5*j*1e-12;
        tt[7][j]=50.0*j*1e-12; //tt[7][j]=0.5*j*1e-12;
    }

    //length-scales:
    double L_val[L_num];
    L_val[0]=10e-9; L_val[1]=50e-9; L_val[2]=100e-9; 
    L_val[3]=500e-9; L_val[4]=1e-6; L_val[5]=5e-6; 
    L_val[6]=10e-6; L_val[7]=50e-6; //L_val[8]=100e-6;

    //number of MC particles for each length-scale:
    int N[L_num];
    N[0]=10000; N[1]=10000; N[2]=10000; N[3]=10000; N[4]=10000;
    N[5]=10000; N[6]=10000; N[7]=10000; N[8]=10000;

    //forward BTE simulations using adjoint MC method:
    double T_eff[L_num];
    for (int i=0;i<L_num;i++) { T_eff[i]= (double) 1.0/N[i];}
    double k_Si=0.0;
    for (int i=0;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    {
        k_Si = k_Si+SD_Si[i]*Dom_Si[i]*V_Si[i]*V_Si[i]/tau_inv_Si[i]*de_dT_Si[i]/3;
    }
    double k_Al=0.0;
    for (int i=0;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    {
        k_Al = k_Al+SD_Al[i]*Dom_Al[i]*V_Al[i]*V_Al[i]/tau_inv_Al[i]*de_dT_Al[i]/3;
    }
    double cumul_base_Si[Nmode_Si_last-Nmode_Si_first+1];
    cumul_base_Si[0] = SD_Si[0]*de_dT_Si[0]*Dom_Si[0];
    double cumul_base_Al[Nmode_Al_last-Nmode_Al_first+1];
    cumul_base_Al[0] = SD_Al[0]*de_dT_Al[0]*Dom_Al[0];
    double cumul_coll_Si[Nmode_Si_last-Nmode_Si_first+1];
    cumul_coll_Si[0] = SD_Si[0]*de_dT_Si[0]*Dom_Si[0]*tau_inv_Si[0];
    double cumul_coll_Al[Nmode_Al_last-Nmode_Al_first+1];
    cumul_coll_Al[0] = SD_Al[0]*de_dT_Al[0]*Dom_Al[0]*tau_inv_Al[0];
    for (int i=1;i<Nmode_Al_last-Nmode_Al_first+1;++i)
    {
        cumul_base_Al[i] = cumul_base_Al[i-1]+SD_Al[i]*de_dT_Al[i]*Dom_Al[i];
        cumul_coll_Al[i] = cumul_coll_Al[i-1]+SD_Al[i]*de_dT_Al[i]*tau_inv_Al[i]*Dom_Al[i];
    }
    for (int i=1;i<Nmode_Si_last-Nmode_Si_first+1;++i)
    {
        cumul_base_Si[i] = cumul_base_Si[i-1]+SD_Si[i]*de_dT_Si[i]*Dom_Si[i];
        cumul_coll_Si[i] = cumul_coll_Si[i-1]+SD_Si[i]*de_dT_Si[i]*tau_inv_Si[i]*Dom_Si[i];
    }
    double RAND1=uni();
    int pos_num=1;
    double pos[1]={0.0};
    for (int pos_ind=0;pos_ind<pos_num;++pos_ind)
    {
        double **T_grid;
       	T_grid=new double*[L_num];
	    for(int i=0;i<L_num;i++) {T_grid[i]=new double[Ntt];}
        for(int i=0;i<L_num;i++) {for (int j=0;j<Ntt;j++) {T_grid[i][j]=0.0;}}
        int N_max=N[0];
        for (int i=1;i<L_num;i++) {if (N[i]>N_max) {N_max=N[i];}}
        int loop_ctrl[L_num];
        for (int i=0;i<L_num;i++) {loop_ctrl[i]=1;}
        for (int counter=1;counter<=N_max;++counter)
        {
            double x0[L_num];
            x0[0]=pos[pos_ind];
            for (int i=1;i<L_num;i++) {x0[i]=x0[0];}
            double t0[L_num];
            t0[0]=0.0;
            for (int i=1;i<L_num;i++) {t0[i]=t0[0];}
            //------------------------------------------------------------------
            double Rx=uni();
            int i1 = 0;
            int i3;
            int i2;         
            i3 = Nmode_Si_last-Nmode_Si_first+1;
            i2 = (int) floor((i1+i3)/2.0);
            while (i3-i1>1)
            {
                if (Rx<cumul_base_Si[i2-1]/cumul_base_Si[Nmode_Si_last-Nmode_Si_first])
                {
                    i3  = i2;
                    i2 = (int) floor((i1+i3)/2.0);
                }
                else
                {
                    i1 = i2;
                    i2 = (int) floor((i1+i3)/2.0);
                }
            }   
            int ind_mod=i3-1;                                 
            //------------------------------------------------------------------
            double R_angle=2*uni()-1;
            double v[L_num];
            v[0]=V_Si[ind_mod]*R_angle;
            for (int i=1;i<L_num;i++) { v[i]=v[0]; }
            int finished = 0;
            int im[L_num];
            for (int i=0;i<L_num;i++) {im[i]=0;}
            double Delta_t[L_num];  
            Delta_t[0]=-tau_Si[ind_mod]*log(1-uni());
            for (int i=1;i<L_num;i++) {Delta_t[i]=Delta_t[0];} 
            while (finished==0)
            {
                double x1[L_num];
                for (int i=0;i<L_num;i++) {x1[i]=x0[i]+ Delta_t[i]*v[i];}
                double t1[L_num];
                for (int i=0;i<L_num;i++) {t1[i]=t0[i]+ Delta_t[i];}
                //**************************************************************
                int traj_exists[L_num];
                for (int i=0;i<L_num;i++) {traj_exists[i]=0;}
                int start_im[L_num];
                for (int i=0;i<L_num;i++) {start_im[i]=im[i];}
                for (int i=0;i<L_num;i++)
                {
                    if (loop_ctrl[i]==1)
                    {
                        while ((im[i]+1)<=Ntt && tt[i][im[i]]<t1[i])
                        {
                            im[i]=im[i]+1;
                            traj_exists[i]=1;
                        }
                        if (traj_exists[i]==1)
                        {
                            double Xpos[im[i]-start_im[i]];
                            for (int j=0;j<im[i]-start_im[i];++j)
                            {
                                Xpos[j]=x0[i]+v[i]*(tt[i][j+start_im[i]]-t0[i]);
                                T_grid[i][j+start_im[i]]=T_grid[i][j+start_im[i]]+ (double) T_eff[i]*cos(2*PI/L_val[i]*Xpos[j]);
                            }
                        }
                    }
                }
                //**************************************************************
                double Rx=uni();
                int i1 = 0;
                int i3 = Nmode_Si_last-Nmode_Si_first+1;
                int i2 = (int) floor((i1+i3)/2.0);
                while (i3-i1>1)
                {
                    if (Rx<cumul_coll_Si[i2-1]/cumul_coll_Si[Nmode_Si_last-Nmode_Si_first])
                    {
                        i3=i2;
                        i2=(int) floor((i1+i3)/2.0);
                    }
                    else
                    {
                        i1=i2;
                        i2=(int) floor((i1+i3)/2.0);
                    }
                }
                ind_mod=i3-1;
                //**************************************************************
                double R=2*uni()-1;
                for (int i=0;i<L_num;i++) {v[i]=V_Si[ind_mod]*R;}
                for (int i=0;i<L_num;i++) {x0[i]=x1[i];}
                for (int i=0;i<L_num;i++) {t0[i]=t1[i];} 
                for (int i=0;i<L_num;i++) {t0[i]=t1[i];}                              
                Delta_t[0]=-tau_Si[ind_mod]*log(1-uni());
                for (int i=1;i<L_num;i++) {Delta_t[i]=Delta_t[0];}
                for (int i=0;i<L_num;i++) {if (t0[i]>tt[i][Ntt-1])
                {
                    loop_ctrl[i]=0;}
                }
                int finish_det={0};
                for (int i=0;i<L_num;++i) {finish_det=finish_det+loop_ctrl[i];}
                if (finish_det==0) {finished = 1;}
            }
            for (int i=0;i<L_num;i++) 
            {
                if (counter>N[i]) 
                {
                    loop_ctrl[i]=0;
                }
                else
                {
                    loop_ctrl[i]=1;
                }
            }
        }
    //reading temperatures for each length-scale:
    ifstream file_res_10n("Si_cos_bkwrd_10n_poly_0n_100p.txt");
    ifstream file_res_50n("Si_cos_bkwrd_50n_poly_0n_200p.txt");
    ifstream file_res_100n("Si_cos_bkwrd_100n_poly_0n_300p.txt");
    ifstream file_res_500n("Si_cos_bkwrd_500n_poly_0n_1200p.txt");
    ifstream file_res_1m("Si_cos_bkwrd_1m_poly_0n_4n.txt");
    ifstream file_res_5m("Si_cos_bkwrd_5m_poly_0n_5n.txt");
    ifstream file_res_10m("Si_cos_bkwrd_10m_poly_0n_5n.txt");
    ifstream file_res_50m("Si_cos_bkwrd_50m_poly_0n_5n.txt");


    double res[L_num][102][2];
    if(file_res_10n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_10n >> res[0][j][i];}} file_res_10n.close();
    }  
    if(file_res_50n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_50n >> res[1][j][i];}} file_res_50n.close();
    } 
    if(file_res_100n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_100n >> res[2][j][i];}} file_res_100n.close();
    } 
    if(file_res_500n.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_500n >> res[3][j][i];}} file_res_500n.close();
    } 
    if(file_res_1m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_1m >> res[4][j][i];}} file_res_1m.close();
    } 
    if(file_res_5m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_5m >> res[5][j][i];}} file_res_5m.close();
    } 
    if(file_res_10m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_10m >> res[6][j][i];}} file_res_10m.close();
    } 
    if(file_res_50m.is_open()) {for(int i=0;i<2;++i) {for(int j=0;j<102;++j)
    {
        file_res_50m >> res[7][j][i];}} file_res_50m.close();
    }
    
    //calculating the objective function:
    for (int i=0;i<L_num;i++) {for (int j=0;j<Ntt;++j)
    {
        y=y+(double) abs(T_grid[i][j]-res[i][j+1][1])/800.0;
    }}
    y=y+0.01*abs(1.438416377372402E+02-k_Si)/1.438416377372402E+02;
    delete [] tt;
    for (int i=0;i<L_num;++i) {delete [] T_grid[i];} delete [] T_grid;
    for (int i=0;i<1399;i++) { if (tau_inv_Si[i]<=0.0) {y=y+1e5;}}
    }
    return y;
}
