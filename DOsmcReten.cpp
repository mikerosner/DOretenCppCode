////////////////////////////////////////////////////////////////////////////////////////
//
// DOsmcReten.cpp
// 
// Originally written by Michael Rosner, and modified by various others in Myung's lab
// Copyright 2009 by Jay I. Myung (Ohio State University: myung.1@osu.edu)
// Version: August 27, 2009
// Compiler: Microsoft Visual C++ 2005 Express Edition
//
// This program implements the design optimization algorithm described in the following
// article:
//
//   Myung, J. I. and Pitt, M. A. (2009). Optimal experimental design. Psychological 
//       Review, 116. 499-518.
//
// The problem to be solved is to find an optimal design defined as a set of N time 
// intervals {t_1, ..., t_N} that best descriminates two retention functions, power and 
// exponential, defined as
//
//   Power (POW): p = a(t+1)^(-b)    (0 < a < 1; 0.001 < b < 3.000)
//   Exponential (EXP): p = a*exp(-b*t)   (0 < a < 1; 0.001 < b < 3.000) 
//
// where p denotes the probability of correct recall as a function of time interval t, 
// and a and b are model parameters.
//
// The optimal design solution is sought over a series of Monte Carlo iterations by 
// simulating an artificial probability distribution. Specifically, a sequential Monte 
// Carlo (SMC, often called 'particle filtering') method with simulated annealing is 
// employed to numerically approximate the solution. In particular, this
// program implements the Resampling-Markov (R-M) algorithm described in page 776 of 
// the following article:
//
//   Amzal, B., Bois, F. Y., Parent, E. and Robert, C. P. (2006). Bayesian-optimal 
//      Design via Interacting Particle Systems. Journal of the American Statistical 
//      Association, 101, 773-785.
//
/////////////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <time.h>
#include <math.h>
#include "randoma.h" // needed to generate random numbers
using namespace std;

#include "FastMatrix.h" // collection of custom-designed fast matrix calculation routines

#define USERTYPE double
typedef MATRIX<USERTYPE> MATRIXTYPE;
#define SN 100  // binomial sample size
#define THIN 300 // thinning parameter of MCMC sampling
#define MCSMP 1000 // number of Monte Carlo samples used in the estimation of FIA complexity

#include "ClientFunctions.h" // collection of functions called by this program


inline void DOsmcReten(fstream& output) //PRIMARY ALGORITHM FUNCTION
{
	const int N = 5;		// number of time intervals 
	const int NP = 50;		// number of particles or MCMC chains
	const int NITER = 300;	// number of SMC iterations

	int iter=0, J, ip, j, i, Jp, Jd;
	USERTYPE u, uc, utilc, alpha, qratio=1, neff=0;

	MATRIXTYPE data(NITER,N+4); // stores the output of each SMC iteration
	MATRIXTYPE xx(NP,N), xxj(NP, N), xxc(1,N), xxf(1,N), tempxxj; // design solutions			
	MATRIXTYPE utiljp(1,NP), utilj(1,NP), util(1,NP); // local utility functions
	MATRIXTYPE theta, parma, parmb; // parameter values drawn from Jeffreys priors
	MATRIXTYPE ppow, pexp, temppow, tempexp; // model predictions

	MATRIXTYPE wgt(1,NP); // normalized weights calculated in the Reweighting Step of the R-M algorithm
	MATRIX<int> mno(1,NP); // resampled draws in the Resample Step of the R-M algorithm	
	MATRIXTYPE accept(1,NP); // acceptance probability estimated in the Metropolis-Hastings Step of the R-M algorithm

	USERTYPE Exp0[2], Pow0[2]; // initial parameter values in maximum likelihood estimation
	Exp0[0]=1;Exp0[1]=.5;Pow0[0]=1;Pow0[1]=.5;
 		
	MATRIXTYPE mc(NP,2), mcc(1,2); // Fisher Information Approximation (FIA) complexity of POW and EXP
	mc=0;mcc=0;
	
	MATRIXTYPE temp1(NP,1), temp2(1, N), temp3(NP, N);  
	MATRIXTYPE onesJ, zerosN(1,N), onesN(1,N), thousandsN(1,N);

	// ----- Intialize variables and matrices
	temp1=1;temp2=1;temp3.Random(&RandomDouble,1);
	zerosN=0;onesN=1;thousandsN=1000;

	xxj.InsertMatrix(temp1*5*temp2+(temp3*2-1)*0.05,0,0);
	iter=0;// intialize iteration number
	Jp=10; // initial number of iid draws (i.e., experiments) from the target distribution

	theta.reset(Jp,2);parma.reset(Jp,N);parmb.reset(Jp,N);ppow.reset(Jp,N);pexp.reset(Jp,N);
	onesJ.reset(Jp,1);onesJ=1;

	// ------------------------------------------------------
	//    1. Intialization Step of the R-M Algorithm
	// ------------------------------------------------------
		ip=0;
		while(ip<NP)
		{
			tempxxj=xxj.GetRow(ip);
			
			// --- Draw from Jeffreys priors
			for(int ii=0; ii<Jp; ii++)
			{
 			  theta.InsertMatrix(jeffpow(tempxxj,1,THIN),ii,0); 
			    // parameter values sampled from Jeffreys prior for POW
            }

			parma=theta.GetColumn(0)*onesN;parmb=theta.GetColumn(1)*onesN;
			ppow=parma.MultMembs((onesJ*tempxxj+1).Pow(-parmb));
			mc[ip][0]=complpow(tempxxj,MCSMP);// FIA complexity of POW

			for(int ii=0; ii<Jp; ii++)
			{
 			  theta.InsertMatrix(jeffexp(tempxxj,1,THIN),ii,0);
                 // parameter values sampled from Jeffreys prior for EXP
            }
			
			parma=theta.GetColumn(0)*onesN;parmb=theta.GetColumn(1)*onesN;
			pexp=parma.MultMembs((-parmb.MultMembs(onesJ*tempxxj)).Exp());
			mc[ip][1]=complexp(tempxxj,MCSMP);// FIX complexity of EXP

			// --- Evaluate local utility function u(d)
            mcc=mc.GetRow(ip);
			u=0;j=0;
			while (j<Jp)
			{
				temppow=ppow.GetRow(j);tempexp=pexp.GetRow(j);
				u=u+log(util_pow(tempxxj,temppow, Exp0, Pow0, mcc, output)+
				    util_exp(tempxxj,tempexp, Exp0, Pow0, mcc, output));
				j++;
			};
			utiljp[0][ip]=exp(u);
			ip++;
		}; //ip
	
		// 
		xxc=xxj.SortByRow().MeanByColumn(); // optimal design d estimated as mean of designs across NP particles 
		data.InsertMatrix(xxc,iter, 0);
		mcc[0][0]=complpow(xxc,MCSMP);mcc[0][1]=complexp(xxc,MCSMP); // FIX complexity given current optimal design
		data[iter][N]=utilestm(xxc,mcc,N,output); // "quick" estimate of global utility U(d) for current optimal design d
	  	data[iter][N+1]=((utiljp).Log()/(double)Jp).Mean(); // mean of log utility function across NP particles
		data[iter][N+2]=Jp;
		data[iter][N+3]=1; // acceptance probability initially set to 1
		iter=iter+1;
		
	// ------------------------------------------------------
	//    2. Reweighting Step of the R-M Algorithm
	// ------------------------------------------------------

	while(iter<NITER)
	{
		if ((iter+1)%1==0) 
  		cout << "iter " << iter << " : ";

        J=Round(2*log(Truncate(double(iter)/10)+1.0))+10; 
	       // anneaing schedule for increasing the value of J where the 'annealing temperature' = 1/J
	    Jd=J-Jp; 
	
        if (Jd > 0) // *** start of Big IF
        {  
	      // --- draw a new experiment if Jd > 0, and update both utility function and design solution	

		  theta.reset(Jd,2);parma.reset(Jd,N);parmb.reset(Jd,N);ppow.reset(Jd,N);	pexp.reset(Jd,N);
		  onesJ.reset(Jd,1);onesJ=1;

    	  ip=0;
		  while(ip<NP)
		  {   
			tempxxj=xxj.GetRow(ip);

			// --- Draw from Jeffreys priors
			for(int ii=0; ii<Jd; ii++)
			{
 			  theta.InsertMatrix(jeffpow(tempxxj,1,THIN),ii,0);
			    // parameter values sampled from Jeffreys prior for POW
			}
			
			parma=theta.GetColumn(0)*onesN;parmb=theta.GetColumn(1)*onesN;
			ppow=parma.MultMembs((onesJ*tempxxj+1).Pow(-parmb));
			mc[ip][0]=complpow(tempxxj,MCSMP);// FIA complexity of POW

			for(int ii=0; ii<Jd; ii++)
			{
 			  theta.InsertMatrix(jeffexp(tempxxj,1,THIN),ii,0);
			    // parameter values sampled from Jeffreys prior for EXP
			}
			
			parma=theta.GetColumn(0)*onesN;parmb=theta.GetColumn(1)*onesN;
			pexp=parma.MultMembs((-parmb.MultMembs(onesJ*tempxxj)).Exp());
			mc[ip][1]=complexp(tempxxj,MCSMP);// FIA complexity of EXP

			// --- Evaluate local utility function
			mcc=mc.GetRow(ip);
			u=0;j=0;
			while (j<Jd)
			{
				temppow=ppow.GetRow(j);tempexp=pexp.GetRow(j);
				u=u+log(util_pow(tempxxj,temppow, Exp0, Pow0, mcc, output)+
				    util_exp(tempxxj,tempexp, Exp0, Pow0, mcc, output));
				j++;
			};
			utilj[0][ip]=exp(u);
			ip++;
		  }; //ip
	   
		  // --- Compute normalized weights
		  wgt=utilj/(utilj.Sum());

 	// ------------------------------------------------------
	//    3. Resampling Step of the R-M Algorithm
	// ------------------------------------------------------
	
		  neff=1./(wgt.MultMembs(wgt).Sum());// effective sample size
		  if (neff < float(NP)*0.95)
		  {
			  mno=multinorm(NP,wgt); // resampling according to a multinomial distribution
		  }
		  else
		  {
			  mno=1;
		  }

		  // --- Update utility function and design solution
		  util=0;xx=0;ip=0;j=0;
		  while(ip<NP)
		  { 
			i=0;
			while (i<mno[0][ip])
			{
				util[0][j]=utiljp[0][ip]*utilj[0][ip];
				xx.InsertMatrix(xxj.GetRow(ip), j, 0);
				j++;
				i++;
			};
			ip++;
		  };

        } 
        else // *** else of Big IF
        { 
     	   // --- keep current utility function and design solution
		   util=utiljp;
		   xx=xxj;
        }; // *** end of big IF

 	// ------------------------------------------------------------
	//    4. Metropolis-Hastings (M-H) Step of the R-M Algorithm
	// ------------------------------------------------------------
	   		
		theta.reset(J,2);parma.reset(J,N);parmb.reset(J,N);ppow.reset(J,N);pexp.reset(J,N);
		onesJ.reset(J,1);onesJ=1;
		ip=0;accept=0;

		//----- Draw candidate designs via M-H sampling
		while (ip < NP)
		{ 
			xxc=xx.GetRow(ip).MultMembs(onesN+randn(N)*.2);// 0.2
			
			if (((xxc > zerosN).Sum() < N) || ((thousandsN > xxc).Sum() < N)) 
				   // check the boundary conditions (i.e., 0 < t_i < 1000 for all i = 1, ..., N)
				{
				  qratio=0.0;
				}
			else
				{
                  // Compute acceptance probability under an asymmetric proposal distribution
				  qratio=0.0;
				  for(int ii=0; ii<N; ii++)
				  {
					u=pow((xx[ip][ii]-xxc[0][ii]),2)*(pow(xx[ip][ii],2)-pow(xxc[0][ii],2))/(2.0*0.04*pow(xx[ip][ii],2)*pow(xxc[0][ii],2));
				    qratio=qratio+log(xx[ip][ii])-log(xxc[0][ii])-u; 	
				  }
				  qratio=exp(qratio);
				 };
			
			// --- Draw from Jeffreys priors
			for(int ii=0; ii<J; ii++)
			{
 			  theta.InsertMatrix(jeffpow(xxc,1,THIN),ii,0);
			     // parameter values sampled from Jeffreys prior for POW
            }

			parma=theta.GetColumn(0)*onesN;
			parmb=theta.GetColumn(1)*onesN;//% J x N
			ppow=parma.MultMembs((onesJ*xxc+1).Pow(-parmb));//% J x N
			mc[ip][0]=complpow(xxc,MCSMP);// FIA complexity of POW

			// --- Draw from Jeffreys priors for EXP
			for(int ii=0; ii<J; ii++)
			{
 			  theta.InsertMatrix(jeffexp(xxc,1,THIN),ii,0);
			     // parameter values sampled from Jeffreys prior for EXP
			}
			
			parma=theta.GetColumn(0)*onesN;
			parmb=theta.GetColumn(1)*onesN;//% J x N
			pexp=parma.MultMembs((-parmb.MultMembs(onesJ*xxc)).Exp());//% J x N
			mc[ip][1]=complexp(xxc,MCSMP);// FIA complexity of EXP

			// FIA complexity of POW
			mcc=mc.GetRow(ip);
			uc=0;j=0;
			while (j<J)
			{
				temppow=ppow.GetRow(j);tempexp=pexp.GetRow(j);
				uc=uc+log(util_pow(xxc,temppow, Exp0, Pow0, mcc, output)+
				      util_exp(xxc,tempexp, Exp0, Pow0, mcc, output));
				j++;
			};
			utilc=exp(uc);// utility of candidate design	
			
			alpha=(qratio*utilc)/util[0][ip];// acceptance probability

			if (alpha>1)
			{
				alpha=1;
			};	

			if (RandomDouble() < alpha)
			{
				xx.InsertMatrix(xxc,ip,0);
				util[0][ip]=utilc;
				accept[0][ip]=1;
				// accept the candidate design
			}   
			
			ip++;
		} // end of ip loop


		xxf=xx.SortByRow().MeanByColumn(); // optimal design (d) estimated as mean of designs across NP particles 
		data.InsertMatrix(xxf,iter, 0);
		mcc[0][0]=complpow(xxf,MCSMP);mcc[0][1]=complexp(xxf,MCSMP); // FIX complexity given current optimal design
		data[iter][N]=utilestm(xxf,mcc,N,output); // "quick" estimate of global utility U(d) for current optimal design d
		data[iter][N+1]=((util).Log()/(double)J).Mean(); // mean of log utility function across NP particles
		data[iter][N+2]=J; // record J value
		data[iter][N+3]=accept.Mean(); //record mean acceptance probability

		// --- Print the current design solution on monitor screen
		cout << "(" << data[iter][0] << ", ";
		for(int ii=1; ii < (N-1); ii++)
		{
			cout << data[iter][ii] << ", ";
		};
			cout <<data[iter][N-1]<< ")   [ "  <<data[iter][N] << ",  " << data[iter][N+1] << " ]" << endl;
		
		xxj=xx;Jp=J;utiljp=util; // update designs and utility function
		iter++;
		};// end of iter loop


	// ----- Save the results in an ASCII output file
	data.Export(output);
};


int main(int argc, char* argv[])
{
	fstream output;
	
	#ifdef DEBUG
	try{
	#endif

	TRandomInit(int(time(NULL)));

	output.open("smcoutput.txt", fstream::out); //initialize output stream
	DOsmcReten(output);output << '\n';output.close();

	#ifdef DEBUG
	}
	catch(char* e)
	{
		output << e;
		cerr << e;
	}
	#endif
    
	return 0;
}
// END of DOsmcReten.cpp
