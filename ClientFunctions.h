////////////////////////////////////////////////////////////////////////////////////////
//
// ClientFunctions.h
//
// Originally written by Michael Rosner, and modified by various others in Myung's lab
// Copyright 2009 by Jay I. Myung (Ohio State University: myung.1@osu.edu)
// Version: August 27, 2009
// Compiler: Microsoft Visual C++ 2005 Express Edition
//
// This program consists of a collection of functions that are called by DOsmcReten.cpp.
//
////////////////////////////////////////////////////////////////////////////////////////

const double PI=3.14159265358979323846;

inline double RandomDouble() //Generate random numbers between 0 and 1
{
	return TRandom();
}

inline double RandNorm()  //Generalte a standard normal random numer
{
	double U1, U2, x;
	U1=RandomDouble();
	U2=RandomDouble();
	x=sqrt(-2*log(U1)) * cos(2*PI*U2);
	return x;
};

inline MATRIXTYPE randn(int np)	//Generate a vector of standard normal random numbers
{
	MATRIXTYPE Result(1,np);
	Result.Random(&RandNorm);
	return Result;
};

inline MATRIXTYPE jeffpow(MATRIXTYPE &t, int n, int TH)
// Simulate Jeffreys prior for POW via an independent Metropolis algorithm
{
	MATRIXTYPE y;
	USERTYPE olda, oldb;
	int i, j;
	MATRIXTYPE X1, X2, X3, X4, X5;
	MATRIXTYPE Iaa;
	MATRIXTYPE Ibb;
	MATRIXTYPE Iab;
	USERTYPE x1, x2, x3;
	USERTYPE fisher;
	USERTYPE oldFI;
	USERTYPE newa,newb, newFI;
	y.reset(n,2);
	olda=RandomDouble();
	oldb=.001+2.999*RandomDouble();
	i=0;

	X1=(t+1).Pow(-oldb);
	X2=-(X1*olda-1);
	X3=(t+1).Log();		
	X4=X1.DivMembs(X2);
	X5=X3.MultMembs(X4);

	Iaa=X4;
	Ibb=X3.MultMembs(X5);
	Iab=-X5;

	
	x1= Iaa.Sum();
	x2= Ibb.Sum();
	x3= pow(Iab.Sum(), 2);
	fisher =  x1 * x2 - x3;

	oldFI = sqrt(fisher);

	while (i<n)
	{
		j=0;
		while (j<TH)
		{
			newa=RandomDouble();
			newb=.001+2.999*RandomDouble();

			X1=(t+1).Pow(-newb);
			X2=-(X1*newa-1);
			X3=(t+1).Log();		
			X4=X1.DivMembs(X2);
			X5=X3.MultMembs(X4);

			Iaa=X4;
			Ibb=X3.MultMembs(X5);
			Iab=-X5;

			fisher = Iaa.Sum() * Ibb.Sum() - pow(Iab.Sum(), 2);
			newFI = sqrt(fisher);

			if (oldFI==0)
			{
				olda=newa;
				oldb=newb;
				oldFI=newFI;
			}
			else
			{
				if (RandomDouble() < (newFI/oldFI))
				{
					olda=newa;
					oldb=newb;
					oldFI=newFI;
				};
			};
			j++;
		}
		y[i][0]=olda;
		y[i][1]=oldb;
		i++;
	}
	return y;
};


inline MATRIXTYPE jeffexp(MATRIXTYPE &t, int n, int TH)
	// Simulate Jeffreys prior for EXP via an independent Metropolis algorithm
{
	MATRIXTYPE y;
	USERTYPE olda, oldb;
	MATRIXTYPE X1, X2, X3, X4;
	int i, j;
	MATRIXTYPE Iaa;
	MATRIXTYPE Ibb;
	MATRIXTYPE Iab;
	USERTYPE fisher;
	USERTYPE oldFI;
	USERTYPE newa,newb, newFI;
	y.reset(n,2);
	i=0;

	olda=RandomDouble();
	oldb=.001+2.999*RandomDouble();

	X1=(t*(-oldb)).Exp();
	X2=-(X1*olda-1);
	X3=	X1.DivMembs(X2);
	X4= t.MultMembs(X3);
	Iaa=X3;
	Ibb=t.MultMembs(X4);
	Iab=-X4;

	fisher = Iaa.Sum() * Ibb.Sum() - pow(Iab.Sum(), 2);
	oldFI = sqrt(fisher);

	while (i<n)
	{
		j=0;
		while (j<TH)
		{
			newa=RandomDouble();
			newb=.001+2.999*RandomDouble();

			X1=(t*(-newb)).Exp();
			X2=-(X1*newa-1);
			X3=	X1.DivMembs(X2);
			X4= t.MultMembs(X3);
			Iaa=X3;
			Ibb=t.MultMembs(X4);
			Iab=-X4;	

			fisher = Iaa.Sum() * Ibb.Sum() - pow(Iab.Sum(), 2);
			newFI = sqrt(fisher);

			if (oldFI==0)
			{
				olda=newa;
				oldb=newb;
				oldFI=newFI;
			}
			else
			{
				if (RandomDouble() < (newFI/oldFI))
				{
					olda=newa;
					oldb=newb;
					oldFI=newFI;
				};
			};
			j++;
		}
		y[i][0]=olda;
		y[i][1]=oldb;
		i++;
	}
	return y;
};

inline USERTYPE complpow(MATRIXTYPE &t, int n)
// Compute FIA complexity of POW
{
	static USERTYPE y, yv;
	yv=0;
	USERTYPE a, b;
	int i;
	MATRIXTYPE X1, X2, X3, X4, X5;
	MATRIXTYPE Iaa;
	MATRIXTYPE Ibb;
	MATRIXTYPE Iab;
	USERTYPE fisher;
		
	i=0;
	while (i<n)
	{
		a=RandomDouble();
		b=.001+2.999*RandomDouble();

		X1=(t+1).Pow(-b);
		X2=-(X1*a-1);
		X3=(t+1).Log();		
		X4=X1.DivMembs(X2);
		X5=X3.MultMembs(X4);

		Iaa=X4;
		Ibb=X3.MultMembs(X5);
		Iab=-X5;

		fisher = Iaa.Sum() * Ibb.Sum() - pow(Iab.Sum(),2);
		yv=yv+sqrt(fisher)/n;
				
		i++;
	}
	y=log(float(SN)/float(2*PI))+log(yv*1.0*2.999);
	return y;
};

inline USERTYPE complexp(MATRIXTYPE &t, int n)
// Compute FIA complexity of POW
{
	static USERTYPE y, yv;
	yv=0;
	USERTYPE a, b;
	int i;
	MATRIXTYPE X1, X2, X3, X4;
	MATRIXTYPE Iaa;
	MATRIXTYPE Ibb;
	MATRIXTYPE Iab;
	USERTYPE fisher;
		
	i=0;
	while (i<n)
	{
		a=RandomDouble();
		b=.001+2.999*RandomDouble();
        
		X1=(t*(-b)).Exp();
		X2=-(X1*a-1);
		X3=	X1.DivMembs(X2);
		X4= t.MultMembs(X3);
		Iaa=X3;
		Ibb=t.MultMembs(X4);
		Iab=-X4;

		fisher = Iaa.Sum() * Ibb.Sum() - pow(Iab.Sum(), 2);
		
		yv=yv+sqrt(abs(fisher))/n;
				
		i++;
	}
	y=log(float(SN)/float(2*PI))+log(yv*1*2.999);
	return y;
};


inline void modelprdpow(USERTYPE *a, MATRIXTYPE &x, MATRIXTYPE &Result) // model prediction of POW
{
	x.AddFast(1, Result);
	Result.PowFast(-a[1], Result);
	Result.MultpFast(a[0], Result);
};

inline void modelprdexp(USERTYPE *a, MATRIXTYPE &x, MATRIXTYPE &Result) // model prediction of EXP
{	
	x.MultpFast(-a[1], Result);
	Result.ExpFast(Result);
	Result.MultpFast(a[0], Result);
};

inline USERTYPE mlikpow(USERTYPE *parm, MATRIXTYPE &x, MATRIX<int> &yobs) // minus log likelihood of POW
{
	static MATRIXTYPE pprd(x.GetHeight(), x.GetWidth());
	static USERTYPE yy;
	static int i;
	static int nc;

	modelprdpow(parm,x,pprd);
	yy=0;
	i=0;
	nc=yobs.GetWidth();
	while (i<nc)
	{
		yy=yy+yobs[0][i]*log(pprd[0][i]+1e-300)+(SN-yobs[0][i])*log(1-pprd[0][i]+1e-300);
		i++;
	}
	return -yy;
};
inline USERTYPE mlikexp(USERTYPE *parm, MATRIXTYPE &x, MATRIX<int> &yobs) // minus log likelihood of EXP
{
	static MATRIXTYPE pprd(x.GetHeight(), x.GetWidth());
	static USERTYPE yy;
	static int i;
	static int nc;

	modelprdexp(parm,x, pprd);
	yy=0;
	i=0;
	nc=yobs.GetWidth();
	while (i<nc)
	{
		yy=yy+yobs[0][i]*log(pprd[0][i]+1e-300)+(SN-yobs[0][i])*log(1-pprd[0][i]+1e-300);
		i++;
	}
	return -yy;
};


struct FMINTYPE	//Return value of fmincon
{
	MATRIXTYPE xmin;
	USERTYPE ymin;
};

inline FMINTYPE fmincon(USERTYPE (*func)(USERTYPE*, MATRIXTYPE&, MATRIX<int>&), USERTYPE *a0, USERTYPE *Amin, 
						USERTYPE *Amax, MATRIXTYPE &x, MATRIX<int> &yobs, int MaxIterations, fstream &output)
// Optimization algorithm for finding maximum likelihood estimates
{
	static int i;
	static bool xDir, yDir;
	static USERTYPE jump[2];
	static USERTYPE current, temp, an;
	FMINTYPE Result;
	Result.xmin.reset(1,2);
	Result.ymin=0;
	static USERTYPE a[2];
	i=0;
	xDir=true;
	yDir=true;
	
	a[0]=a0[0];
	a[1]=a0[1];
	
	jump[0]=(Amax[0]-Amin[0])/10;
	jump[1]=(Amax[1]-Amin[1])/10;
	current=func(a, x, yobs);
	while((jump[0] > 0.00001 && jump[1] > 0.00001) || (i<MaxIterations*10))
	{
			an=a[0];					//x axis shift
			if (xDir==true)
			{	//move positive
				a[0]+=jump[0];
				if(a[0]<Amax[0])
				{
					temp=func(a, x, yobs);

					if (temp<current)
					{
						current=temp;
						jump[0]*=2;
					}
					else
					{
						a[0]=an;
						jump[0]/=2;
						xDir=false;
					};
				}
				else
				{
					a[0]=an;
					jump[0]/=2;
					xDir=false;
				};
			}
			else //move negative
			{
				a[0]-=jump[0];
				if(a[0]>Amin[0])
				{
					temp=func(a, x, yobs);

					if (temp<current)
					{
						current=temp;
						jump[0]*=2;
					}
					else
					{
						a[0]=an;
						jump[0]/=2;
						xDir=true;
					};
				}
				else
				{
					a[0]=an;
					jump[0]/=2;
					xDir=true;
				};
			};						//end xaxis shift
	
			an=a[1];				//yaxis shift
			if (yDir==true)
			{
				a[1]+=jump[1];
				if(a[1]<Amax[1])
				{
					temp=func(a, x, yobs);

					if (temp<current)
					{
						current=temp;
						jump[1]*=2;
					}
					else
					{
						a[1]=an;
						jump[1]/=2;
						yDir=false;
					};
				}
				else
				{
					a[1]=an;
					jump[1]/=2;
					yDir=false;
				};
			}
			else
			{
				a[1]-=jump[1];
				if(a[1]>Amin[1])
				{
					temp=func(a, x, yobs);

					if (temp<current)
					{
						current=temp;
						jump[1]*=2;
					}
					else
					{
						a[1]=an;
						jump[1]/=2;
						yDir=true;
					};
				}
				else
				{
					a[1]=an;
					jump[1]/=2;
					yDir=true;
				};
			};						//end yaxis shift
			
		i++;
	};
	Result.xmin[0][0]=a[0];
	Result.xmin[0][1]=a[1];
	Result.ymin=current;
	return Result;
};

inline MATRIX<int> binornd(int NUM, MATRIXTYPE &P) // Generate binomial random numbers
{
	MATRIX<int> Result(P.GetHeight(), P.GetWidth());
	int j=0, l, n;
	n=P.Getn();
	while (j<n)
	{
		l=0;
		Result[0][j]=0;
		while (l<NUM)
		{
			if (RandomDouble()<P[0][j])
			{
				Result[0][j]++;
			}
			l++;
		};
		j++;
		
	};
	return Result;
};


inline MATRIX<int> multinorm(int NUM, MATRIXTYPE wgt) // Generate multinomial random numbers
{
	int i=0;
	MATRIX<int> Result(1, wgt.GetWidth());
	MATRIX<int> z;
	MATRIXTYPE onesNC(1,wgt.GetWidth());
	MATRIXTYPE onesN(NUM,1);
	MATRIXTYPE randN(NUM,1);
	MATRIXTYPE cprb;
	onesNC=1;
	onesN=1;

	cprb=wgt.CumSum();
	randN.Random(&RandomDouble,1,0);
	z=((randN*onesNC) > (onesN*cprb)).SumByRow();
	Result=0;
	
	while (i<NUM)
	{
		Result[0][z[0][i]]++;
		i++;
	};
	return Result;
};


inline USERTYPE util_pow(MATRIXTYPE &x, MATRIXTYPE &pobs, USERTYPE *Exp0, USERTYPE *Pow0, 
						 MATRIXTYPE &mcc, fstream &output) // Compute the local utility function for POW
{
	static MATRIX<int> yobs;
	static FMINTYPE EXP;
	static FMINTYPE POW;
	static USERTYPE y;

	static USERTYPE Amin[2];
	Amin[0]=0;
	Amin[1]=0.001;
	static USERTYPE Amax[2];
	Amax[0]=1;
	Amax[1]=3;

	yobs=binornd(SN,pobs);
	EXP=fmincon(&mlikexp,Exp0,Amin,Amax,x,yobs, 70, output);
	POW=fmincon(&mlikpow,Pow0,Amin,Amax,x,yobs, 70, output); 

	Exp0[0]=EXP.xmin[0][0];	Exp0[1]=EXP.xmin[0][1];
	Pow0[0]=POW.xmin[0][0];	Pow0[1]=POW.xmin[0][1];

	y=0.5*(EXP.ymin + mcc[0][1] > POW.ymin + mcc[0][0])+1e-6;
	return y;
};


inline USERTYPE util_exp(MATRIXTYPE &x, MATRIXTYPE &pobs, USERTYPE *Exp0, USERTYPE *Pow0, 
						 MATRIXTYPE &mcc, fstream &output)  // Compute the local utility function for EXP
{
	static MATRIX<int> yobs;
	static FMINTYPE EXP;
	static FMINTYPE POW;
	static USERTYPE y;

	static USERTYPE Amin[2];
	Amin[0]=0;
	Amin[1]=0.001;
	static USERTYPE Amax[2];
	Amax[0]=1;
	Amax[1]=3;

	yobs=binornd(SN,pobs);	
	EXP=fmincon(&mlikexp,Exp0,Amin,Amax,x,yobs, 70, output);
	POW=fmincon(&mlikpow,Pow0,Amin,Amax,x,yobs, 70, output); 

	Exp0[0]=EXP.xmin[0][0];	Exp0[1]=EXP.xmin[0][1];
	Pow0[0]=POW.xmin[0][0];	Pow0[1]=POW.xmin[0][1];

	y=.5*(POW.ymin + mcc[0][0] > EXP.ymin + mcc[0][1])+1e-6;
	return y;
};

inline USERTYPE utilestm(MATRIXTYPE &xx, MATRIXTYPE &mcc, int N, fstream &output) 
// Estimate the global utility value U(d) given design solution d
{
	const int NN = 250;
	static USERTYPE yy=0;

	static MATRIXTYPE datastore(2*NN, 1);
	datastore=0;
	
	static MATRIX<int> yobs;
	static MATRIXTYPE theta, thetapow, thetaexp;
	static MATRIXTYPE parma(1,N), parmb(1,N), ppow(1,N), pexp(1,N);
	static MATRIXTYPE onesN(1,N);
	static FMINTYPE EXP;
	static FMINTYPE POW;
	onesN=1;

	static USERTYPE a0[2];
	static USERTYPE Amin[2];
	Amin[0]=0;
	Amin[1]=0;
	static USERTYPE Amax[2];
	Amax[0]=1;
	Amax[1]=3;

	theta.reset(1,2);thetapow.reset(NN,2);thetaexp.reset(NN,2);
	thetapow=jeffpow(xx,NN,THIN);thetaexp=jeffexp(xx,NN,THIN);

	for(int ii=0; ii<NN; ii++)
	{		
		theta=thetapow.GetRow(ii);
		parma=theta.GetColumn(0)*onesN;
		parmb=theta.GetColumn(1)*onesN;
		ppow=parma.MultMembs((xx+1).Pow(-parmb));

		yobs=binornd(SN,ppow);
 
			a0[0]=RandomDouble();a0[1]=2*RandomDouble();
			POW=fmincon(&mlikpow,a0,Amin,Amax,xx,yobs,70,output);
			EXP=fmincon(&mlikexp,a0,Amin,Amax,xx,yobs,70,output); 
					
		datastore[ii][0]=(POW.ymin + mcc[0][0] < EXP.ymin + mcc[0][1]);
		
		theta=thetaexp.GetRow(ii);
		parma=theta.GetColumn(0)*onesN;
		parmb=theta.GetColumn(1)*onesN;
		pexp=parma.MultMembs((-parmb.MultMembs(xx)).Exp());

		yobs=binornd(SN,pexp);

			a0[0]=RandomDouble();a0[1]=2*RandomDouble();
			POW=fmincon(&mlikpow,a0,Amin,Amax,xx,yobs, 70, output);
			EXP=fmincon(&mlikexp,a0,Amin,Amax,xx,yobs, 70, output); 
				
		datastore[(NN+ii)][0]=(POW.ymin + mcc[0][0] > EXP.ymin +mcc[0][1]);	
	
	};

	MATRIXTYPE aa1=datastore.SubSet(0,NN,0,datastore.GetWidth()).MeanByColumn();
	MATRIXTYPE aa2=datastore.SubSet(NN,2*NN,0,datastore.GetWidth()).MeanByColumn();
    MATRIXTYPE aa3=(aa1+aa2)/2;
	
	yy=aa3[0][0];
	return yy;
};
// END of ClienFunctions.h
