//////////////////////////////////////////////////////////////////////////////////////
//
// FastMatrix.h
//
// This programs consists of various custom-designed functions for fast matrix 
//    calculations.
// 
// Written by Michael Rosner 
// Copyright 2009 by Jay I. Myung (Ohio State University: myung.1@osu.edu)
//    All rights reserved. Permission to use, copy and distribute content is granted
//    for non-commercial use only, provided that the above copyright notice appears
//    in all copies. 
// Version: August 27, 2009
// Compiler: Microsoft Visual C++ 2005 Express Edition
//
///////////////////////////////////////////////////////////////////////////////////////


template <class T> class MATRIX 
//A somewhat undynamic Matrix class whose Width and Height must be defined before it is used. 
// Contains many common functions
{
private:
	///<summary>
	///representation of matrix class
	///</summary>
	T* Cells;
	int Width;
	int Height;
	int n;
	typedef double (*Ftype)(double, double);

	void ImportHelper(fstream& input, int num, int w, bool firstline) //HELPER FUNCTION
	{
		T token;
		char c;
		c=input.peek();
		while ((c=='\t') || (c==' '))
		{
			c=input.get();
			c=input.peek();
		};
		
		if (firstline)
		{	
			if (c=='\n')
			{
				firstline=false;
			}
			else
				w++;
			input >> token;
			num++;
		}
		else
		{
			input >> token;
			num++;
		};
		if (!input.eof())
			ImportHelper(input, num, w, firstline); 
		else
		{
			Cells=new T[num];
			n=num;
			Width=w;
			Height=n/w;
		}
		Cells[num-1] = token;

	};

	void Squarify()	//INCREASES WIDTH OR HEIGHT TO MAKE MATRIX SQUARE IN SHAPE
	{
		if (Width==Height)
			return;
		if (Width>Height)
		{
			MATRIX X(Width, Width);
			X=0;
			X.InsertMatrix((*this),0,0);
			(*this) &= X;
		}
		else
		{
			MATRIX X(Height, Height);
			X=0;
			X.InsertMatrix((*this),0,0);
			(*this) &= X;
		}
	};

	void StrassenHelper(MATRIX &A, MATRIX &B, MATRIX &C) //BULK OF STRASSEN'S RECURSIVE ALGORITHM
	{
		int H2, W2, Hc, Wc;
		MATRIX Ai, Bi;
		bool Hb, Wb, Lb;

		if (A.Width<500 || B.Width<500 || A.Height<500)
		{
			C&=A*B;
			return;
		};
		Hc=A.Height;
		Wc=B.Width;

		A.Squarify();
		B.Squarify();

		Hb=(A.Height%2==1);
		Wb=(B.Width%2==1);
		Lb=(B.Height%2==1);

		if (Wb && Hb)
		{
			Ai.reset(A.Height+1, B.Width+1);
			Bi.reset(A.Height+1, B.Width+1);

			Bi=0;
			Ai=0;
			
			Bi.InsertMatrix(B, 0, 0);
			Ai.InsertMatrix(A, 0, 0);
			
		}
		else
			if (Hb)
			{
				Ai.reset(A.Height+1, B.Width);
				Bi.reset(A.Height+1, B.Width);

				Bi=0;
				Ai=0;
			
				Bi.InsertMatrix(B, 0, 0);
				Ai.InsertMatrix(A, 0, 0);
			}
			else
				if (Wb)
				{
					Ai.reset(A.Height, B.Width+1);
					Bi.reset(A.Height, B.Width+1);	

					Bi=0;
					Ai=0;
			
					Bi.InsertMatrix(B, 0, 0);
					Ai.InsertMatrix(A, 0, 0);
				}
				else
				{
					Ai=A;
					Bi=B;
				};
		H2=Ai.Height/2;
		W2=Ai.Width/2;
		MATRIX A11(Ai.SubSet(0, W2, 0, H2)), A12(Ai.SubSet(W2, Ai.Width, 0, H2)), A21(Ai.SubSet(0, W2, H2, Ai.Height)), A22(Ai.SubSet(W2, Ai.Width, H2, Ai.Height));
		MATRIX B11(Bi.SubSet(0, W2, 0, H2)), B12(Bi.SubSet(W2, Bi.Width, 0, H2)), B21(Bi.SubSet(0, W2, H2, Bi.Height)), B22(Bi.SubSet(W2, Bi.Width, H2, Bi.Height));
		MATRIX C11, C12, C21, C22;
		MATRIX M1, M2, M3, M4, M5, M6, M7;

		(*this).StrassenHelper((A11+A22), (B11+B22), M1);
		(*this).StrassenHelper((A21+A22), B11, M2);
		(*this).StrassenHelper(A11, (B12-B22), M3);
		(*this).StrassenHelper(A22, (B21-B11), M4);
		(*this).StrassenHelper((A11+A12), B22, M5);
		(*this).StrassenHelper((A21-A11), (B11+B12), M6);
		(*this).StrassenHelper((A12-A22), (B21+B22), M7);
		
		C11&= M1+M4-M5+M7;
		C12&= M3+M5;
		C21&= M2+M4;
		C22&= M1-M2+M3+M6;

		MATRIX Ci(C11, C12, C21, C22);
		C&=Ci.SubSet(0, Wc, 0, Hc);

	};

	void SortHelper(T* A, int L, int R) //RECURSIVE CALL FOR QUICKSORT
	{
		if (L>=R)
			return;
		int index=L;
		int l=L;
		T temp;
		while (l < R)
		{
			if (A[l]<A[R])
			{
				temp=A[l];
				A[l]=A[index];
				A[index]=temp;
				index++;
			};
			l++;
		};
		temp=A[R];
		A[R]=A[index];
		A[index]=temp;
		SortHelper(A, L, index-1);
		SortHelper(A, index+1, R);
	};

	void SortHelperWRanks(T* A, int* Ranks, int L, int R) //RECURSIVE CALL FOR QUICKSORT
	{
		if (L>=R)
			return;
		int index=L;
		int l=L;
		int tempRank;
		T temp;
		while (l < R)
		{
			if (A[l]<A[R])
			{
				temp=A[l];
				A[l]=A[index];
				A[index]=temp;

				tempRank=Ranks[l];
				Ranks[l]=Ranks[index];
				Ranks[index]=tempRank;
				index++;
			};
			l++;
		};
		temp=A[R];
		A[R]=A[index];
		A[index]=temp;

		tempRank=Ranks[R];
		Ranks[R]=Ranks[index];
		Ranks[index]=tempRank;

		SortHelperWRanks(A, Ranks, L, index-1);
		SortHelperWRanks(A, Ranks, index+1, R);
	};

public:
	MATRIX()				//CONSTRUCTOR
	{
		Width=0;
		Height=0;
		n=0;
	};

	MATRIX(MATRIX const &X)			//CONSTRUCTOR
	{
		Width=X.Width;
		Height=X.Height;
		n=X.n;
		int h=0;
		Cells = new T[n];

		while (h < n)
		{
			Cells[h] = X.Cells[h];
			h++;
		};
	};

	MATRIX(int H, int W)		//CONSTRUCTOR
	{
		Width = W;
		Height = H;
		n=H*W;
		
		if (W==0 || H==0)
			return;
		Cells = new T[n];
	};

	T& operator ()(int h, int w)
	{
		#ifdef DEBUG
				if (h>=Height)
					throw("h to large");
				if (w>=Width)
					throw("w to large");
		#endif
		return Cells[h*Width+w];
	};

	MATRIX(MATRIX &UpperLeft, MATRIX &UpperRight, MATRIX &LowerLeft, MATRIX &LowerRight) //CONSTRUCTOR-CREATES ONE MATRIX FROM 4 EQUAL-SIZED MATRICES; USED IN STRASSEN'S
	{
		Width=2*UpperLeft.Width;
		Height=2*UpperLeft.Height;
		n=Width*Height;
		Cells=new T[n];
		int h=0, Hx, Wx;

		while (h<n)
		{
			Hx=h/Width;
			Wx=h%Width;

			if ((Hx<UpperLeft.Height) && (Wx<UpperLeft.Width))
			{
				Cells[h]=UpperLeft.Cells[Hx*UpperLeft.Width+Wx];
			}
			else
				if (!(Hx<UpperLeft.Height) && Wx<(UpperLeft.Width))
				{
					Cells[h]=LowerLeft.Cells[(Hx-UpperLeft.Height)*LowerLeft.Width+Wx];
				}
				else
					if ((Hx<UpperLeft.Height) && !(Wx<UpperLeft.Width))
					{
						Cells[h]=UpperRight.Cells[Hx*UpperRight.Width+Wx-UpperLeft.Width];
					}
					else
						if (!(Hx<UpperLeft.Height) && !(Wx<UpperLeft.Width))
						{
							Cells[h]=LowerRight.Cells[(Hx-UpperLeft.Height)*LowerRight.Width+Wx-UpperLeft.Width];
						};
			h++;
		};
	};

	~MATRIX()					//DESTRUCTOR
	{
		if (n>0)
		{
			delete[] Cells;
		};
	};

	void clear() //DELETE ALL VALUES
	{
		if (n>0)
		{
			delete[] Cells;
			n=0;
			Width=0;
			Height=0;
		};
	};

	void reset(int H, int W) //REASSIGN THE WIDTH AND HEIGHT OF MATRIX-CLEARS CURRENT VALUES
	{
		(*this).clear();
		Width = W;
		Height = H;
		n=H*W;
		if (W==0 || H==0)
			return;
		Cells = new T[n];
	};

	T* operator [] (int H)  //ACCESSOR
	{

		#ifdef DEBUG
				if (H>=n)
					throw("Error at \"T& operator [] (int H)\": H is greater than array size");
		#endif

		return Cells+H*Width;
	};

	T& Current()	//RETURNS FIRST ELEMENT
	{
		return *Cells;
	};

	void operator = (MATRIX const &X) //ASSIGNMENT OPERATOR; MATRIX = MATRIX X
	{
		int h=0;
		if (Width==0 || Height==0)
		{
			Width=X.Width;
			Height=X.Height;
			n=X.n;
			int h=0;
			Cells = new T[n];
	
			while (h < n)
			{
				Cells[h] = X.Cells[h];
				h++;
			};
			return;
		};

		#ifdef DEBUG
				if (Width!=X.Width || Height!=X.Height)
					throw("Error at \"void operator = (MATRIX const &X)\": Matrices are not of equal size");
		#endif

		while (h < n)
		{
			Cells[h] = X.Cells[h];
			h++;
		};
	};

	void operator = (const T X) //ASSIGNMENT OPERATOR; MATRIX(:,:) = X
	{
		int h=0;

		#ifdef DEBUG
				if (Width==0 || Height==0)
					throw("Error at \"void operator = (T X)\": Cannot assign a constant to an undeclared matrix.");
		#endif

		while (h < n)
		{
			Cells[h] = X;
			h++;
		};
	};

	void ImportArray(const T *X) //ASSIGNMENT OPERATOR; MATRIX(:,:) = X
	{
		int h=0;

		#ifdef DEBUG
				if (Width==0 || Height==0)
					throw("Error at \"void operator = (const T *X)\": Cannot assign a constant to an undeclared matrix.");
		#endif

		while (h < n)
		{
			Cells[h] = X[h];
			h++;
		};
	};

	void operator &= (MATRIX &X) //SWAP THE VALUES OF TWO MATRICES
	{
		if (this->Cells==X.Cells)
		{
			return;
		}
		else
		{
			T* temp;

			Width=X.Width^Width;
			X.Width=X.Width^Width;
			Width=X.Width^Width;

			Height=X.Height^Height;
			X.Height=X.Height^Height;
			Height=X.Height^Height;

			temp=Cells;
			Cells=X.Cells;
			X.Cells=temp;

			n=Height*Width;
			X.n=X.Height*X.Width;
		};
	};	

	void Export(fstream& output) //TO FILE
	{
		int w, h=0;
		while (h < Height)
		{
			w=0;
			if (h>0)
				output << '\n';
			while(w < Width)
			{
				output << Cells[Width*h+w] << '\t';
				w++;
			};
			h++;
		};
	};

	void Import(fstream& input) //FROM FILE
	{
		if (n>0)
			return;
		ImportHelper(input, 0, 0, true);
	};

	MATRIX MultMembs(MATRIX const &X) //MULTIPICATION
	{
		int h=0;
		MATRIX Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX MultMembs(MATRIX const &X)\": Unequal Sizes");
				};
		#endif
		
		while(h < n)
		{
			Result.Cells[h] = Cells[h] * X.Cells[h];
			h++;
		};
		return Result;
	};

	MATRIX DivMembs(MATRIX const &X) //DIVISION
	{
		int h=0;
		MATRIX Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX DivMembs(MATRIX const &X)\": Unequal Sizes");
				};
		#endif
		
		while(h < n)
		{
			Result.Cells[h] = Cells[h] / X.Cells[h];
			h++;
		};
		return Result;
	};

	MATRIX operator ^ (int x) //EXPONENT-ONLY TAKES INTEGERS
	{
		int l=1;
		MATRIX Result;

		while (l<x)
		{
			Result &= (*this) * (*this);
			l++;
		};
		return Result;
	};

	T Det()
	{
		if (Height==1 && Width==1)
			return Cells[0];

		double Determinant=0;
		bool sign=true;
		MATRIX temp(Height-1,Width-1);
		for(int i=0;i<Width;i++)
		{
			for (int n=0;n<i;n++)
			{
				for (int p=1;p<Height;p++)
				{
					temp(p-1,n)=(*this)(p,n);
				};
			};

			for (int n=i+1;n<Width;n++)
			{
				for (int p=1;p<Height;p++)
				{
					temp(p-1,n-1)=(*this)(p,n);
				};
			};

			if (sign=true)
			{
				Determinant+=temp.Det();
				sign=false;
			}
			else
			{
				Determinant-=temp.Det();
				sign=true;
			};
		};
		return abs(Determinant);
	};

	MATRIX diag()
	{
		MATRIX temp;
		
		if (Width==1 || Height==1)
		{
			if (Width==1)
				temp.reset(Height,Height);
			else
				temp.reset(Width,Width);

			temp=0;
			for(int i=0;i<n;i++)
			{
				temp(i,i)=Cells[i];
			};
		}
		else
			if (Width<Height)
			{
				temp.reset(1,Width);
				for(int i=0;i<Width;i++)
				{
					temp(0,i)=(*this)(i,i);
				};
			}
			else
			{
				temp.reset(1,Height);
				for(int i=0;i<Height;i++)
				{
					temp(0,i)=(*this)(i,i);
				};
			};
		return temp;
	};

	MATRIX repmat(int y, int x)
	{
		MATRIX Result(y*Height, x*Width);
		for(int i=0;i<y;i++)
			for(int j=0;j<x;j++)
				for(int i1=0;i1<Height;i1++)
					for(int j1=0;j1<Width;j1++)
					{
						Result.Cells[i1*Result.Width+i*Result.Width*Height+j1+j*Width]=Cells[i1*Width+j1];
					};
		return Result;
	};

	MATRIX Pow(double x)	//A(h,w)^x for all h,w
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = pow(Cells[h], x);
			h++;
		};
		return Result;
	}

	MATRIX Abs()	//abs(A(h,w)) for all h,w
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = abs(Cells[h]);
			h++;
		};
		return Result;
	}

	MATRIX Pow(MATRIX const &X)		//A(h,w)^B(h,w) for all h,w
	{
		int h=0;
		MATRIX Result(Height,Width);
		if (X.Width!=Width || X.Height!=Height)
		{
			Result.clear();
			return Result;
		};
		while (h < n)
		{
			Result.Cells[h] = pow(Cells[h], X.Cells[h]);
			h++;
		};
		return Result;
	}

	MATRIX Exp() //COMPUTES e^A[x,y] for all x and y
	{
		int i=0;
		MATRIX Result(Height, Width);
		while (i<n)
		{
			Result.Cells[i]=exp(Cells[i]);
			i++;
		};
		return Result;
	};	

	MATRIX Log() //COMPUTES ln(A[x,y]) for all x and y
	{
		int i=0;
		MATRIX Result(Height, Width);
		while (i<n)
		{
			Result.Cells[i]=log(Cells[i]);
			i++;
		};
		return Result;
	}
	MATRIX operator * (MATRIX const &X) //MATRIX MULTIPICATION -WORKS
	{
		int w, h=0, l, W=X.Width, H=Height, offset2, offset3=0;
		T tempA;
		MATRIX Result(H,W);

		#ifdef DEBUG
				if (Width!=X.Height)
					throw("Error at \"MATRIX operator * (MATRIX const &X)\": Inappropriate sizes for linear multiplication");
		#endif

		while (h < H)
		{
			w=0;
			while(w < W)
			{
				l=0;
				tempA=0;
				offset2=w;
				while (l < Width)
				{
					tempA+=Cells[l] * X.Cells[offset2];
					offset2+=X.Width;
					l++;
				};
				Result.Cells[w] = tempA;
				w++;
			};
			Cells+=Width;
			Result.Cells+=X.Width;
			h++;
		};
		Result.Cells-=Result.n;
		Cells-=n;
		return Result;
	};
//////////////////////////////////////////////Faster simple arithmatic
	void ExpFast(MATRIX &Result) //COMPUTES e^A[x,y] for all x and y
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void ExpFast(MATRIX const &Result\": Result's size is not equal to This");
				};
		#endif

		static int i;
		i=0;
		while (i<n)
		{
			Result.Cells[i]=exp(Cells[i]);
			i++;
		};
	};	

	void LogFast(MATRIX &Result) //COMPUTES ln(A[x,y]) for all x and y
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void LogFast(MATRIX const &Result)\": Result's size is not equal to This");
				};
		#endif

		static int i;
		i=0;
		while (i<n)
		{
			Result.Cells[i]=log(Cells[i]);
			i++;
		};
	}
	void PowFast(T x, MATRIX &Result)	//A(h,w)^x for all h,w
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void PowFast(double x, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = pow(Cells[h], x);
			h++;
		};
	}

	void PowFast(MATRIX const &X, MATRIX &Result)		//A(h,w)^B(h,w) for all h,w
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void PowFast(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void PowFast(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = pow(Cells[h], X.Cells[h]);
			h++;
		};
	}

	void MultpFast(MATRIX const &X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void MultpFast(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void MultpFast(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};

				if (X.Cells == Result.Cells)
				{
					throw("Error at \"void MultpFast(MATRIX const &X, MATRIX &Result)\": Result may not be the same matrix as X");
				};

				if (this->Cells == Result.Cells)
				{
					throw("Error at \"void MultpFast(MATRIX const &X, MATRIX &Result)\": Result may not be the same matrix as This");
				};
		#endif

		static int w, h, l, offset2, offset3;
		static T temp;
		h=0;
		offset3=0;
		while (h < Result.Height)
		{
			w=0;
			while(w < Result.Width)
			{
				l=0;
				temp=0;
				offset2=w;
				while (l < this->Width)
				{
					temp+=this->Cells[l] * X.Cells[offset2];
					offset2+=X.Width;
					l++;
				};
				Result.Cells[w] = temp;
				w++;
			};
			this->Cells+=this->Width;
			Result.Cells+=Result.Width;
			h++;
		};
		Result.Cells-=Result.n;
		this->Cells-=this->n;
	};

	void MultMembsFast(MATRIX const &X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void MultMembsFast(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void MultMembsFast(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while(h < n)
		{
			Result.Cells[h] = Cells[h] * X.Cells[h];
			h++;
		};
	};

	void DivMembFast(MATRIX const &X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void DivMembFast(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void DivMembFast(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while(h < n)
		{
			Result.Cells[h] = Cells[h] / X.Cells[h];
			h++;
		};
	};

	void AddFast(MATRIX const &X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void AddFast(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void AddFast(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = Cells[h] + X.Cells[h];
			h++;
		};
	};

	void AddFast(T X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void AddFast(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif


		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = Cells[h] + X;
			h++;
		};
	};

	void SubtFast(MATRIX const &X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void SubtFast(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void SubtFast(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = Cells[h] - X.Cells[h];
			h++;
		};
	};

	void SubtFast(T X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void SubtFast(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = Cells[h] - X;
			h++;
		};
	};

	void MultpFast(T X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void MultpFast(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = Cells[h] * X;
			h++;
		};
	};

	void DivFast(T X, MATRIX &Result) //done
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void DivFast(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			Result.Cells[h] = Cells[h] / X;
			h++;
		};
	};

	void NegFast(MATRIX &Result) //MAKE THE MATRIX NEGATIVE
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void NegFast(MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		
		while (h < n)
		{
			Result.Cells[h] = -Cells[h];
			h++;
		};
	};

	void GreaterThan(MATRIX const &X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void GreaterThan(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void GreaterThan(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]>X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void LessThan(MATRIX const &X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void LessThan(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void LessThan(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]<X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void GreaterThanEqual(MATRIX const &X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void GreaterThanEqual(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void GreaterThanEqual(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]>=X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void LessThanEqual(MATRIX const &X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void LessThanEqual(MATRIX const &X, MATRIX &Result)\": Result's size is not equal to This");
				};

				if (X.Width!=Width || X.Height!=Height)
				{
					throw("Error at \"void LessThanEqual(MATRIX const &X, MATRIX &Result)\": X's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]<=X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void GreaterThan(T X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void GreaterThan(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]>X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void LessThan(T X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void LessThan(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]<X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void GreaterThanEqual(T X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void GreaterThanEqual(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]>=X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void LessThanEqual(T X, MATRIX &Result) //MATRIX GREATER-THAN
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void LessThanEqual(T X, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h < n)
		{
			if (Cells[h]<=X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
	};

	void TransPoseFast(MATRIX &Result) 
	{

		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void TransPoseFast(MATRIX &Result)\": Result's size is not equal to This");
				};

				if (this->Cells == Result.Cells)
				{
					throw("Error at \"void TransPoseFast(MATRIX &Result)\": Result may not be the same matrix as This");
				};
		#endif

		static int w, h, offset1, offset2;
		h=0;
		offset1=0;

		while (h < Height)
		{
			w=0;
			offset2=h;
			while(w < Width)
			{
				Result.Cells[offset2] = Cells[w];
				offset2+=Width;
				w++;
			};
			Cells+=Width;
			h++;
		};
		Cells-=n;
	};

	void SortByColumnFast(MATRIX &Result) //RETURNS A SORTED VERSION OF THIS
	{
		#ifdef DEBUG
				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void SortByColumnFast(MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		this->TransPoseFast(Result);
		Result.SortSelfByRow();
		Result = Result.TransPose();
	};

	void GetRowFast(int index, MATRIX &Result) //RETURNS A SINGLE-ROW MATRIX FROM CURRENT MATRIX
	{

		#ifdef DEBUG
				if (index>=Height)
				{
					Result.clear();
					throw("Error at \"void GetRowFast(int index, MATRIX &Result)\": index is greater than height of matrix");
				};

				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void GetRowFast(int index, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int w;
		w=0;
		while (w<Width)
		{
			Result.Cells[w] = Cells[Width*index+w];
			w++;
		};
	};

	void GetColumnFast(int index, MATRIX &Result) //RETURNS A SINGLE-COLUMN MATRIX FROM CURRENT MATRIX
	{

		#ifdef DEBUG
				if (index>=Width)
				{
					Result.clear();
					throw("Error at \"void GetColumnFast(int index, MATRIX &Result)\": index is greater than width of matrix");
				};

				if (Result.Width!=Width || Result.Height!=Height)
				{
					throw("Error at \"void GetColumnFast(int index, MATRIX &Result)\": Result's size is not equal to This");
				};
		#endif

		static int h;
		h=0;
		while (h<Height)
		{
			Result.Cells[h] = Cells[h*Width+index];
			h++;
		};
	};
//////////////////////////////////////////////End Faster Simple Arithmatic
	MATRIX Strassen (MATRIX &X) //STRASSEN'S MULTIPICATION
	{
		MATRIX Ci;
		StrassenHelper((*this), X, Ci);
		return Ci;
	};

	/*MATRIX operator / (MATRIX &X) //MATRIX DIVISION-NOT IMPLIMENTED
	{
		return X;
	};*/
	MATRIX operator + (MATRIX const &X) //MATRIX ADDITION
	{
		int h=0;
		MATRIX Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX operator + (MATRIX const &X)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			Result.Cells[h] = Cells[h] + X.Cells[h];
			h++;
		};
		return Result;
	};

	MATRIX operator - (MATRIX const &X) //MATRIX SUBTRACTION
	{
		int h=0;
		MATRIX Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX operator - (MATRIX const &X)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			Result.Cells[h] = Cells[h] - X.Cells[h];
			h++;
		};
		return Result;
	};

	MATRIX<int> operator > (MATRIX const &X) //MATRIX GREATER-THAN
	{
		int h=0;
		MATRIX<int> Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX<int> operator > (MATRIX const &X)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			if (Cells[h]>X.Cells[h])
				Result[0][h]=1;
			else
				Result[0][h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator < (MATRIX const &X) //MATRIX LESS-THAN
	{
		int h=0;
		MATRIX<int> Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX<int> operator < (MATRIX const &X)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			if (Cells[h]<X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator >= (MATRIX const &X) //MATRIX >=
	{
		int h=0;
		MATRIX<int> Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX<int> operator >= (MATRIX const &X)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			if (Cells[h]>=X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator <= (MATRIX const &X) //MATRIX <=
	{
		int h=0;
		MATRIX<int> Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX<int> operator <= (MATRIX const &X)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			if (Cells[h]<=X.Cells[h])
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator > (T X) //ect...
	{
		int h=0;
		MATRIX<int> Result(Height,Width);
		while (h < n)
		{
			if (Cells[h]>X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator < (T X) //ect...
	{
		int h=0;
		MATRIX<int> Result(Height,Width);
		while (h < n)
		{
			if (Cells[h]<X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator >= (T X) //ect...
	{
		int h=0;
		MATRIX<int> Result(Height,Width);
		while (h < n)
		{
			if (Cells[h]>=X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX<int> operator <= (T X) //ect...
	{
		int h=0;
		MATRIX<int> Result(Height,Width);
		while (h < n)
		{
			if (Cells[h]<=X)
				Result[h]=1;
			else
				Result[h]=0;
			h++;
		};
		return Result;
	};

	MATRIX operator / (T X) //MATRIX DIVISION BY SINGLE NUMBER
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = Cells[h] / X;
			h++;
		};
		return Result;
	};

	MATRIX operator * (T X) //MATRIX MULTIPICATION BY SINGLE NUMBER
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = Cells[h] * X;
			h++;
		};
		return Result;
	};

	MATRIX operator + (T X) //MATRIX ADDITION BY SINGLE NUMBER
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = Cells[h] + X;
			h++;
		};
		return Result;
	};

	MATRIX operator - (T X) //MATRIX SUBTRACTION BY SINGLE NUMBER
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = Cells[h] - X;
			h++;
		};
		return Result;
	};

	MATRIX operator - () //MAKE THE MATRIX NEGATIVE
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = -Cells[h];
			h++;
		};
		return Result;
	};

	MATRIX RunFunc(MATRIX const &X, Ftype Func) //RUNS USER-DEFINED FUNC(MATRIX, MATRIX X)
	{
		int h=0;
		MATRIX Result(Height,Width);

		#ifdef DEBUG
				if (X.Width!=Width || X.Height!=Height)
				{
					Result.clear();
					throw("Error at \"MATRIX RunFunc(MATRIX const &X, Ftype Func)\": X's size is not equal to This");
				};
		#endif

		while (h < n)
		{
			Result.Cells[h] = Func(Cells[h],X.Cells[h]);
			h++;
		};
		return Result;
	};

	MATRIX RunFunc(T X, Ftype Func) //RUNS USER-DEFINED FUNC(MATRIX, X)
	{
		int h=0;
		MATRIX Result(Height,Width);
		while (h < n)
		{
			Result.Cells[h] = Func(Cells[h],X);
			h++;
		};
		return Result;
	};

	MATRIX TransPose() //TRANSPOSES THE MATRIX -WORKS
	{
		int w, h=0, offset1=0, offset2;
		MATRIX<T> Result(Width,Height);
		while (h < Height)
		{
			w=0;
			offset2=h;
			while(w < Width)
			{
				Result.Cells[offset2] = Cells[w];
				offset2+=Result.Width;
				w++;
			};
			Cells+=Width;
			h++;
		};
		Cells-=n;
		return Result;
	};

	MATRIX SortByColumn() //RETURNS A SORTED VERSION OF THIS
	{
		MATRIX Result(Height, Width);
		Result=(((*this).TransPose()).SortByRow()).TransPose();
		return Result;
	};

	MATRIX SortByRow() //RETURNS A SORTED VERSION OF THIS
	{
		MATRIX Result(*this);
		int l=0, w=0;

		while (l<Height)
		{
			SortHelper(Result.Cells, w, w+Width-1);
			w+=Width;
			l++;
		};
		return Result;
	};

	void SortRowsByCol(int c)
	{
		MATRIX Result(Height, Width);

		double* temp = new double[Height];
		int* ranks = new int[Height];

		for (int i=0; i<Height; i++)
		{
			temp[i]=Cells[i*Width+c];
			ranks[i]=i;
		};

		SortHelperWRanks(temp, ranks, 0, Height-1);

		for (int i=0; i<Height; i++)
		{
			for(int j=0; j<Width; j++)
			{
				Result.Cells[i*Result.Width+j]=Cells[ranks[i]*Width+j];
			};
		};

		Result&=*this;
		delete[] temp;
		delete[] ranks;
	};
	/*void SortSelfByColumn() //SORTS THIS-NOT IMPLIMENTED
	{
		
	};
	*/
	void SortSelfByRow()
	{
		int l=0, w=0;

		while (l<Height)
		{
			SortHelper(Cells, w, w+Width-1);
			w+=Width;
			l++;
		};
	}

	MATRIX GetRow(int index) //RETURNS A SINGLE-ROW MATRIX FROM CURRENT MATRIX
	{

		#ifdef DEBUG
				if (index>=Height)
				{
					throw("Error at \"MATRIX GetRow(int index)\": index is greater than height of matrix");
				};
		#endif

		int w=0;
		MATRIX Result(1, Width);
		while (w<Width)
		{
			Result.Cells[w] = Cells[Width*index+w];
			w++;
		};
		return Result;
	};

	MATRIX GetColumn(int index) //RETURNS A SINGLE-COLUMN MATRIX FROM CURRENT MATRIX
	{

		#ifdef DEBUG
				if (index>=Width)
				{
					throw("Error at \"MATRIX GetColumn(int n)\": index is greater than width of matrix");
				};
		#endif

		static int h;
		h=0;
		MATRIX Result(Height, 1);
		while (h<Height)
		{
			Result.Cells[h] = Cells[h*Width+index];
			h++;
		};
		return Result;
	};

	void InsertMatrix(MATRIX const &X, int y, int x) //INSERTS MATRIX A AT y,x
	{
		static int h;
		h=0;
		
		#ifdef DEBUG
				if ((X.Width>Width) || (X.Height>Height))
					throw("Error at \"void InsertMatrix(MATRIX const &X, int y, int x)\": X will not fit in This");
				if (y<0 || x<0)
					throw("Error at \"void InsertMatrix(MATRIX const &X, int y, int x)\": x,y may not be negative");
				if (x+X.Width>Width)
					throw("Error at \"void InsertMatrix(MATRIX const &X, int y, int x)\": X is to far to the right");
				if (y+X.Height>Height)
					throw("Error at \"void InsertMatrix(MATRIX const &X, int y, int x)\": X is to far down");
		#endif

		while (h<X.n)
		{
			Cells[(h/X.Width+y)*Width+h%X.Width+x] = X.Cells[h];
			h++;
		};
	};

	MATRIX SubSet(int RowLowBnd, int RowUpBnd, int ColLowBnd, int ColUpBnd) //RETURNS A SUBSET OF MATRIX
	{
		int w,h=0, offset1=0, offset2, offset3=0, subwidth, subheight;
		subwidth=ColUpBnd-ColLowBnd;
		subheight=RowUpBnd-RowLowBnd;
		MATRIX Result(subheight, subwidth);

		#ifdef DEBUG
				if (ColUpBnd<ColLowBnd || RowUpBnd<RowLowBnd || ColLowBnd<0 || ColUpBnd>Width || RowLowBnd<0 || RowUpBnd>Height)
				{
					throw("Error at \"MATRIX SubSet(int ColLowBnd, int ColUpBnd, int RowLowBnd, int RowUpBnd)\": These bounds do not define a valid subset of This");
				};
		#endif
		
		offset2=RowLowBnd*Width+ColLowBnd;
		while (h < (subheight))
		{
			w=0;
			while(w < (subwidth))
			{
				Result.Cells[offset1+w] = Cells[offset2+w];
				w++;
			};
			offset1=offset1+subwidth;
			offset2=offset2+Width;
			h++;
		};
		return Result;
	};

	void Random(double (*RandFunc)(), double Multiply = 1, double Add = 0) //ASSIGNS RANDOM VALUES TO THE MATRIX
	{
		int j=0;
		while (j<n)
		{
			Cells[j] = RandFunc() * Multiply + Add;
			j++;
		};		

	};

	int GetWidth() //RETURNS WIDTH
	{
		return Width;
	};

	int GetHeight() //RETURNS HEIGHT
	{
		return Height;
	};

	int Getn() //RETURNS n
	{
		return n;
	};


	T Mean() //GETS THE MEAN OF ALL VALUES IN THE MATRIX
	{
		static int h;
		T sum;
		sum=0;
		h=0;
		
		while (h < n)
		{
			sum = sum + Cells[h];
			h++;
		};
		return sum/(T)n;
	};

	MATRIX MeanByRow() //GETS THE MEAN OF EACH ROW
	{
		int w, h=0, offset;
		T sum;
	
		MATRIX Means(1, Height);
	
		while (h < Height)
		{
			sum=0;
			w=0;
			offset=Width*h;
			while(w < Width)
			{
				sum = sum + Cells[offset+w];
				w++;
			};
			Means.Cells[h] = sum/(T)Width;
			h++;
		};
		return Means;
	};

	MATRIX MeanByColumn() //GETS THE MEAN OF EACH COLUMN
	{
		int w=0, h;
		T sum;
	
		MATRIX Means(1, Width);
	
		while (w < Width)
		{
			sum=0;
			h=0;
			while(h < Height)
			{
				sum = sum + Cells[Width*h+w];
				h++;
			};
			Means.Cells[w] = sum/(double)Height;
			w++;
		};
		return Means;
	};	

	T VarAll() //GETS THE VARIANCE OF ALL VALUES IN THE MATRIX
	{
		int h=0;
		T sum=0, mean;
		while (h < n)
		{
			sum = sum + Cells[h];
			h++;
		};
		mean = sum/n;
		sum = 0;
		h = 0;
		while (h < n)
		{
			sum = sum + (Cells[h]-mean) * (Cells[h]-mean);
			h++;
		};
		return sum/n;
	};

	T VarAll(T mean) //GETS THE VARIANCE OF ALL VALUES IN THE MATRIX, PROVIDED A MEAN
	{
		int h=0;
		T sum=0;
		while (h < n)
		{
			sum = sum + (Cells[h]-mean) * (Cells[h]-mean);
			h++;
		};
		return sum/n;
	};

	T Sum()
	{
		int l=0;
		T sum=0;
		while(l<n)
		{
			sum=sum+Cells[l];
			l++;
		};
		return sum;
	};

	MATRIX SumByRow()
	{
		int l, h=0;
		T sum;
		MATRIX Result(Height,1);
		while (h<Height)
		{
			l=0;
			sum=0;
			while(l<Width)
			{
				sum=sum+Cells[l];
				l++;
			};
			Cells+=Width;
			Result.Cells[h]=sum;
			h++;
		}
		Cells-=n;
		return Result;
	};

	MATRIX SumByColumn()
	{
		int l, w=0, offset1=0;
		T sum;
		MATRIX Result(1,Width);
		while (w<Width)
		{
			l=0;
			sum=0;
			offset1=w;
			while(l<Height)
			{
				sum=sum+Cells[offset1];
				offset1+=Width;
				l++;
			};
			Result.Cells[w]=sum;
			w++;
		}
		return Result;
	};

	MATRIX CumSum()
	{
		int l=1;
		MATRIX Result(Height,Width);
		Result.Cells[0]=Cells[0];
		while(l<n)
		{
			Result.Cells[l]=Result.Cells[l-1]+Cells[l];
			l++;
		};
		return Result;
	};

	MATRIX CumSumByRow()
	{
		int l, h=0;
		MATRIX Result(Height,Width);
		while (h<Height)
		{
			l=1;
			*(Result.Cells)=*Cells;
			while(l<Width)
			{
				Result.Cells[l]=Result.Cells[l-1]+Cells[l];
				l++;
			};
			Cells+=Width;
			Result.Cells+=Width;
			h++;
		}
		Cells-=n;
		Result.Cells-=n;
		return Result;
	};

	MATRIX CumSumByColumn()
	{
		int l, w=0, offset1=0;
		MATRIX Result(Height,Width);
		while (w<Width)
		{
			l=1;
			Result.Cells[w]=Cells[w];
			offset1=w+Width;
			while(l<Height)
			{
				Result[offset1]=Result[offset1-Width]+Cells[offset1];
				offset1+=Width;
				l++;
			};
			w++;
		}
		return Result;
	};

	T Min()
	{
		T min;
		min=Cells[0];
		for(int i=1; i<n; i++)
		{
			if (min>Cells[i])
			{
				min=Cells[i];
			};
		};
		return min;
	};

	T Max()
	{
		T max;
		max=Cells[0];
		for(int i=1; i<n; i++)
		{
			if (max<Cells[i])
			{
				max=Cells[i];
			};
		};
		return max;
	};

	MATRIX Reshape(int H, int W) //redefine Height and Width, n remaining the same
	{
		MATRIX Result(H,W);
		if (H*W!=n)
			return Result;
		Result.ImportArray(Cells);
		return Result;
	};
};//END OF MATRIX CLASS

template <class T> class CUBE
{
private:
	T *Cells;
	int Height;
	int Width;
	int Depth;
	int n;
	
public:
	CUBE()
	{
		n=0;
		Height=0;
		Width=0;
		Depth=0;
	};
	
	CUBE(int H,int W,int D)
	{
		n=H*W*D;
		Height=H;
		Width=W;
		Depth=D;

		if (n>0)
			Cells=new T[n];
	};

	~CUBE()
	{
		if (n>0)
		{
			delete[] Cells;
		};
	};

	void clear() //DELETE ALL VALUES
	{
		if (n>0)
		{
			delete[] Cells;
			n=0;
			Width=0;
			Height=0;
		};
	};

	void reset(int H, int W, int D) //REASSIGN THE WIDTH AND HEIGHT OF MATRIX-CLEARS CURRENT VALUES
	{
		(*this).clear();
		Width = W;
		Height = H;
		Depth = D;
		n=H*W;
		if (W==0 || H==0)
			return;
		Cells = new T[n];
	};
	
	MATRIX<T> reshape(int H, int W)
	{
		MATRIX<T> NewMat(H, W);
		NewMat.ImportArray(Cells);
		return NewMat;
	};
	
	void operator = (CUBE const &X) //ASSIGNMENT OPERATOR; MATRIX = MATRIX X
	{
		int h=0;
		if (Width==0 || Height==0 || Depth==0)
		{
			Width=X.Width;
			Height=X.Height;
			Depth=X.Depth;
			n=X.n;
			int h=0;
			Cells = new T[n];
	
			while (h < n)
			{
				Cells[h] = X.Cells[h];
				h++;
			};
			return;
		};

		#ifdef DEBUG
				if (Width!=X.Width || Height!=X.Height)
					throw("Error at \"void operator = (MATRIX const &X)\": Matrices are not of equal size");
		#endif

		while (h < n)
		{
			Cells[h] = X.Cells[h];
			h++;
		};
	};
	
	void operator = (const T X) //ASSIGNMENT OPERATOR; MATRIX(:,:) = X
	{
		int h=0;

		#ifdef DEBUG
				if (Width==0 || Height==0)
					throw("Error at \"void operator = (T X)\": Cannot assign a constant to an undeclared matrix.");
		#endif

		while (h < n)
		{
			Cells[h] = X;
			h++;
		};
	};
	
	void operator &= (CUBE &X) //SWAP THE VALUES OF TWO MATRICES
	{
		if (this->Cells==X.Cells)
		{
			return;
		}
		else
		{
			T* temp;

			Width=X.Width^Width;
			X.Width=X.Width^Width;
			Width=X.Width^Width;

			Height=X.Height^Height;
			X.Height=X.Height^Height;
			Height=X.Height^Height;

			Depth=X.Depth^Depth;
			X.Depth=X.Depth^Depth;
			Depth=X.Depth^Depth;

			temp=Cells;
			Cells=X.Cells;
			X.Cells=temp;

			n=Height*Width;
			X.n=X.Height*X.Width;
		};
	};	

	void GetMatrix(int x, int dimension, MATRIX<T> &Result)//0=Height; 1=Width; 2=Depth
	{
		switch (dimension)
		{
		case 0:
			for (int w=0;w<Width;w++)
			{
				for (int d=0;d<Depth;d++)
				{
					Result[d][w]=Cells[d*Width*Height+x*Width+w];
				};
			};
			break;
		case 1:
			for (int h=0;h<Height;h++)
			{
				for (int d=0;d<Depth;d++)
				{
					Result[h][d]=Cells[d*Width*Height+h*Width+x];
				};
			};
			break;
		case 2:
			for (int h=0;h<Height;h++)
			{
				for (int w=0;w<Width;w++)
				{
					Result[h][w]=Cells[x*Width*Height+h*Width+w];
				};
			};
			break;
		default:
			break;
		};
	};

	T& operator () (int H, int W, int D)
	{
		return Cells[D*Width*Height+H*Width+W];
	};

	void InsertMatrix(MATRIX<T> &Input, int x, int dimension)
	{
		switch (dimension)
		{
		case 0://height
			for (int w=0;w<Width;w++)
			{
				for (int d=0;d<Depth;d++)
				{
					Cells[d*Width*Height+x*Width+w]=Input[d][w];
				};
			};
			break;
		case 1://width
			for (int h=0;h<Height;h++)
			{
				for (int d=0;d<Depth;d++)
				{
					Cells[d*Width*Height+h*Width+x]=Input[h][d];
				};
			};
			break;
		case 2://depth
			for (int h=0;h<Height;h++)
			{
				for (int w=0;w<Width;w++)
				{
					Cells[x*Width*Height+h*Width+w]=Input[h][w];
				};
			};
			break;
		default:
			break;
		};
	};

	void InsertVector(MATRIX<T> &Input, int x, int y, int dimension)
	{
		switch (dimension)
		{
		case 0:
			for (int h=0;h<Height;h++)
			{
				Cells[y*Width*Height+h*Width+x]=Input[0][h];
			};
			break;
		case 1:
			for (int w=0;w<Width;w++)
			{
				Cells[y*Width*Height+x*Width+w]=Input[0][w];
			};
			break;
		case 2:
			for (int d=0;d<Depth;d++)
			{
				Cells[d*Width*Height+x*Width+y]=Input[0][d];
			};
			break;
		default:
			break;
		};
	};
};