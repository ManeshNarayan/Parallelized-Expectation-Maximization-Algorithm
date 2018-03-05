// At the moment the matrix operations are working fine. Below concerns need to be addressed.
// Not urgent.



// fgetm has been commented out right now. Meant for writing matrix to file.
// The tag "//?//" is used to to look out for changes that may need to be reviewd. Primarily to root out the verify calls
// The tag "//////" is used for the major edits done to the code by commenting out or including new code.
// Check the issue in destructor ~matrix. Removed delete for now
// For some reason all matrix allocations need to be carried out before other matrix operations. Need to fix this error
		// Handled due to the allocate function called initially in EM.
// Matrix mult is failing due to allocation of memory to matrix c in mult function.
		// Temp fix. Requires matrix to be allocated prior to calling mult on the matrix.
// Main currently used in this file is meant for testing only. To be removed later.
// Locking using afxmt has been removed. Can be addressed while parallelizing using private variable for c if the matrix c allocation is included in the mult function (removed in fix. Check above)







#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "mat.h"

//////
/*
// //verify will allow you to check for errors in the matrix operations
// I've undefined it for improved performance

#ifdef DEBUG
	#define //verify(a){  if (!a)   printf("ERROR"); exit(-1); }
#else 
	#define //verify(a)
#endif
*/



/*
///////
// The matrix class provides naive methods for matrix manipulation 
class matrix
{
protected:   
   int rows, cols, max_row, max_col;
   double *cells;
public:
   
   matrix ();
   ~matrix ();
   //matrix ( matrix &a );   
   //matrix (int x, int y);   
   matrix &operator=(matrix &a);               // perform a deep copy of the data in matrix a

   void allocate(int x, int y);            // create memory for x rows and y columns
   double * allocate(int x, int y, double * memory);
   void cleanup();                         // free the allocated memory
 
   void init();                              // zero out matrix structure 
   void ident(int size);                     // create an identity matrix 
   
   void add( matrix &a,  matrix &b);         // c = a + b
   void sub( matrix &a,  matrix &b);         // c = a - b
   void mult( matrix &a,  matrix &b);        // c = a * b 
   void trans( matrix &a);                   // c = a`        transpose a 
   double trace() ;                            // trace a 
   void sweep(matrix &a, int k);             // c = sweep(a)  sweep operator
   void inv( matrix &a);                     // c = a^-1      inverse a 
   void scamult(double scalar,  matrix &a);  // c = scalar * a 
   double determinant();
   double determinantOfMatrix(int n);
   void getCofactor(matrix temp, int p, int q, int n);
   double & pos(int row, int col);             // data at a(row, col)
   double * getAddr() { return cells;}
   // Made static to allow invoking without class object
   static double absEucDist(matrix a, matrix b);   // returns the absolute value of Euclidean distance

   void printm();                                        // print on screen(stdout)
   void fprintm(FILE * outfile);                         // print to outfile 
//?//   void   fgetm(FILE *infile, int rows, int cols); // read from ascii file 
   inline int get_rows() { return rows;}   
};
*/

// returns a matrix structure initialized to all 0's 
void matrix::init()
{
   int i, j;
   rows = max_row;
   cols = max_col;	
   cells[max_row*max_col] = 0;
   for(i = 0; i < max_row; i++)
      for(j = 0; j < max_col; j++)
	     pos(i,j) = 0;   
      	//cells[i*max_col + j];
}



// returns a size x size identity matrix 
void matrix::ident(int size)
{
   int i, j;
   //verify(max_row >= size);
   //verify(max_col >= size);
   cells[max_row*max_col] = 0;
   rows = size;
   cols = size;
   for(i = 0; i < size; i++)
   {
	   for(j = 0; j < size; j++)
		   pos(i,j) = 0;
	   pos(i,i) = 1;
   }
}


// returns a matrix structure multiplied by a scalar 
void matrix::scamult(double scalar,  matrix &a)
{
	int i, j;	
	rows = a.rows;
	cols = a.cols;
	for(i = 0; i < a.rows; i++)
		for(j = 0; j < a.cols; j++)
			pos(i,j) = scalar * a.pos(i,j) ;	

}




// adds to matrixes (this = a + b)
void matrix::add( matrix &a,  matrix &b)
{
   int i, j;	
   //verify(((a.rows + a.cols > 0) || (b.rows + b.cols > 0)));   
   if(a.rows!=b.rows || a.cols!=b.cols){
		printf("Error : Matrix dimensions don't match\n");
		return;
	}
	rows = a.rows;
	cols = a.cols;
      
  	//verify (((a.rows == b.rows) && (a.cols == b.cols)));
   
	for(i = 0; i < a.rows; i++)
		for(j = 0; j < a.cols; j++)
			pos(i,j) = a.pos(i,j) + b.pos(i,j);
}



// the difference of 2 matrixes (this = a - b)
void matrix::sub( matrix &a,  matrix &b)
{
	int i, j;	
	//verify(a.rows == b.rows);
   //verify(a.cols == b.cols);
	if(a.rows!=b.rows || a.cols!=b.cols){
		printf("Error : Matrix dimensions don't match\n");
		return;
	}
	rows = a.rows;
	cols = a.cols;
	for(i = 0; i < a.rows; i++)
		for(j = 0; j < a.cols; j++)
			//?// May need to make it "b.pos(i,j%b.cols)" to allow subtraction from all points
			pos(i,j) = a.pos(i,j) - b.pos(i,j);
}





//matrix c;
// If this code is run using threads (-mpi_nsot) the 
// commented lines using the variable "Lock" and "single"
// need to be included to prevent the threads from
// contending over the shared memory in the global variable ::c



////// Can make c a local var. No need to lock

//?//#include <afxmt.h>
//?//CCriticalSection Lock;

//?//
/*
// Mulitplies to matrixes (this = a * b)
void matrix ::mult( matrix &a,  matrix &b)
{
	//?//CSingleLock single(&Lock);
	int i, j, k;
	//?//single.Lock();
	//verify(a.cols == b.rows);
	
	////// Can make c a local var. No need to lock	
	c.init();
	c.rows = a.rows;
	c.cols = b.cols;
	for(i = 0; i < a.rows; i++)
		for(j = 0; j < b.cols; j++)
			for(k=0; k < a.cols; k++)
			{
				c.pos(i,j) += (a.pos(i,k) * b.pos(k,j));
				printf("%d %d\n",i,j);
			}

	*this = c;
	//?//single.Unlock();
}
*/

////// This fix is to remove the need to use a separate matrix c
// Mulitplies to matrixes (this = a * b)
void matrix::mult( matrix &a,  matrix &b)
{
	//?//CSingleLock single(&Lock);
	int i, j, k;
	//?//single.Lock();
	//verify(a.cols == b.rows);
	
	////// Can make c a local var. No need to lock	
	//c.init();
	if(a.cols!=b.rows){
		printf("Error matmult\n");
		return;
	}
	rows = a.rows;
	cols = b.cols;
	init();
	for(i = 0; i < a.rows; i++)
		for(j = 0; j < b.cols; j++)
			for(k=0; k < a.cols; k++)
			{
				pos(i,j) += (a.pos(i,k) * b.pos(k,j));
			}

	//?//single.Unlock();
}


// matrix transpose
void matrix::trans( matrix &a )
{
	int i, j;	
	
	rows = a.cols;
	cols = a.rows;
	//verify((a.rows <= max_col) || (a.cols <= max_row));
	for(i = 0; i < a.rows; i++)
		for(j = 0; j < a.cols; j++)
			pos(j,i) = a.pos(i,j);
}



// matrix trace
double matrix::trace() 
{
	int i;
	double trace;   

	//verify(rows == cols);
	
	for(i = 0, trace = 0.0; i < rows; i++)
   {      
		trace += pos(i,i);
   }
	return trace;
}


// matrix sweep
void matrix :: sweep(matrix &a, int k)
{
	double  d, b;
	int i, j;

	rows = a.rows;
	cols = a.cols;		
	init();

   d = a.pos(k,k);
	for(j = 0; j < a.rows; j++)
			a.pos(k,j) = a.pos(k,j)/d;
	for(i = 0; i < a.rows; i++)
	{
		b = a.pos(i,k);
		if( i!=k)
			for(j = 0; j < a.rows; j++)
				if(j==k)
						pos(i,j) = -b/d;
				else
						pos(i,j) = a.pos(i,j) - b *a.pos(k,j);
	}
}


					
// inverse			
void matrix::inv( matrix &a)
/* This method assumes a is a positive-definite, symmetric matrix.
	Returns the inverse of matrix a useing the sweep operator.
*/
{
	double  d, b;
	int i, j, k;
		
	//verify(a.rows == a.cols);	
	*this = a;
	
	for(k = 0; k < rows; k++)
	{
		d = pos(k,k);
		for(j = 0; j < rows; j++)
			pos(k,j) = pos(k,j)/d;
		for(i = 0; i < rows; i++)
			if( i!=k)
			{
				b = pos(i,k);
				for(j = 0; j < rows; j++)
					if(j !=k )
						pos(i,j) -= b * pos(k,j);
				pos(i,k) = -b/d;
			}
		pos(k,k) = 1.0/d;
	}
}



// get the data at the rowth row and colth column
double &matrix::pos(int row, int col)
{
	//verify( (row <= rows) && (col <= cols) );	
	//verify(cells[max_row*max_col] == 0);
   return cells[row * cols + col];
}



// prints a matrix to stdout */
void matrix::printm()
{
	fprintm(stdout);
}

// prints a matrix to outfile
void matrix::fprintm(FILE * outfile)
{
	int i, j;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
				//fprintf(outfile, "%20.10f", pos(i,j));
				fprintf(outfile, "%lf", pos(i,j));
		}
		
		fprintf(outfile, "\n");
	}
   fprintf(outfile, "\n");
}


// the constructor
matrix::matrix()
{
   cells = NULL;
}


// performs a deep copy
// (otherwise only the pointer is copied
// may lead to serious problems
matrix &matrix::operator=(matrix &a)
{
   if (&a==this)
      return a;

   //verify(max_row >= a.rows); // //verify(max_row * max_col >= a.rows * a.cols);
   //verify(max_col >= a.cols);

   rows = a.rows;
   cols = a.cols;
   for (int i = 0; i < a.rows; i++)
      for (int j = 0; j < a.cols; j++)
         pos(i,j) = a.pos(i,j);
   return *this;
}

// destructor
matrix::~matrix()
{
   if (cells != NULL) cleanup();
}

// frees memory
void matrix::cleanup()
{
   //?//delete cells;
	cells = NULL;
}

// allocate memory as needed
void matrix::allocate( int row, int col )
{
   cells = new double[row * col];
   max_row = row;
   max_col = col;
   rows = row; 
   cols = col;
}



// allocate memory from a pre-defined buffer
double * matrix::allocate( int row, int col , double * mem)
{
   cells = mem;
   max_row = row;
   max_col = col;
   rows = row; 
   cols = col;
   return mem + row * col;
}


void matrix::getCofactor(matrix temp, int p, int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp.pos(i,j++) = pos(row,col);
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}


double matrix::determinantOfMatrix(int n)
{
    double D = 0; // Initialize result
 
    //  Base case : if matrix contains single element
    if (n == 1)
        return pos(0,0);
 
    matrix temp; // To store cofactors
    double * addr = new double[n*n];
    temp.allocate(n,n,addr);
 
    int sign = 1;  // To store sign multiplier
 
     // Iterate for each element of first row
    for (int i = 0; i < n; i++)
    {
        // Getting Cofactor of mat[0][i]
        getCofactor(temp, 0, i, n);
        D += sign * pos(0,i) * temp.determinantOfMatrix(n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}

double matrix::determinant()
{
	if(rows == cols)
		return determinantOfMatrix(rows);
	else
	{
		printf("Error : Not square matrix\n");
		return (double)0;
	}
}

double matrix::absEucDist(matrix a, matrix b)
{

	if(a.cols != 1 || b.cols != 1)
	{
		printf("Error : Euclidian distance attempted on non vector\n");
		return 0;
	}
	if(a.rows != b.rows)
	{
		printf("Error : Euclidean distance attempted on incompatible vectors\n");
		return 0;
	}

	double val = 0;
	double temp;

	for(int i = 0; i < a.rows; i++)
	{
		temp = a.pos(i,0) - b.pos(i,0);
		val += temp*temp;
	}

	val = sqrt(val);
	return val;
}



/* reads a matrix from the next rows rows of infile(an ascii file).
	there should be cols numbers in each row.
	infile is left pointing to the next row in the file to allow several
	matrices to be input in a single file or to allow access to several
	rows of a matrix.
*/
	//////
	/*
void matrix ::fgetm(FILE *infile, int row, int col)
{
	int i, j;	
   static char * str = (sizeof(double)==8)? "%lf":"%f";
	rows = row;
	cols = col;
	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
      {
			fscanf(infile, str, &pos(i,j));
      }
		fscanf(infile,"\n");
	}
   //verify(cells[max_row*max_col] == 0);
}
*/


/*

////// This main function is just meant for testing. Can be removed later.
int main(){
	printf("hello\n");
	matrix a,b,c,d;
	int memReq = 0;
	memReq = 3*3 + 3*3 + 1*3 +3*3;
	double * memspace = new double[memReq];
	printf("%p\n",(void*)memspace);
	memspace = a.allocate(3,3,memspace);
	printf("%p\n",a.getAddr());
	memspace = c.allocate(3,3,memspace);
	printf("%p\n",(void*)memspace);
	memspace = b.allocate(1,3,memspace);
	printf("%p\n",(void*)memspace);
	memspace = d.allocate(3,1,memspace);
	printf("%p\n",(void*)memspace);
	a.pos(0,0) = -2;
	a.pos(0,1) = 2;
	a.pos(0,2) = -3;
	a.pos(1,0) = 0;
	a.pos(1,1) = 2;
	a.pos(1,2) = -4;
	a.pos(2,0) = 0;
	a.pos(2,1) = 0;
	a.pos(2,2) = 4.5;
	printf("%lf\n",a.determinant());
	b.pos(0,0) = -1;
	b.pos(0,1) = 1;
	b.pos(0,2) = 0;
	d.pos(0,0) = -1;
	d.pos(1,0) = 1;
	d.pos(2,0) = 0;
	c.mult(b,a);
	a.mult(c,d);
	//c.printm();
	return 0;
}
*/
