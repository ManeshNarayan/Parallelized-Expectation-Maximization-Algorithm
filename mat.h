


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
 
   void init();                              // zero out matrix structure */
   void ident(int size);                     // create an identity matrix */
   
   void add( matrix &a,  matrix &b);         // c = a + b */
   void sub( matrix &a,  matrix &b);         // c = a - b */
   void mult( matrix &a,  matrix &b);        // c = a * b */
   void trans( matrix &a);                   // c = a`        transpose a */
   double trace() ;                            // trace a */
   void sweep(matrix &a, int k);             // c = sweep(a)  sweep operator */
   void inv( matrix &a);                     // c = a^-1      inverse a */
   void scamult(double scalar,  matrix &a);  // c = scalar * a */
   double determinant();
   double determinantOfMatrix(int n);
   void getCofactor(matrix temp, int p, int q, int n);
   double & pos(int row, int col);             // data at a(row, col) */
   double * getAddr() { return cells;}
   /* matrix i/o functions */


   /*vector functions*/
   // Made static to allow invoking without class object
   static double absEucDist(matrix a, matrix b);   // returns the absolute value of Euclidean distance

   void printm();                                        // print on screen(stdout) */
   void fprintm(FILE * outfile);                         // print to outfile */
//?//   void   fgetm(FILE *infile, int rows, int cols); // read from ascii file */
   inline int get_rows() { return rows;}   
};
