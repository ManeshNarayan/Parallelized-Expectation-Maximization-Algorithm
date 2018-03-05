
#include "mat.h"

class Expectation_Maximization
{
private:

// important values
   int I;               // Holds the number of points
   int I_loc;           // Holds the number of points on local system
   int D;               // Holds the number of dimensions for each point
   int C;               // Holds the number of clusters to initialize
   int myrank;          // Holds rank of process
   int totalProcs;      // Holds number of porcs
   int * I_loc_list;    // Array holding number of points in each process
   int * Disp;          // Array holding displacement using which to send points

   // Matrices used to input values
   matrix * pts;          // holds the points. I*D matrix. Needs to hold only I/p.

   // Matrices used to hold gaussian values
   matrix * mean;       // holds the weighted mean of attribute j for source c. C*D matrix.
   matrix * covar;      // holds the covariance values for all attr pairs per source. C*D*D matrix.
   matrix Prob_xi_c;    // Probability of point xi given source c. I*C matrix. Needs to hold only I/p. Allocate reducibly for each source.
   matrix Prob_c_xi;    // Probability of source c given point xi. I*C matrix. Needs to hold only I/p. Allocate reducibly for each point.
   matrix w_i_c;        // matrix holds the importance/weight of a point i for source c. I*C matrix. Needs to hold only I/p.
   double * Prob_c;       // matrix holds the prior probability for source c.

   // Arrays used for reducing values
   double * red_Prob_c_xi;
   double * red_Prob_xi_c;

   // Variables to check convergence
   //int trial;
   //double loglik_old;
   //double loglik;
   double temp;
   //bool converged;
   double epsilon;


// internal use values
   //matrix covar_inv;  // Holds the inverse of covariance matrix.

// temporary variables
   /*
   matrix vec;
   matrix vec2;
   matrix vecT;
   matrix vecT_covarInv;
   matrix vecT_covarInv_vec;
   matrix red_covar;
   */

   
   matrix * mean_old;
   //double * loc_red_Prob_c_xi;





protected:
   
//   void pass_1();
//   void pass_2();
//   void pass_3();
//   void pass_4();
//   void pass_5();
   void getInitData();
   void printCluster();
   int loc_count(int procNum);
   //double pass_6();
   

// function to spread the data to all processors appropriately
   //void getInitData();


// functions to allocate memory


public:
/*
   Expectation_Maximization(char *, void*);   // Initializes the variables from the files
   ~Expectation_Maximization();               // cleans up memory
*/
   void run();                                // runs the EM algorithm

};


inline int Expectation_Maximization::loc_count(int procNum)
{
   int count = I / totalProcs;
   if ( procNum < I % totalProcs ) {
      count++;
   }
   return count;
}




// Global variables to be used by EM class
/*matrix vec;
matrix vec2;
matrix vecT;
matrix vecT_covarInv;
matrix vecT_covarInv_vec;
*/
//#pragma omp threadprivate(vec,vecT,vec2,vecT_covarInv,vecT_covarInv_vec)
