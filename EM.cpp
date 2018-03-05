// Mat sizes in comments are correct. Verify allocation based on column vector
// Need to add functionality to allocate weight as 1 or 0 for initial few steps. To do later. Optimization according to armadillo
   // approximates to K means allowing faster computations to get start locations.

// I - no of points
// D - no of dim/attr
// C - no of sources/gaussians
// p - no of mpi nodes


//#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <mpi.h>
#include <omp.h>
//#include "mat.h"
#include "EM.h"
#include <ctime>
#include <chrono>
#define  MINCHANGE 0.001
#define numThreads 64
using namespace std;


bool converged;
int trial;

//?// Make changes here after coding individual functions
void Expectation_Maximization::run()
{  
	
   //omp_set_num_threads(4);
   //int a = omp_get_num_threads();
   //cout<<"Num of threads is"<<a<<endl;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &totalProcs);
   getInitData();
   chrono::system_clock::time_point start;
   chrono::system_clock::time_point end;
   if(myrank == 0)
   {
   	start = chrono::system_clock::now();
   }
   double centroidShift = 0;
   #pragma omp parallel num_threads(numThreads) shared(trial,converged)
   {  //cout<<converged;

      matrix covar_inv;
      covar_inv.allocate(D,D);
      matrix vec;
      matrix vec2;
      matrix vecT;
      matrix vecT_covarInv;
      matrix vecT_covarInv_vec;
      matrix red_covar;
      vec.allocate(D,1);
      vec2.allocate(D,1);
      vecT.allocate(1,D);
      vecT_covarInv.allocate(1,D);
      vecT_covarInv_vec.allocate(1,1);
      red_covar.allocate(D,D); 

      while(!converged&&trial<20){
         //#pragma omp barrier
         //printCluster();
         if(myrank == 0){
            #pragma omp master
            {
               cout<<"**************Start trial #"<<trial+1<<"**************"<<endl;
            }
         }
         //#pragma omp barrier

         // Pass_1
         {  
            //red_Prob_xi_c to be zeroed out in later pass, as efficiency necessitates C in outer loop and I_loc in inner
            
            for(int i = 0; i < C; i++)
            {  

               // Error check for determinant
               double det;
         //      #pragma omp critical
         //      {
                  det = covar[i].determinant();
         //      }

               if(det<=0)
               {  
                  #pragma omp master
                  {
                     cout<<"non positive det in clus "<<i<<"with val = "<<det<<endl;   
                  }
                  MPI_Abort(MPI_COMM_WORLD,1);
                  exit(0);
               }
               
               covar_inv.inv(covar[i]);

               #pragma omp for
               for(int j = 0; j < I_loc; j++)
               {

                  // Matrix computations
                  vec.sub(pts[j],mean[i]);
                  vecT.trans(vec);
                  vecT_covarInv.mult(vecT,covar_inv);
                  vecT_covarInv_vec.mult(vecT_covarInv,vec);

                  // Value computation
                  double val = vecT_covarInv_vec.pos(0,0);
                  val = (val/2)*(-1);
                  val = exp(val);
                  val = val/(sqrt(pow(2*M_PI,D)*det));

                  // Result store
                  Prob_xi_c.pos(j,i) = val;
                    
                  // Reduction for next pass // No need for all reduce
                  red_Prob_xi_c[j] += val*Prob_c[i];         

               }
               
            }

         }



         // Need to do trial increment here, otherwise extra barrier needed at while entry to ensure all threads enter
         #pragma omp single
         {
            trial++;
         }



         //#pragma omp barrier
         // pass_2();
         {
            #pragma omp for
            for(int i = 0; i < C; i++)
            {
               //loc_red_Prob_c_xi[i] = 0;
               red_Prob_c_xi[i] = 0;
               for(int j = 0; j < I_loc; j++)
               {
                  Prob_c_xi.pos(i,j) = (Prob_xi_c.pos(j,i)*Prob_c[i])/red_Prob_xi_c[j];

                  // Reduction for next pass // Need all reduce
                  red_Prob_c_xi[i] += Prob_c_xi.pos(i,j);
               }
               
            }

            // Zero out for next iteration
            //#pragma omp master
            #pragma omp for
            for(int j = 0; j < I_loc; j++)
            {
               red_Prob_xi_c[j] = 0;
            }   

         }





         //#pragma omp barrier
         // pass_3();
         {
            //#pragma omp parallel num_threads(numThreads)
            {

               #pragma omp for
               for(int i = 0; i < C; i++)
               {
                  vec.init();
                  //vec.printm();
                  for(int j = 0; j < I_loc; j++){

                     // Computes importance of point xi for cluster C
                     w_i_c.pos(j,i) = Prob_c_xi.pos(i,j)/red_Prob_c_xi[i];
                     //cout<< w_i_c.pos(j,i)<<"\t"<<endl;
                     //pts[j].printm();

                     // Compute contribution of a point on node to mean of cluster C
                     vec2.scamult(w_i_c.pos(j,i),pts[j]);
                     //vec2.printm();
                     // Sums contributions for all points on the node to mean of cluster
                     vec.add(vec,vec2);
                     //vec.printm();
                     //cout<<endl<<endl;
                  }
                  
                  // Value store
                  mean[i] = vec;
                  //mean[i].printm();

                  // Computes priors for next iteration
                  Prob_c[i] = red_Prob_c_xi[i]/I;
                  //cout<<Prob_c[i]<<endl;
                     
               }
               
            }
         }




         //#pragma omp barrier
         // pass_4();
         {
            #pragma omp for
            for(int i = 0; i < C; i++){
               covar[i].init();
            }
            #pragma omp for
            for(int i = 0; i < C; i++)
            {
               //#pragma omp barrier
              //#pragma omp parallel num_threads(numThreads)
               //{
        

               //#pragma omp for
               for(int j = 0; j < I_loc; j++)
               {
                  vec.sub(pts[j], mean[i]);
                  vecT.trans(vec);
                  vec.scamult(w_i_c.pos(j, i), vec);
                  red_covar.mult(vec,vecT);
                  //#pragma omp critical
                  //{
                     covar[i].add(covar[i],red_covar);
                  //}
                  //#pragma omp barrier
               }
               //}
            }
            
         }



         //#pragma omp barrier
         //MPI_Allreduce(MPI_IN_PLACE, covar[0].getAddr(), C*D*D, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         //#pragma omp single
         {
            //pass_5();
            {

               //double centroidShift = 0;
               if(trial!=0)
               {

                  #pragma omp for reduction(+:centroidShift)
                  for(int i = 0; i < C; i++)
                  {
                     centroidShift += matrix::absEucDist(mean[i],mean_old[i]);
                     mean_old[i] = mean[i];
                  }
               }
               

               #pragma omp single
               {
                  if(myrank == 0)
                  {
                     cout<<"CentroidShift = "<<centroidShift<<endl;
                  }

                  if(centroidShift>epsilon||trial==0)
                  {
                     converged = false;
                  }
                  else{
                     converged = true;
                  }
                  centroidShift = 0;
               }
               
               
            }

         }
      }
      
   }
   //printCluster();
   //cout<<"Hello"<<endl<<endl<<endl<<"******"<<endl;
   if(myrank == 0)
   {
	   end = chrono::system_clock::now();
      chrono::duration<double> elapsed = end-start;
   	cout<<endl<<"Time Taken = "<<elapsed.count();
   }
   
}


void Expectation_Maximization::getInitData(){

   double* initParam = new double[5];
   double * pts_temp = NULL;
   double maxval;
   double minval;

   if(myrank == 0)
   {
      cin>>I;     // I value or no of pts
      cin>>D;     // D value or dimensions per pt
      cin>>C;     // C value or no of cluster

      I_loc_list = new int[totalProcs];
      Disp = new int[totalProcs];
      I_loc_list[0] = loc_count(0)*D;
      Disp[0] = 0;
      for(int i = 1; i < totalProcs; i++)
      {
         I_loc_list[i] = loc_count(i)*D;
         Disp[i] = I_loc_list[i-1] + Disp[i-1];
      }

      minval = INT_MAX;
      maxval = INT_MIN;
      pts_temp = new double[I*D];
      for(int i=0; i<I; i++)
      {
         for(int j=0; j<D; j++)
         {
            cin>>pts_temp[i*D + j];   // Column matrix
            if(pts_temp[i*D + j] > maxval)
            {
               maxval = pts_temp[i*D + j];
            }
            if(pts_temp[i*D + j] < minval)
            {
               minval = pts_temp[i*D + j];
            }
         }
      }
      initParam[0] = I;
      initParam[1] = D;
      initParam[2] = C;
      initParam[3] = maxval;
      initParam[4] = minval;      
   }

   // Sends parameters for allocations
   MPI_Bcast(initParam, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
   I = (int)initParam[0];
   D = (int)initParam[1];
   C = (int)initParam[2];
   I_loc = loc_count(myrank);
   maxval = initParam[3];
   minval = initParam[4];

   long long memReq = 0;
   double * memSpace;

   //I_loc = I;

   ////// When using mpi the value of I will have to be changed for just pts on system for non proc0

   w_i_c.allocate(I_loc,C);
   Prob_xi_c.allocate(I_loc,C);
   Prob_c_xi.allocate(C,I_loc);
   //covar_inv.allocate(D,D);
   //red_covar.allocate(D,D);

   Prob_c = new double[C];

   memReq = I_loc*D;
   memSpace = new double[memReq];
   pts = new matrix[I_loc];
   for(int i = 0; i < I_loc; i++){
      memSpace = pts[i].allocate(D,1,memSpace);
   }
   
   memReq = C*D;          // Weighted mean (mean)
   memSpace = new double[memReq];
   mean = new matrix[C];
   mean_old = new matrix[C];
   for(int i = 0; i < C; i++){
      memSpace = mean[i].allocate(D,1,memSpace);
      mean_old[i].allocate(D,1);
   }

   memReq = C*D*D;        // Covariance values (covar)
   memSpace = new double[memReq];
   covar = new matrix[C];
   for(int i=0; i<C; i++){
      memSpace = covar[i].allocate(D,D,memSpace);
      // Initialize values
      covar[i].ident(D);
   }

   // Arrays used to hold the summed value of the probabilities over all xi
   red_Prob_c_xi = new double[C];
   //loc_red_Prob_c_xi = new double[C];
   red_Prob_xi_c = new double[I_loc];


   for(int i=0; i<C; i++){
      // Initialize values
      Prob_c[i]= 1/(double)C;
      red_Prob_c_xi[i] = 0;
   }

   for(int i = 0; i < I_loc; i++){
      red_Prob_xi_c[i] = 0;
   }

   
   
   /*
   // Initialize mean to same random values on all machines.
   memSpace = mean[0].getAddr();
   srand(1);
   double range = maxval - minval;
   
   //?// remove later. Kept for testing 
   if(range == 0){
      cout<<endl<<"DAMN"<<endl;
   }

   for(int i=0; i<C*D; i++){
      memSpace[i] = (double)(rand()%(int)ceil(range)) + minval;
   }
   */


   MPI_Scatterv(pts_temp, I_loc_list, Disp, MPI_DOUBLE, pts[0].getAddr(), I_loc*D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   delete pts_temp;

   if(myrank == 0)
   {
      for(int i = 0; i < C; i++)
      {
         mean[i] = pts[i];
      }
   }

   MPI_Bcast(mean[0].getAddr(), C*D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   //MPI_Bcast(initParam, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   trial = -1;
   //loglik = 0;
   epsilon = 0.001;
   converged = false;  
   //* converged = false;
}

// P(xi|C) = 1/(sqrt((2pi^D)*det(covar)))*exp((-1/2)*(xi-mc)T*(covar_inv)*(xi-mc))
//void Expectation_Maximization::pass_1()



// P(C|xi) = P(xi|C)P(C)/sumAcrossClusters(P(xi|C)P(C))
//void Expectation_Maximization::pass_2()


// w_i_c = P(C|xi)/sumAcrossPoints(P(C|xi))
// P(C) = sumAcrossPoints(P(C|xi))/I
//void Expectation_Maximization::pass_3()



// covar(c,j,k) = sumAcrossPoints(w_i_c*(xi_j - mean[c]_j)*(xi_k - mean[c]_k))
//void Expectation_Maximization::pass_4()


/*
// Check for convergence using loglikelihood formula
void Expectation_Maximization::pass_5(){
   loglik = 0;
   for(int j = 0; j < I_loc; j++){
      temp = 0;
      for(int i = 0; i < C; i++){
         temp += Prob_xi_c.pos(j,i)*Prob_c[i];
      }
      loglik += log(temp);
   }
//   cout<<loglik_old<<"  "<<loglik<<endl;
   //?// Need to call allreduce here on loglik
   if(trial == 0){
      converged = false;
   }
   else{
      if(fabs(loglik_old - loglik) > epsilon){
         converged = false;
      }
      else{
         converged = true;
      }
   }
   loglik_old = loglik;
}
*/

//void Expectation_Maximization::pass_5()


void Expectation_Maximization::printCluster()
{  
   ///*
   if(myrank == 0){
      for(int i = 0; i < C; i++)
      {
         cout<<"mean"<<i<<" = ";
         mean[i].printm();
      }
      for(int i = 0; i < C; i++)
      {
         cout<<"covar"<<i<<" = ";
         covar[i].printm();
      }
      for(int i = 0; i < C; i++)
      {
         cout<<"Prob"<<i<<" = "<<Prob_c[i]<<endl<<endl;
      }   
   }
   
   //*/
   for(int j = 0; j < I_loc; j++)
   {
      int clus = -1;
      double maxProb = -1;
      for(int i = 0; i < C; i++)
      {
         //cout<<"Point"<<j<<" Cluster"<<i<<" Prob"<<Prob_xi_c.pos(i,j)<<endl;
         
         if(Prob_c_xi.pos(i,j)>maxProb)
         {
            maxProb = Prob_c_xi.pos(i,j);
            clus = i;
         }
         //cout<<"Point"<<j<<" Cluster"<<i<<" Prob"<<Prob_c_xi.pos(i,j)<<endl;
         
      }
      cout<<"Point"<<j<<" Cluster"<<clus<<" Prob"<<maxProb<<endl;
   }
}

