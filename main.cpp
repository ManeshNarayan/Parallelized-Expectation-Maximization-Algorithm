#include <iostream>
#include "EM.h"
#include <mpi.h>
#include <ctime>

using namespace std;

//int myrank;

int main(int argc, char* argv[]){

	//MPI_Init(&argc,&argv);
	int requested, provided;
	requested = MPI_THREAD_FUNNELED;
	MPI_Init_thread(&argc,&argv,requested,&provided);
	if(provided<requested){
		cout<<endl<<"Insufficient level of parallelism provided. Aborting"<<endl;
		MPI_Abort(MPI_COMM_WORLD,1);
		exit(0);
	}
	//MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	Expectation_Maximization A;
	A.run();
	MPI_Finalize();
   	return 0;
}
