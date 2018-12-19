#include <mpi.h>
#include <math.h>
#include <iostream>
#include <chrono>


#ifdef __cplusplus
extern "C" {
#endif

  int check2DHeat(double** H, long n, long rank, long P, long k); //this assumes array of array and grid block decomposition

#ifdef __cplusplus
}
#endif

//Use similarily as the genA, genx from matmult assignment.
double genH0(long row, long col, long n) {
  double val = (double)(col == (n/2));
  return val;
}

//get the overlapping values between processes needed for computation
void communicate(double** H_prev, long row, long col, int pid, long size, long blocks, double *top, double *bottom, double* left, double *right, int nbprocess) {
 	MPI_Request *request;
	MPI_Status *status;
  
	for (long i=0; i<size; i++) {
			left[i] = H_prev[i][0];
			right[i] = H_prev[i][size-1];
		for (long j=0; j<size; j++) {
			top[j] = H_prev[0][j];
			bottom[j] = H_prev[size-1][j];
			}
		}
	if (nbprocess > 1) {
	if (row == 0 && col == 0) {				
		request  = new MPI_Request[4];
  		status = new MPI_Status[4];
  
  		MPI_Irecv(top, size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(left, size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Isend(bottom, size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[2]);;
  		MPI_Isend(right, size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[3]);
  		MPI_Waitall(4, request, status); }

	else if (row == 0 && col==(blocks-1)) {
		request  = new MPI_Request[4];
		status = new MPI_Status[4];
		
		MPI_Irecv(top, size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(right, size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Isend(bottom, size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[2]);;
  		MPI_Isend(left, size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[3]);
  		MPI_Waitall(4, request, status); }
			
	else if (row == 0) {
		request  = new MPI_Request[6];
  		status = new MPI_Status[6];
  	
		MPI_Irecv(top,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(left,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Irecv(right,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[2]);
  		MPI_Isend(bottom,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[3]);;
  		MPI_Isend(right,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[4]);
  		MPI_Isend(left,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[5]);
  		MPI_Waitall(6, request, status); }
			
	else if (row == blocks-1 && col ==0) {
		request  = new MPI_Request[4];
		status = new MPI_Status[4];

		MPI_Irecv(bottom,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(left,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[1]);
		MPI_Isend(top,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[2]);;
		MPI_Isend(right,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[3]);
		MPI_Waitall(4, request, status); }
  				
	else if (row == blocks-1 && col == blocks-1) {
 		request  = new MPI_Request[4];
  		status = new MPI_Status[4];

  		MPI_Irecv(bottom,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(right,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Isend(top,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[2]);;
  		MPI_Isend(left,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[3]);
  		MPI_Waitall(4, request, status); }
  				
	else if (row == blocks-1) {
 		request  = new MPI_Request[6];
  		status = new MPI_Status[6];
  
  		MPI_Irecv(bottom,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(left,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Irecv(right,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[2]);
  		MPI_Isend(top,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[3]);;
  		MPI_Isend(right,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[4]);
  		MPI_Isend(left,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[5]);
  		MPI_Waitall(6, request, status); }
			
	else if (col == 0) {
 		request  = new MPI_Request[6];
  		status = new MPI_Status[6];

		MPI_Irecv(top,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(bottom,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Irecv(left,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[2]);
  		MPI_Isend(bottom,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[3]);;
  		MPI_Isend(right,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[4]);
  		MPI_Isend(top,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[5]);
  		MPI_Waitall(6, request, status); }
  				
  	else if (col == blocks-1) {
  		request  = new MPI_Request[6];
  		status = new MPI_Status[6];
  	
  		MPI_Irecv(top,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(bottom,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Irecv(right,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[2]);
  		MPI_Isend(bottom,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[3]);;
  		MPI_Isend(left,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[4]);
  		MPI_Isend(top,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[5]);
  		MPI_Waitall(6, request, status); }
  				
  	else {
   		request  = new MPI_Request[8];
  		status = new MPI_Status[8];
  	
  		MPI_Irecv(top, size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[0]);
  		MPI_Irecv(bottom,size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[1]);
  		MPI_Irecv(left,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[2]);
  		MPI_Irecv(right,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[3]);
  		MPI_Isend(bottom,size, MPI_DOUBLE, pid+blocks, 0, MPI_COMM_WORLD, &request[4]);
  		MPI_Isend(right,size, MPI_DOUBLE, pid+1, 0, MPI_COMM_WORLD, &request[5]);
  		MPI_Isend(left,size, MPI_DOUBLE, pid-1, 0, MPI_COMM_WORLD, &request[6]);
  		MPI_Isend(top, size, MPI_DOUBLE, pid-blocks, 0, MPI_COMM_WORLD, &request[7]);
  		MPI_Waitall(8, request, status); } }
} // end communicate()


//Calculate the heat
double compute(double **H_prev, double **H, long size, double *top, double *bottom, double *left, double *right) {
	double k = 1/static_cast<double>(5);

	for(int i=0;i<size;i++){
      for(int j=0;j<size;j++){

        if(i==0 && j==0){
          H[i][j]=(k)*(H_prev[i+1][j]+H_prev[i][j+1]+H_prev[i][j]+right[i]+bottom[j]); } //top-left corner element
          
        else if(i==0 && j==size-1){
          H[i][j]=(k)*(H_prev[i][j-1]+H_prev[i][j+1]+H_prev[i][j]+left[i]+bottom[j]); } //top-right corner element
          
        else if(i==0){
          H[i][j]=(k)*(H_prev[i][j-1]+H_prev[i+1][j]+H_prev[i][j+1]+bottom[j]+H_prev[i][j]); } //other elements in top row
        
        //top row computation done
        
        else if(i==size-1 && j==0){
          H[i][j]=(k)*(H_prev[i][j+1]+H_prev[i-1][j]+H_prev[i][j]+right[i]+top[j]); } //bottom-left corner element
          
        else if(i==size-1 && j==size-1){
          H[i][j]=(k)*(H_prev[i][j-1]+H_prev[i-1][j]+H_prev[i][j]+left[i]+top[j]); } //bottom-right corner element
          
        else if(i==size-1){
          H[i][j]=(k)*(H_prev[i][j-1]+H_prev[i-1][j]+H_prev[i][j+1]+top[j]+H_prev[i][j]); } //other elements in bottom row
          
        //bottom row computation done
          
        else if(j==0){
          H[i][j]=(k)*(H_prev[i][j+1]+H_prev[i+1][j]+H_prev[i-1][j]+right[i]+H_prev[i][j]); } //other elements in leftmost column
                   
        else if(j==size-1){
          H[i][j]=(k)*(H_prev[i][j-1]+H_prev[i+1][j]+H_prev[i-1][j]+left[i]+H_prev[i][j]); } //other elements in rightmost column
          
        else{
          H[i][j]=(k)*(H_prev[i-1][j]+H_prev[i][j-1]+H_prev[i][j]+H_prev[i+1][j]+H_prev[i][j+1]); } //core elements
      
    }
 }
}// end compute()
    
int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <N> <K>"<<std::endl;
    return -1; }

  // declare and init command line params
  long N, K;
  N = atol(argv[1]);
  K = atol(argv[2]);
  
  //MPI Initialization
  MPI_Init(NULL, NULL);

  //getting number of processes 
  int nbprocess;
  MPI_Comm_size(MPI_COMM_WORLD, &nbprocess);

  //getting the rank of process
  int pid;
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  
  //Block init
  long blocks = sqrt(nbprocess);
  long size = N/blocks;
  long row = pid/blocks;
  long col = pid%blocks;
  long row_start = row * size;
  long col_start = col * size;
  long row_end = row_start+size;
  long col_end = col_start+size;
  
  // use double for heat 2d 
  double** H = new double *[size];
  double** H_prev = new double *[size];
  
  //intializing matrix
  for (long i=row_start; i<row_end; i++) {
  	H_prev[i-row_start] = new double[size];
  	H[i-row_start] = new double[size];
  	for (long j=col_start; j<col_end; j++) {
  		H_prev[i-row_start][j-col_start] = genH0(i,j,N);
  		H[i-row_start][j-col_start] = 0; }	}
		
  //uncomment below to print the block sized-matrices
  /*for(int i=0;i<size;i++){
  	for(int j=0;j<size;j++){
  		std::cout<<H_prev[i][j]<<" ";}
  	std::cout<<pid<<std::endl;}*/
  	
  double *top = new double[size];
  double *bottom = new double[size];
  double *right = new double[size];
  double *left = new double[size];
  	
  //start time
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  
  MPI_Barrier(MPI_COMM_WORLD);
  //write code here
  for (long iter=0; iter<K; iter++) {
  	double **temp;
    communicate(H_prev, row, col, pid, size, blocks, top, bottom, left, right, nbprocess);
    compute(H_prev, H, size, top, bottom, left, right);
    
    /*if (check2DHeat(H, N, pid, nbprocess, iter)) {
     	//std::cout<<"oh no it is (maybe) wrong"<<std::endl; 
     std::cout<<"rank: "<<pid<<" is incorrect"<<std::endl;
   }*/
    
  	temp = H_prev;
  	H_prev = H;
  	H = temp;
  }

  //end time and evaluating the duration of computation
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  if (pid ==0) {
  		std::cerr<<elapsed_seconds.count()<<std::endl; }
  MPI_Finalize();
  
  //freeing matrix
  for (long i =0; i< size; i++) {
  	delete[] H[i];
  	delete[] H_prev[i];}
  delete[] H;
  delete[] H_prev;
  delete[] top;
  delete[] bottom;
  delete[] right;
  delete[] left;
  return 0;
}
