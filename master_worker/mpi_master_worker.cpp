#include <mpi.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <tuple>

#ifdef __cplusplus
extern "C" {
#endif

float f1(float x, int intensity);
float f2(float x, int intensity);
float f3(float x, int intensity);
float f4(float x, int intensity);

#ifdef __cplusplus
}
#endif

#define INITIAL_WORK_REQ 1
#define QUIT 1


//This function chooses the function based on the input arguments and returns the function result
float func_selector(int func_id, float x, int intensity) {

    switch (func_id) {
        case 1 : return f1(x, intensity);
        case 2 : return f2(x, intensity);
        case 3 : return f3(x, intensity);
        case 4 : return f4(x, intensity);
        default : return -1;

        }
    }

float computeIntegral(int start, int end, int function_id, int intensity, float a, float b, long n) {
	float result = 0.0;
	//std::cout<<"From CI: start "<<start<<"   end "<<end<<std::endl;
	float width = (b - a)/float(n);
	for(int i = start; i < end; i++) {
		float x = (a + (i + 0.5) * width);
        //To implement the function selection based on the input
        float func = func_selector(function_id, x, intensity);
        result = result + (width * func);
    }
    //std::cout<<"From CI: result "<<result<<std::endl;
    return result;
}

std::tuple<int, int> getData(int req_id, long n, int nbprocess) {
	nbprocess = nbprocess-1;
	int gran = n/(nbprocess);
	//std::cout<<"In GD "<<req_id<<std::endl;
	int start_ptr = req_id * gran;
	int end_ptr = start_ptr + gran;
	//std::cout<<start_ptr<<" "<<end_ptr<<std::endl;
	
	if ((n%nbprocess !=0) && (end_ptr > n)){
		end_ptr = n; }
	
	return std::make_tuple(start_ptr, end_ptr); 
}
	
float master(long n, int nbprocess) {
	int is_inital_req = 0;
	MPI_Status status;
	float final_result = 0.0;
	int req_id = -1;
	int work_sent = 0;
	int start, end=0;
	float result=0.0;
	
	//Send work to all processes
	for (int i=1; i<nbprocess; i++) {
		if (end < n) {
			req_id++;
			work_sent++;
			//std::cout<<req_id<<std::endl;
			//get start and end pointers
			std::tie(start, end) = getData(req_id, n, nbprocess);
			//std::cout<<"ID#"<<i<<" "<<start<<" "<<end<<std::endl;
			int work[2]={0};
			work[0] = start;
			work[1] = end;
			MPI_Send(work, 2, MPI_INT, i, 0, MPI_COMM_WORLD);
			//std::cout<<"sent"<<i<<std::endl; 
			}
			
		else {
			MPI_Send(0, 0, MPI_INT, i, QUIT, MPI_COMM_WORLD); }

		}
		
	while(work_sent != 0) {
		MPI_Status status;
		//Receive the result
		MPI_Recv(&result, 1, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		int id = status.MPI_SOURCE;
		//get result
		final_result += result;
		work_sent--;
		
		//check if there is any more work else send QUIT signal
		if (end < n) {
			req_id++;
			work_sent++;
			//get start and end pointers
			std::tie(start, end) = getData(req_id, n, nbprocess);
			//std::cout<<"ID#"<<id<<"start "<<start<<"   end "<<end<<std::endl;
			int work[2] = {0};
			work[0] = start;
			work[1] = end;
			MPI_Send(work, 2, MPI_INT, id, 0, MPI_COMM_WORLD); }
		
		else {
			MPI_Send(0, 0, MPI_INT, id, QUIT, MPI_COMM_WORLD); }
	
		}
			
	
	//std::cout<<"final result "<<final_result<<std::endl;
	//return the final result
	return final_result;
}

void worker(int function_id, int intensity, float a, float b, long n) {
	float result = 0.0;
	int work[2] = {0};
	MPI_Status status;
	while(1) {
		//Receive work if any from master process or receive quit signal
		MPI_Recv(work, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		int tag = status.MPI_TAG;
		if (tag != QUIT) { 
			int start = work[0];
			int end = work[1];
			//std::cout<<"From worker: start "<<start<<"   end "<<end<<std::endl;
			//compute the result
			result = computeIntegral(start, end, function_id, intensity, a, b, n);
	
			//Send computed result
			MPI_Send(&result, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD); }
		else {
		return; }
		}
	}

int main (int argc, char* argv[]) {

  if (argc < 6) {
    std::cerr<<"usage: mpirun "<<argv[0]<<" <functionid> <a> <b> <n> <intensity>"<<std::endl;
    return -1;
  }
  
   //Initializing the variables
    int function_id, intensity;
    long int n;
    float a, b;
    float result = 0.0;
    //MPI initialization
	MPI_Init(NULL, NULL);

    //getting number of processes into nbprocess
	int nbprocess;
	MPI_Comm_size(MPI_COMM_WORLD, &nbprocess);
    //getting the rank of process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //reading the input values into corresponding variables
    sscanf(argv[1], "%i", &function_id);
    sscanf(argv[2], "%f", &a);
    sscanf(argv[3], "%f", &b);
    sscanf(argv[4], "%ld", &n);
    sscanf(argv[5], "%i", &intensity);
    
    using namespace std::chrono;
    //high_resolution_clock: clock specified by the chrono library for evaluating the time 
    //time_point: relative measurement of time
    high_resolution_clock::time_point start_time = high_resolution_clock::now();
    
    if (rank == 0) {
    	result = master(n, nbprocess);
    	//std::cout<<result<<std::endl;
    	}
    else {
    	worker(function_id, intensity, a, b, n);}
    
    
  	high_resolution_clock::time_point end_time = high_resolution_clock::now();
    //evaluating the duration 
    duration<double> time_period = duration_cast<duration<double>>(end_time - start_time);
    if (rank == 0) {
    //print the integration result onto the console screen
    std::cout << result << std::endl;
    //print the time period for computing the integral onto stderr
    std::cerr << time_period.count() << std::endl;
    }
    
    
  	MPI_Finalize();
  	return 0;
}
