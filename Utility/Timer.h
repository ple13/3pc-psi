#ifndef _TIMER_H__
#define _TIMER_H__

#include <sstream>
#include <iostream>
#include <sys/time.h>
#include <chrono>
#include <ctime>

class Timer 
{
public:
	struct timeval diff, startTV, endTV;
	
// 	std::chrono::time_point<std::chrono::system_clock> start, end;
	
	Timer()
	{
		gettimeofday(&startTV, NULL); 
// 		start = std::chrono::system_clock::now();
	}
	
	void Tick(std::string str)
	{
		gettimeofday(&endTV, NULL); 
		timersub(&endTV, &startTV, &diff);
		
		printf("**%s time: %f seconds\n", str.c_str(), diff.tv_sec + (float)diff.tv_usec/1000000.0);
		
		startTV = endTV;
		
// 		end = std::chrono::system_clock::now();
// 		std::chrono::duration<double> elapsed_seconds = end-start;
// 		std::time_t end_time = std::chrono::system_clock::to_time_t(end);
// 
// 		std::cout << "finished computation at " << std::ctime(&end_time) << "elapsed time: " << elapsed_seconds.count() << "s\n";
// 		
// 		start = end;
	}
};

#endif

//     auto start = std::chrono::system_clock::now();
//     // Some computation here
//     auto end = std::chrono::system_clock::now();
// 
//     std::chrono::duration<double> elapsed_seconds = end-start;
//     std::time_t end_time = std::chrono::system_clock::to_time_t(end);
// 
//     std::cout << "finished computation at " << std::ctime(&end_time)
//               << "elapsed time: " << elapsed_seconds.count() << "s\n";