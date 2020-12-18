#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "Machine.hpp"
#include "Graph.hpp"
#include "Utility/Timer.h"

using namespace std;

int main(int argc, char** argv) 
{
	if (argc < 8)
		printf("Error: More arguments are needed. \n");
	int totalMachines = atoi(argv[1]);
	int machineId = atoi(argv[2]);
	int peerBasePort = atoi(argv[3]);
	int party = atoi(argv[4]);
	int partyBasePort = atoi(argv[5]);
	int itemLength = atoi(argv[6]);
	int test = atoi(argv[7]);
	float rate = atof(argv[8]);

	std::cout << totalMachines << " " << machineId << " " << party << std::endl;
	
	Machine * machine = new Machine(totalMachines, machineId, party);
	machine->connecting(peerBasePort, partyBasePort);
	
	cout << "Networking Done." << endl;
	
	usleep(5000);
	
	Graph * graph = new Graph(machineId, party, itemLength, machine, rate);
	
	Timer t;
	graph->GraphComputation(test);
	t.Tick("GraphComputation()");

	delete machine;
}