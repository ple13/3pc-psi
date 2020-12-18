#include <string>
#include <vector>
#include <iostream>
#include <functional>
#include "BooleanCircuit.h"
#include "CircuitBuilder.h"
#include "LayeredBooleanCircuit.h"
#include "RandomCircuitBuilder.h"
#include "LayeredBooleanCircuitBuilder.h"


using namespace Circuit;

int main() {
	  BooleanCircuit *circuit;
	  LayeredBooleanCircuit *lc;
	  
	  int inputSize;
	  int outputSize;
	  int circuitSize;

	  inputSize = 256;
	  outputSize = 128;
	  circuitSize = 4000;
	  
	  std::cout << "Generate a random circuit" << std::endl;
	  RandomCircuitBuilder::Sample(inputSize, outputSize, circuitSize, circuit, lc);
	  std::cout << "Done and terminate&" << std::endl;
}