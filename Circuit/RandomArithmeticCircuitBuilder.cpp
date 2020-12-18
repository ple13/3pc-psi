#include <ctime>
#include <random>
#include <chrono>
#include <cassert>
#include <iostream>
#include "ArithmeticGate.h"
#include "ArithmeticCircuit.h"
#include "InvalidOperationException.h"
#include "RandomArithmeticCircuitBuilder.h"
#include "../Utility/CryptoUtility.h"

namespace Circuit
{
	ArithmeticCircuit *RandomArithmeticCircuitBuilder::SampleArithmeticCircuit(int inputSize, int outputSize, int size)
	{
		assert(inputSize >= 2);
		assert(outputSize >= 1);
		assert(size >= inputSize + outputSize);
		
		unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::uniform_int_distribution<int> distribution(0,4);
		
		std::vector<IGate *> gates;

		for (int i = inputSize; i < size; i++)
		{
			int gateChoice = distribution(generator);

			int range = std::min(i, size - outputSize);
			
// 			std::cout << "gate choice: " << gateChoice << std::endl;
			
			if (gateChoice == 0)
			{
				int leftWire = sampleIndex(range);
				int rightWire = sampleFreshIndex(range,leftWire);
				IGate *tempVar = new AdditionGate(leftWire, rightWire, i);
				gates.push_back(tempVar);
			}
			else if (gateChoice == 1)
			{
				int leftWire = sampleIndex(range);
				int rightWire = sampleFreshIndex(range,leftWire);
				IGate *tempVar = new SubtractionGate(leftWire, rightWire, i);
				gates.push_back(tempVar);
			}
			else if (gateChoice == 2)
			{
				int leftWire = sampleIndex(range);
				int rightWire = sampleFreshIndex(range,leftWire);
				IGate *tempVar = new MultiplicationGate(leftWire, rightWire, i);
				gates.push_back(tempVar);
			}
			else if (gateChoice == 3)
			{
				int64_t constant = (2 << PRECISION_BIT_LENGTH);
				int inputWire = sampleIndex(range);
				IGate *tempVar = new ConstantAdditionGate(constant, inputWire, i);
				gates.push_back(tempVar);
			}
			else if (gateChoice == 4)
			{
				int64_t constant = (2 << PRECISION_BIT_LENGTH);
				int inputWire = sampleIndex(range);
				IGate *tempVar = new ConstantMultiplicationGate(constant, inputWire, i);
				gates.push_back(tempVar);
			}
			else if (gateChoice == 5)
			{
				int leftWire = sampleIndex(range);
				int rightWire = sampleFreshIndex(range,leftWire);
				IGate *tempVar = new DivisionGate(leftWire, rightWire, i);
				gates.push_back(tempVar);
			}
			else
			{
				throw InvalidOperationException();
			}
		}
		
		std::cout << "Before return an arithmetic circuit" << std::endl;
		
		ArithmeticCircuit *ret = new ArithmeticCircuit(size - inputSize, size, inputSize, outputSize, gates);
		
		std::cout << "Circuit params: " << ret->InputLength << " " << ret->OutputLength << " " << ret->NumGates << " " << ret->NumWires <<  " " << ret->Gates.size() << std::endl;
		
// 		std::cout << ret->ToString() << std::endl;
		
		std::cout << "Return an arithmetic circuit" << std::endl;
		
		return ret;
	}
	
	int RandomArithmeticCircuitBuilder::sampleIndex(int range) 
	{
		unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine gen(seed);
		std::uniform_int_distribution<int> dist(0, range - 1);
		return dist(gen);
	}
	
	int RandomArithmeticCircuitBuilder::sampleFreshIndex(int range, int previousIndex)
	{
		int result = 0;
		do
		{
			result = sampleIndex(range);
		} while (result == previousIndex);

		return result;
	}
}
