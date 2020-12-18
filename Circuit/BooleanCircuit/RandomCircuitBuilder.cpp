#include "RandomCircuitBuilder.h"
#include "BooleanCircuit.h"
#include "BooleanGate.h"
#include "InvalidOperationException.h"
#include <random>
#include <cassert>
#include <iostream>
#include <ctime>
#include <chrono>

namespace Circuit
{
	BooleanCircuit *RandomCircuitBuilder::SampleBooleanCircuit(int inputSize, int outputSize, int size)
	{
		assert(inputSize >= 2);
		assert(outputSize >= 1);
		assert(size >= inputSize + outputSize);
		
		unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::uniform_int_distribution<int> distribution(0,2);
		
		std::vector<IGate *> gates;

		for (int i = inputSize; i < size; i++)
		{
// 			std::cout << i << std::endl;
			int gateChoice = distribution(generator);

			int range = std::min(i, size - outputSize);
// 			std::cout << "range = " << range << std::endl;
			
			if (gateChoice == 0)
			{
				NotGate *tempVar = new NotGate(sampleIndex(range), i);
				gates.push_back(tempVar);
			}
			else if (gateChoice == 1)
			{
				int leftWire = sampleIndex(range);
				int rightWire = sampleFreshIndex(range,leftWire);
				XorGate *tempVar2 = new XorGate(leftWire, rightWire, i);
				gates.push_back(tempVar2);
			}
			else if (gateChoice == 2)
			{
				int leftWire = sampleIndex(range);
				int rightWire = sampleFreshIndex(range, leftWire);
				AndGate *tempVar3 = new AndGate(leftWire, rightWire, i);
				gates.push_back(tempVar3);

			}
			else
			{
				throw InvalidOperationException();
			}
		}

		return new BooleanCircuit(size - inputSize, size, inputSize, outputSize, gates);
	}
	
	int RandomCircuitBuilder::sampleIndex(int range) 
	{
		unsigned seed =  std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine gen(seed);
		std::uniform_int_distribution<int> dist(0, range - 1);
		return dist(gen);
	}
	
	int RandomCircuitBuilder::sampleFreshIndex(int range, int previousIndex)
	{
		int result = 0;
		do
		{
			result = sampleIndex(range);
// 			std::cout << previousIndex << " " << result << " " << range << std::endl;
		} while (result == previousIndex);

		return result;
	}
}
