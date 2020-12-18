#ifndef RANDOM_CIRCUIT_BUILDER_H__
#define RANDOM_CIRCUIT_BUILDER_H__

#pragma once

#include <cassert>
#include <vector>
#include <functional>
#include <iostream>
#include "BooleanCircuit.h"
#include "LayeredBooleanCircuitBuilder.h"


namespace Circuit
{
	class RandomCircuitBuilder
	{
	public:
		static void Sample(int inputSize, int outputSize, int circuitSize, BooleanCircuit* circuit, LayeredBooleanCircuit* lc)
		{
			assert(inputSize >= 2);
			assert(outputSize >= 1);
			assert(circuitSize >= inputSize + outputSize);

			auto rcb = new RandomCircuitBuilder();
			auto lcb = new LayeredBooleanCircuitBuilder();
			
			std::cout << "Sample a Boolean Circuit" << std::endl;
			circuit = rcb->SampleBooleanCircuit(inputSize, outputSize, circuitSize);
			
			std::cout << "Create layered Circuit" << std::endl;
			lc = lcb->CreateLayeredCircuit(circuit);
		}


	public:
		BooleanCircuit *SampleBooleanCircuit(int inputSize, int outputSize, int size);
		int sampleIndex(int range);
		int sampleFreshIndex(int range, int previousIndex);
	};
}

#endif
