#ifndef RANDOM_ARITHMETIC_CIRCUIT_BUILDER_H__
#define RANDOM_ARITHMETIC_CIRCUIT_BUILDER_H__

#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <functional>
#include "ArithmeticCircuit.h"
#include "LayeredArithmeticCircuitBuilder.h"


namespace Circuit
{
	class RandomArithmeticCircuitBuilder
	{
	public:
		static void Sample(int inputSize, int outputSize, int circuitSize, ArithmeticCircuit* &circuit, LayeredArithmeticCircuit* &lc)
		{
			assert(inputSize >= 2);
			assert(outputSize >= 1);
			assert(circuitSize >= inputSize + outputSize);

			auto rcb = new RandomArithmeticCircuitBuilder();
			auto lcb = new LayeredArithmeticCircuitBuilder();
			
			std::cout << "Sample a Arithmetic Circuit" << std::endl;
			circuit = rcb->SampleArithmeticCircuit(inputSize, outputSize, circuitSize);
			
			
			std::cout << "Circuit params in .h: " << circuit->InputLength << " " << circuit->OutputLength << " " << circuit->NumGates << " " << circuit->NumWires << " " << circuit->Gates.size() << std::endl;
			std::cout << "Create layered Circuit" << std::endl;
			lc = lcb->CreateLayeredCircuit(circuit);
		}


	public:
		ArithmeticCircuit *SampleArithmeticCircuit(int inputSize, int outputSize, int size);
		int sampleIndex(int range);
		int sampleFreshIndex(int range, int previousIndex);
	};
}

#endif
