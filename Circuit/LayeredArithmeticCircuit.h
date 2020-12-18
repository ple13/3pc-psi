#ifndef LAYERED_ARITHMETIC_CIRCUIT_H__
#define LAYERED_ARITHMETIC_CIRCUIT_H__

#pragma once

#include <string>
#include <vector>

#include "LayeredCircuit.h"
#include "Circuit.h"
#include "ArithmeticGate.h"

namespace Circuit
{
	/// <summary>
	/// A class used to represnet a layer of a circuit.
	/// </summary>
	class ArithmeticLayer
	{
	public:
		std::vector<AdditionGate *> AddGates;
		std::vector<SubtractionGate *> SubGates;
		std::vector<MultiplicationGate *> MulGates;
		std::vector<DotProductGate *> DotGates;
		std::vector<DivisionGate *> DivGates;
		std::vector<NaryAdditionGate *> NAddGates;
		std::vector<NaryDotGate *> NDotGates;
		std::vector<ConstantAdditionGate *> CAddGates;
		std::vector<ConstantMultiplicationGate *> CMulGates;

		/// We may not need to delete data here
		virtual ~ArithmeticLayer() {}

		ArithmeticLayer(std::vector<AdditionGate *> AddGates, 
				std::vector<SubtractionGate *> SubGates,
				std::vector<MultiplicationGate *> MulGates, 
				std::vector<DotProductGate *> DotGates,
				std::vector<DivisionGate *> DivGates,
				std::vector<NaryAdditionGate *> NAddGates,
				std::vector<NaryDotGate *> NDotGates,
				std::vector<ConstantAdditionGate *> CAddGates, 
				std::vector<ConstantMultiplicationGate *> CMulGates);

		std::string ToString();
	};

	class LayeredArithmeticCircuit : public LayeredCircuit<ArithmeticLayer>
	{
	public:
		LayeredArithmeticCircuit(std::vector<ArithmeticLayer *> &layers, int depth, int numWires, int numMulGates, int inputLength, int outputLength);

		std::string ToString();
	};
}

#endif
