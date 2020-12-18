#ifndef _ARITHMETIC_CIRCUIT_H__
#define _ARITHMETIC_CIRCUIT_H__

#include <string>
#include "Circuit.h"

namespace Circuit
{
	class ArithmeticCircuit : public ICircuit<IGate>
	{
	public:
		ArithmeticCircuit(int numGates, int numWires, int inputLength, int outputLength, const std::vector<IGate *>& gates);
		ArithmeticCircuit(int numGates, int numWires, int inputLength, int outputLength, std::vector<IGate *>&& gates);

		std::string ToString();
	};
}

#endif