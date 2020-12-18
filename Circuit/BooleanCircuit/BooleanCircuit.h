#ifndef BOOLEAN_CIRCUIT_H__
#define BOOLEAN_CIRCUIT_H__

#pragma once
#include "Circuit.h"
#include <string>
#include <vector>


namespace Circuit
{
	class BooleanCircuit : public Circuit<IGate>
	{
	public:
		BooleanCircuit(int numGates, int numWires, int inputLength, int outputLength, std::vector<IGate *> &gates);

		std::string ToString();
	};
}

#endif
