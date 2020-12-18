#include <string>
#include "BooleanGate.h"
#include "BooleanCircuit.h"


namespace Circuit
{

	BooleanCircuit::BooleanCircuit(int numGates, int numWires, int inputLength, int outputLength, std::vector<IGate *> &gates) : Circuit<IGate *>(numGates, numWires, inputLength, outputLength, gates)
	{
	}

	std::string BooleanCircuit::ToString()
	{
		std::string builder;

		for (int i = 0; i < NumGates; i++)
		{
			IGate *gate = this->operator[](i);

			if (nullptr != dynamic_cast<AndGate *>(dynamic_cast<BinaryGate *>(gate)))
			{
				builder = builder + "AND" + "|" + std::to_string(((AndGate *)gate)->LeftWire) + "|" + std::to_string(((AndGate *)gate)->RightWire) + "|" + std::to_string(((AndGate *)gate)->OutputWire) + "\n";
			}
			else if (nullptr != dynamic_cast<XorGate *>(gate))
			{
				builder = builder + "XOR" + "|" + std::to_string(((XorGate *)gate)->LeftWire) + "|" + std::to_string(((XorGate *)gate)->RightWire) + "|" + std::to_string(((XorGate *)gate)->OutputWire) + "\n";
			}
			else if (nullptr != dynamic_cast<NotGate *>(gate))
			{
				builder = builder + "NOT" + "|" + std::to_string(((NotGate *)gate)->InputWire) + "|" + std::to_string(((NotGate *)gate)->OutputWire) + "\n";
			}
		}

		return builder;
	}
}
