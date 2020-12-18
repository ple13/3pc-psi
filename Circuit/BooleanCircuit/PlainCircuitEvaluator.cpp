#include "PlainCircuitEvaluator.h"
#include "BooleanCircuit.h"
#include "Gate.h"
#include "BooleanGate.h"
#include <boost/dynamic_bitset.hpp>
#include <cassert>

namespace Circuit
{

	boost::dynamic_bitset<> *PlainBooleanCircuitEvaluator::EvaluateBooleanCircuit(BooleanCircuit *circuit, boost::dynamic_bitset<> *input)
	{
		assert(circuit->InputLength == input->size());

		boost::dynamic_bitset<> *wireValues = new boost::dynamic_bitset<>(*input);

		for (int i = 0;i < circuit->NumGates;i++)
		{
			IGate *g = circuit->operator[](i);
			if (dynamic_cast<BinaryGate*>(g) != nullptr)
			{
				bool leftValue = wireValues->operator[](((BinaryGate *)g)->LeftWire);
				bool rightValue = wireValues->operator[](((BinaryGate *)g)->RightWire);
				int outputWire = ((BinaryGate *)g)->OutputWire;

				if (dynamic_cast<AndGate*>(g) != nullptr)
				{
					wireValues->operator[](outputWire) = leftValue & rightValue;
// 					wireValues[outputWire] = leftValue & rightValue;
				}
				else if (dynamic_cast<XorGate*>(g) != nullptr)
				{
					wireValues->operator[](outputWire) = leftValue ^ rightValue;
// 					wireValues[outputWire] = leftValue ^ rightValue;
				}
				else
				{
					throw InvalidOperationException();
				}
			}
			else if (dynamic_cast<UnitaryGate*>(g) != nullptr)
			{
				wireValues->operator[](((UnitaryGate *)g)->OutputWire) = !wireValues->operator[](((UnitaryGate *)g)->InputWire);
			}
		}

		return wireValues;
	}
}
