#include <string>
#include "ArithmeticGate.h"
#include "ArithmeticCircuit.h"


namespace Circuit
{

	ArithmeticCircuit::ArithmeticCircuit(int numGates, int numWires, int inputLength, int outputLength, const std::vector<IGate *>& gates) : ICircuit<IGate>(numGates, numWires, inputLength, outputLength, gates)
	{
	}
	
	ArithmeticCircuit::ArithmeticCircuit(int numGates, int numWires, int inputLength, int outputLength, std::vector<IGate *>&& gates) : ICircuit<IGate>(numGates, numWires, inputLength, outputLength, gates)
	{
	}

	std::string ArithmeticCircuit::ToString()
	{
		std::string builder;

		for (int i = 0; i < NumGates; i++)
		{
			IGate *gate = this->operator[](i);

			if (nullptr != dynamic_cast<AdditionGate *>(gate))
			{
				builder = builder + "ADD" + "|" + std::to_string(((AdditionGate *)gate)->LeftWire) 
							  + "|" + std::to_string(((AdditionGate *)gate)->RightWire) 
							  + "|" + std::to_string(((AdditionGate *)gate)->OutputWire) 
							  + "\n";
			}
			else if (nullptr != dynamic_cast<SubtractionGate *>(gate))
			{
				builder = builder + "SUB" + "|" + std::to_string(((SubtractionGate *)gate)->LeftWire) 
							  + "|" + std::to_string(((SubtractionGate *)gate)->RightWire) 
							  + "|" + std::to_string(((SubtractionGate *)gate)->OutputWire) 
							  + "\n";
			}
			else if (nullptr != dynamic_cast<MultiplicationGate *>(gate))
			{
				builder = builder + "MUL" + "|" + std::to_string(((MultiplicationGate *)gate)->LeftWire) 
							  + "|" + std::to_string(((MultiplicationGate *)gate)->RightWire) 
							  + "|" + std::to_string(((MultiplicationGate *)gate)->OutputWire) 
							  + "\n";
			}
			else if (nullptr != dynamic_cast<DotProductGate *>(gate))
			{
				builder = builder + "DOT" + "|" + std::to_string(((DotProductGate *)gate)->LeftWire) 
							  + "|" + std::to_string(((DotProductGate *)gate)->RightWire) 
							  + "|" + std::to_string(((DotProductGate *)gate)->OutputWire) 
							  + "\n";
			}
			else if (nullptr != dynamic_cast<DivisionGate *>(gate))
			{
				builder = builder + "DIV" + "|" + std::to_string(((DivisionGate *)gate)->LeftWire) 
							  + "|" + std::to_string(((DivisionGate *)gate)->RightWire) 
							  + "|" + std::to_string(((DivisionGate *)gate)->OutputWire) 
							  + "\n";
			}
			else if (nullptr != dynamic_cast<ConstantAdditionGate *>(gate))
			{
				builder = builder + "CADD" + "|(" + std::to_string(((ConstantAdditionGate *)gate)->constant) 
							   + ")|" + std::to_string(((ConstantAdditionGate *)gate)->InputWire) 
							   + "|" + std::to_string(((ConstantAdditionGate *)gate)->OutputWire) 
							   + "\n";
			}
			else if (nullptr != dynamic_cast<ConstantMultiplicationGate *>(gate))
			{
				builder = builder + "CMUL" + "|(" + std::to_string(((ConstantMultiplicationGate *)gate)->constant) 
							   + ")|" + std::to_string(((ConstantMultiplicationGate *)gate)->InputWire) 
							   + "|"  + std::to_string(((ConstantMultiplicationGate *)gate)->OutputWire) 
							   + "\n";
			}
			else if (nullptr != dynamic_cast<NaryAdditionGate *>(gate))
			{
				NaryAdditionGate * g = (NaryAdditionGate *)gate;
				
				builder = builder + "NADD" + "|";
				for(int idx = 0; idx < g->InputWires.size(); idx++)
				{
					builder = builder + std::to_string(g->InputWires[idx]) + "|";
				}
				builder = builder + std::to_string(g->OutputWire) + "\n";
			}
			else if (nullptr != dynamic_cast<NaryDotGate *>(gate))
			{
				NaryDotGate * g = (NaryDotGate *)gate;
				
				builder = builder + "NDOT" + "|";
				for(int idx = 0; idx < g->InputWires.size(); idx++)
				{
					builder = builder + std::to_string(g->InputWires[idx]) + "|";
				}
				builder = builder + std::to_string(g->OutputWire) + "\n";
			}
		}

		return builder;
	}
}
