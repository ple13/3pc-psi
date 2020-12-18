#include <algorithm>
#include "LayeredBooleanCircuit.h"
#include "BooleanGate.h"

namespace Circuit
{

	BooleanLayer::BooleanLayer(std::vector<AndGate *> &AndGates, std::vector<XorGate *> &XorGates, std::vector<NotGate *> &NotGates) : AndGates(AndGates), XorGates(XorGates), NotGates(NotGates)
	{
// 		MyDebug::NonNull(AndGates, XorGates, NotGates);
	}

	std::string BooleanLayer::ToString()
	{
		std::string builder;
		
		for (auto& ag : AndGates)
		{
			builder = builder + "And" + "|" + ag->ToString() + "\n";
		}

		for  (auto &xg : XorGates)
		{
			builder = builder + "XOR" + "|" + xg->ToString() + "\n";
		}

		for  (auto &ng : NotGates)
		{
			builder = builder + "NOT" + "|" + ng->ToString() + "\n";
		}

		return builder;
	}

	LayeredBooleanCircuit::LayeredBooleanCircuit(std::vector<BooleanLayer*> &layers, int depth, int numWires, int numAndGates, int inputLength, int outputLength) : LayeredCircuit<BooleanLayer *>(layers, depth, numWires, numAndGates, inputLength, outputLength)
	{
	}

	std::string LayeredBooleanCircuit::ToString()
	{
		std::string builder;
		
		builder += "\n";//		f("");
		builder += "---------------------------------------\n"; //f("---------------------------------------");

		for (int i = 0; i < Depth; i++)
		{
			BooleanLayer *layer = this->operator[](i);
			builder = builder + "Layer " + std::to_string(i) + "\n"; //f("Layer " + std::to_wstring(i));
			builder = builder + layer->ToString() + "\n"; //builder->append(layer->ToString());
			builder += "---------------------------------------\n"; //f("---------------------------------------");
		}

		return builder;
	}
}
