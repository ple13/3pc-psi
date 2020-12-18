#include <cassert>
#include <iostream>
#include "LayeredBooleanCircuitBuilder.h"
#include "BooleanGate.h"
#include "BooleanCircuit.h"
#include "LayeredBooleanCircuit.h"
#include "Gate.h"


namespace Circuit
{

	LayeredBooleanCircuitBuilder::LayeredBooleanCircuitBuilder()
	{
		andGates = std::unordered_map<int, std::vector<AndGate*>>();
		xorGates = std::unordered_map<int, std::vector<XorGate*>>();
		notGates = std::unordered_map<int, std::vector<NotGate*>>();
		mapWireToLayer = std::unordered_map<int, int>();
	}

	LayeredBooleanCircuit *LayeredBooleanCircuitBuilder::CreateLayeredCircuit(BooleanCircuit *circuit)
	{
// 		MyDebug::NonNull(circuit);
		std::cout << "LayeredBooleanCircuitBuilder::CreateLayeredCircuit" << std::endl;
		std::cout << "Num gates: " << circuit->NumGates << std::endl;
		
		for (int i = 0; i < circuit->NumGates; i++)
		{
// 			std::cout << "i: " << i << std::endl;
			IGate *igate = circuit->operator[](i);
// 			std::cout << igate << std::endl;
			
			if (dynamic_cast<UnitaryGate *>(igate) != NULL)
			{
// 				std::cout << "if igate" << std::endl;
				PrepareUnitaryGate((UnitaryGate *)igate, i);
			}
			else if (NULL != dynamic_cast<BinaryGate *>(igate))
			{
// 				std::cout << "else igate" << std::endl;
				PrepareBinaryGate((BinaryGate *)igate, i);
			}
		}
		
// 		std::cout << "Find max depth " << std::endl;
		
		for (auto &kv : andGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : xorGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : notGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}

		circuitDepth++;
		
		std::cout << "circuit depth: " << circuitDepth << std::endl;
		
		std::vector<BooleanLayer*> layers(circuitDepth);

		for (int layer = 0; layer < circuitDepth; layer++)
		{
			auto AndGateArray = CreateGateList<AndGate*>(andGates, layer);
			auto XorGateArray = CreateGateList<XorGate*>(xorGates, layer);
			auto NotGateArray = CreateGateList<NotGate*>(notGates, layer);
			
			layers[layer] = new BooleanLayer(AndGateArray, XorGateArray, NotGateArray);
		}

		return new LayeredBooleanCircuit(layers, circuitDepth, circuit->NumWires, numAndGates, circuit->InputLength,circuit->OutputLength);
	}

template<typename T>
	std::vector<T> LayeredBooleanCircuitBuilder::CreateGateList(std::unordered_map<int, std::vector<T>> &dict, int i)
	{
// 		MyDebug::NonNull(dict);

		if (dict.find(i) == dict.end())
		{
			return std::vector<T>();
		}
		else
		{
			return dict[i];
		}
	}

	int LayeredBooleanCircuitBuilder::MaxDepth(int depth, int wire)
	{
		if (mapWireToLayer.find(wire) == mapWireToLayer.end())
		{
			return depth;
		}
		else
		{
			return std::max(depth, mapWireToLayer[wire]);
		}
	}

	void LayeredBooleanCircuitBuilder::PrepareBinaryGate(BinaryGate *bg, int index)
	{
// 		std::cout << "LayeredBooleanCircuitBuilder::PrepareBinaryGate" << std::endl;
// 		MyDebug::NonNull(bg);

		int layer = MaxDepth(INPUT_LAYER, bg->LeftWire);
		layer = MaxDepth(layer, bg->RightWire);
		layer++;

		mapWireToLayer[bg->OutputWire] = layer;

		if (dynamic_cast<AndGate*>(bg) != nullptr)
		{
			numAndGates++;
			AddToLayer<AndGate*>(andGates, layer, (AndGate *)bg);
		}
		else if (dynamic_cast<XorGate*>(bg) != nullptr)
		{
			AddToLayer<XorGate*>(xorGates, layer, (XorGate *)bg);
		}
	}

	void LayeredBooleanCircuitBuilder::PrepareUnitaryGate(UnitaryGate *ug, int index)
	{
// 		std::cout << "LayeredBooleanCircuitBuilder::PrepareUnitaryGate" << std::endl;
		assert(dynamic_cast<NotGate*>(ug) != nullptr);
// 		MyDebug::NonNull(ug);

		int layer = MaxDepth(INPUT_LAYER, ug->InputWire);
		layer++;

		mapWireToLayer[ug->OutputWire] = layer;

		if (dynamic_cast<NotGate*>(ug) != nullptr)
		{
			AddToLayer<NotGate*>(notGates, layer, (NotGate *)ug);
		}
	}

template<typename T>
	void LayeredBooleanCircuitBuilder::AddToLayer(std::unordered_map<int, std::vector<T>> &layerDict, int layer, T t)
	{
// 		MyDebug::NonNull(layerDict, t);

		if (layerDict.find(layer) == layerDict.end())
		{
			layerDict.emplace(layer, std::vector<T>());
		}

		layerDict[layer].push_back(t);
	}
}
