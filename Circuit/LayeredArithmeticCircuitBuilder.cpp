#include <cassert>
#include <iostream>
#include "Gate.h"
#include "ArithmeticGate.h"
#include "ArithmeticCircuit.h"
#include "LayeredArithmeticCircuit.h"
#include "LayeredArithmeticCircuitBuilder.h"


namespace Circuit
{

	LayeredArithmeticCircuitBuilder::LayeredArithmeticCircuitBuilder()
	{
		addGates  = std::unordered_map<int, std::vector<AdditionGate *> >();
		subGates  = std::unordered_map<int, std::vector<SubtractionGate *> >();
		mulGates  = std::unordered_map<int, std::vector<MultiplicationGate *> >();
		dotGates  = std::unordered_map<int, std::vector<DotProductGate *> >();
		divGates  = std::unordered_map<int, std::vector<DivisionGate *> >();
		nAddGates = std::unordered_map<int, std::vector<NaryAdditionGate *> >();
		nDotGates = std::unordered_map<int, std::vector<NaryDotGate *> >();
		cAddGates = std::unordered_map<int, std::vector<ConstantAdditionGate *> >();
		cMulGates = std::unordered_map<int, std::vector<ConstantMultiplicationGate *> >();
		
		mapWireToLayer = std::unordered_map<int, int>();
	}

	LayeredArithmeticCircuit *LayeredArithmeticCircuitBuilder::CreateLayeredCircuit(ArithmeticCircuit *circuit)
	{
		assert(circuit != nullptr);
		
		std::cout << "LayeredArithmeticCircuitBuilder::CreateLayeredCircuit" << std::endl;
		
		for (int i = 0; i < circuit->NumGates; i++)
		{
			IGate *igate = circuit->operator[](i);
			
			if (dynamic_cast<UnitaryGate *>(igate) != NULL)
			{
				PrepareUnitaryGate((UnitaryGate *)igate, i);
			}
			else if (dynamic_cast<BinaryGate *>(igate) != NULL)
			{
				PrepareBinaryGate((BinaryGate *)igate, i);
			}
			else if (dynamic_cast<NaryGate *>(igate) != NULL)
			{
				PrepareNaryGate((NaryGate *)igate, i);
			}
		}
		
		for (auto &kv : addGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : subGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : mulGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : dotGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : divGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : cAddGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : cMulGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		for (auto &kv : nAddGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}

		for (auto &kv : nDotGates)
		{
			circuitDepth = std::max(circuitDepth, kv.first);
		}
		
		circuitDepth++;
		
		std::vector<ArithmeticLayer *> layers(circuitDepth);

		for (int layer = 0; layer < circuitDepth; layer++)
		{
			auto AddGateArray = CreateGateList<AdditionGate>(addGates, layer);
			auto SubGateArray = CreateGateList<SubtractionGate>(subGates, layer);
			auto MulGateArray = CreateGateList<MultiplicationGate>(mulGates, layer);
			auto DotGateArray = CreateGateList<DotProductGate>(dotGates, layer);
			auto DivGateArray = CreateGateList<DivisionGate>(divGates, layer);
			auto NAddGateArray = CreateGateList<NaryAdditionGate>(nAddGates, layer);
			auto NDotGateArray = CreateGateList<NaryDotGate>(nDotGates, layer);
			auto CAddGateArray = CreateGateList<ConstantAdditionGate>(cAddGates, layer);
			auto CMulGateArray = CreateGateList<ConstantMultiplicationGate>(cMulGates, layer);
			
			layers[layer] = new ArithmeticLayer(AddGateArray, SubGateArray, MulGateArray, DotGateArray, DivGateArray, NAddGateArray, NDotGateArray, CAddGateArray, CMulGateArray);
		}
		
		LayeredArithmeticCircuit *ret = new LayeredArithmeticCircuit(layers, circuitDepth, circuit->NumWires, numMulGates, circuit->InputLength,circuit->OutputLength);
		return ret;
	}

	template<typename T>
	std::vector<T *> LayeredArithmeticCircuitBuilder::CreateGateList(std::unordered_map<int, std::vector<T *> > &dict, int i)
	{
		if (dict.find(i) == dict.end())
		{
			return std::vector<T *>();
		}
		else
		{
			return dict[i];
		}
	}

	int LayeredArithmeticCircuitBuilder::MaxDepth(int depth, int wire)
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

	void LayeredArithmeticCircuitBuilder::PrepareBinaryGate(BinaryGate *bg, int index)
	{
		assert(bg != nullptr);
		
		int layer = MaxDepth(INPUT_LAYER, bg->LeftWire);
		layer = MaxDepth(layer, bg->RightWire);
		layer++;

		mapWireToLayer[bg->OutputWire] = layer;

		if (dynamic_cast<AdditionGate *>(bg) != nullptr)
		{
			AddToLayer<AdditionGate>(addGates, layer, (AdditionGate *)bg);
		}
		else if (dynamic_cast<SubtractionGate *>(bg) != nullptr)
		{
			AddToLayer<SubtractionGate>(subGates, layer, (SubtractionGate *)bg);
		}
		else if (dynamic_cast<MultiplicationGate *>(bg) != nullptr)
		{	
			numMulGates++;
			AddToLayer<MultiplicationGate>(mulGates, layer, (MultiplicationGate *)bg);
		}
		else if (dynamic_cast<DotProductGate *>(bg) != nullptr)
		{	
			numMulGates++;
			AddToLayer<DotProductGate>(dotGates, layer, (DotProductGate *)bg);
		}
		else if (dynamic_cast<DivisionGate *>(bg) != nullptr)
		{	
			numMulGates++;
			AddToLayer<DivisionGate>(divGates, layer, (DivisionGate *)bg);
		}
	}

	void LayeredArithmeticCircuitBuilder::PrepareUnitaryGate(UnitaryGate *ug, int index)
	{
		assert(ug != nullptr);
		
		int layer = MaxDepth(INPUT_LAYER, ug->InputWire);
		layer++;

		mapWireToLayer[ug->OutputWire] = layer;

		if (dynamic_cast<ConstantAdditionGate *>(ug) != nullptr)
		{
			AddToLayer<ConstantAdditionGate>(cAddGates, layer, (ConstantAdditionGate *)ug);
		}
		else if (dynamic_cast<ConstantMultiplicationGate *>(ug) != nullptr)
		{
			AddToLayer<ConstantMultiplicationGate>(cMulGates, layer, (ConstantMultiplicationGate *)ug);
		}
	}
	
	void LayeredArithmeticCircuitBuilder::PrepareNaryGate(NaryGate *ng, int index)
	{
		assert(ng != nullptr);
		int layer = INPUT_LAYER;
		
		if (dynamic_cast<NaryAdditionGate *>(ng) != nullptr)
		{
			NaryAdditionGate *gate = (NaryAdditionGate *)ng;
			for(int idx = 0; idx < gate->InputWires.size(); idx++)
			{
				layer = MaxDepth(layer, gate->InputWires[idx]);
			}
			layer++;
			mapWireToLayer[gate->OutputWire] = layer;
			AddToLayer<NaryAdditionGate>(nAddGates, layer, gate);
		}
		else if (dynamic_cast<NaryDotGate *>(ng) != nullptr)
		{
			NaryDotGate *gate = (NaryDotGate *)ng;
			for(int idx = 0; idx < gate->InputWires.size(); idx++)
			{
				layer = MaxDepth(layer, gate->InputWires[idx]);
			}
			layer++;
			mapWireToLayer[gate->OutputWire] = layer;
			AddToLayer<NaryDotGate>(nDotGates, layer, gate);
		}
	}

	template<typename T>
	void LayeredArithmeticCircuitBuilder::AddToLayer(std::unordered_map<int, std::vector<T*> > &layerDict, int layer, T *t)
	{
		assert(nullptr != t);
		
		if (layerDict.find(layer) == layerDict.end())
		{
			layerDict.emplace(layer, std::vector<T *>());
		}

		layerDict[layer].push_back(t);
	}
}
