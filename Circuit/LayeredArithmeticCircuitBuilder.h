#ifndef LAYERED_ARITHMETIC_CIRCUIT_BUILDER_H__
#define LAYERED_ARITHMETIC_CIRCUIT_BUILDER_H__

#pragma once

#include <unordered_map>
#include <vector>

#include "ArithmeticCircuit.h"
#include "LayeredArithmeticCircuit.h"

namespace Circuit
{
	/// <summary>
	/// Creates a layered representation of a circuit from a <c>Circuit</c>.
	/// </summary>
	class LayeredArithmeticCircuitBuilder
	{
	private:
		std::unordered_map<int, std::vector<AdditionGate *> > addGates;
		std::unordered_map<int, std::vector<SubtractionGate *> > subGates;
		std::unordered_map<int, std::vector<MultiplicationGate *> > mulGates;
		std::unordered_map<int, std::vector<DotProductGate *> > dotGates;
		std::unordered_map<int, std::vector<DivisionGate *> > divGates;
		std::unordered_map<int, std::vector<NaryAdditionGate *> > nAddGates;
		std::unordered_map<int, std::vector<NaryDotGate *> > nDotGates;
		std::unordered_map<int, std::vector<ConstantAdditionGate *> > cAddGates;
		std::unordered_map<int, std::vector<ConstantMultiplicationGate *> > cMulGates;
		
		std::unordered_map<int, int> mapWireToLayer;
		
		int circuitDepth = 0;
		int numMulGates = 0;
		
		static constexpr int INPUT_LAYER = -1;

		/// <summary>
		/// Creates a layered representation of a circuit from a <c>Circuit</c>.
		/// </summary>
	public:
		LayeredArithmeticCircuitBuilder();

		/// <summary>
		/// Creation of layered circuit from circuit.
		/// </summary>
		LayeredArithmeticCircuit *CreateLayeredCircuit(ArithmeticCircuit *circuit);

		/// <summary>
		/// Looks for all gates of a given type in a given layer. 
		/// </summary>
	private:
		template<typename T>
		std::vector<T *> CreateGateList(std::unordered_map<int, std::vector<T*> > &dict, int i);

		/// <summary>
		/// Max between depth or Layer[wire].  
		/// </summary>
		int MaxDepth(int depth, int wire);

		/// <summary>
		/// Prepares a binary gate to be added to the layered circuit.
		/// </summary>
		void PrepareBinaryGate(BinaryGate *bg, int index);

		/// <summary>
		/// Prepares a unitary gate to be added to the layered circuit.
		/// </summary>
		void PrepareUnitaryGate(UnitaryGate *ug, int index);

		void PrepareNaryGate(NaryGate *g, int index);
		
		/// <summary>
		/// Prepares a gate to be added to a given layer.
		/// </summary>
		
		template<typename T>
		void AddToLayer(std::unordered_map<int, std::vector<T*> > &layerDict, int layer, T *t);
	};
}

#endif
