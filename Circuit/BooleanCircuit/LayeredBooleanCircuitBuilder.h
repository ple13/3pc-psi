#ifndef LAYERED_BOOOLEAN_CIRCUIT_BUILDER_H__
#define LAYERED_BOOOLEAN_CIRCUIT_BUILDER_H__

#pragma once

#include <unordered_map>
#include <vector>

#include "BooleanCircuit.h"
#include "LayeredBooleanCircuit.h"

namespace Circuit
{
	/// <summary>
	/// Creates a layered representation of a circuit from a <c>Circuit</c>.
	/// </summary>
	class LayeredBooleanCircuitBuilder
	{
	private:
		std::unordered_map<int, std::vector<AndGate*>> andGates;
		std::unordered_map<int, std::vector<XorGate*>> xorGates;
		std::unordered_map<int, std::vector<NotGate*>> notGates;
		std::unordered_map<int, int> mapWireToLayer;
		int circuitDepth = 0;
		int numAndGates = 0;
		static constexpr int INPUT_LAYER = -1;

		/// <summary>
		/// Creates a layered representation of a circuit from a <c>Circuit</c>.
		/// </summary>
	public:
		LayeredBooleanCircuitBuilder();

		/// <summary>
		/// Creation of layered circuit from circuit.
		/// </summary>
		LayeredBooleanCircuit *CreateLayeredCircuit(BooleanCircuit *circuit);

		/// <summary>
		/// Looks for all gates of a given type in a given layer. 
		/// </summary>
	private:
		template<typename T>
		std::vector<T> CreateGateList(std::unordered_map<int, std::vector<T>> &dict, int i);

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


		/// <summary>
		/// Prepares a gate to be added to a given layer.
		/// </summary>
		template<typename T>
		void AddToLayer(std::unordered_map<int, std::vector<T>> &layerDict, int layer, T t);
	};
}

#endif
