#ifndef LAYERED_BOOLEAN_CIRCUIT_H__
#define LAYERED_BOOLEAN_CIRCUIT_H__

#pragma once

#include "LayeredCircuit.h"
#include "Circuit.h"
#include <string>
#include <vector>
#include <functional>
#include "Circuit.h"
#include "BooleanGate.h"

namespace Circuit
{
	/// <summary>
	/// A class used to represnet a layer of a circuit.
	/// </summary>
	class BooleanLayer : public ILayer
	{
	public:
		std::vector<AndGate *> /* *const*/ AndGates;
		std::vector<XorGate *> /* *const*/ XorGates;
		std::vector<NotGate *> /* *const*/ NotGates;

		virtual ~BooleanLayer()
		{
			for(int idx = 0; idx < AndGates.size(); idx++) {
				delete AndGates[idx];
			}
			
			for(int idx = 0; idx < XorGates.size(); idx++) {
				delete XorGates[idx];
			}
			
			for(int idx = 0; idx < NotGates.size(); idx++) {
				delete NotGates[idx];
			}
		}

		BooleanLayer(std::vector<AndGate*> &AndGates, std::vector<XorGate*> &XorGates, std::vector<NotGate*> &NotGates);

		std::string ToString();
	};

	class LayeredBooleanCircuit : public LayeredCircuit<BooleanLayer *>
	{
	public:
		LayeredBooleanCircuit(std::vector<BooleanLayer *> &layers, int depth, int numWires, int numAndGates, int inputLength, int outputLength);

		std::string ToString();
	};
}

#endif
