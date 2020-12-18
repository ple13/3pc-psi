#ifndef CIRCUIT_H__
#define CIRCUIT_H__

#pragma once

#include <vector>
#include <type_traits>
#include "Gate.h"

namespace Circuit
{
	/// <summary>
	/// Represents a circuits using the number of gates, the number of wires, and the gates.
	/// </summary>
	template<typename T>
	class ICircuit
	{
	public:
		const int InputLength;
		const int OutputLength;
		const int NumGates;
		const int NumWires;
		std::vector<T *> Gates;
		
		/// <summary>
		/// Returns the ith gate.
		/// </summary>
	public:
		/// <summary>
		/// Initializes a circuit.
		/// </summary>
		ICircuit(int numGates, int numWires, int inputLength, int outputLength, const std::vector<T *>& gates) : InputLength(inputLength), OutputLength(outputLength), NumGates(numGates), NumWires(numWires), Gates(gates)
		{
		}
		
		ICircuit(int numGates, int numWires, int inputLength, int outputLength, std::vector<T *>&& gates) : InputLength(inputLength), OutputLength(outputLength), NumGates(numGates), NumWires(numWires), Gates(gates)
		{
		}
		
		~ICircuit()
		{
			for(int idx = 0; idx < Gates.size(); idx++)
			{
				delete Gates[idx];
			}
		}
		
		T * operator [](int i)
		{
			return Gates[i];
		}
	};
}

#endif
