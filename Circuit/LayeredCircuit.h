#ifndef LAYERED_CIRCUIT_H__
#define LAYERED_CIRCUIT_H__

#pragma once

#include <vector>
#include <type_traits>


namespace Circuit
{
	class ILayer
	{
	};

	/// <summary>
	/// Representation of a circuit divided into layers.
	/// </summary>
	template<typename T>
	class LayeredCircuit : public ILayer
	{
	private:
		std::vector<T *> Layers;
	public:
		const int Depth;
		const int NumWires;
		const int NumMulGates;
		const int InputLength;
		const int OutputLength;

		/// <summary>
		/// Representation of a circuit divided into layers.
		/// </summary>
		virtual ~LayeredCircuit()
		{
			for(int idx = 0; idx < Layers.size(); idx++)
			{
				delete Layers[idx];
			}
		}

		LayeredCircuit(std::vector<T *> layers, int depth, int numWires, int numMulGates, int inputLength, int outputLength) : Layers(std::move(layers)), Depth(depth), NumWires(numWires), NumMulGates(numMulGates), InputLength(inputLength), OutputLength(outputLength)
		{
		}

		/// <summary>
		/// Returns the ith layer of a circuit.
		/// </summary>
		T * operator [](int layerID)
		{
			return Layers[layerID];
		}
	};
}

#endif
