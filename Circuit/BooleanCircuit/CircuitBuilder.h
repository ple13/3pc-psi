#ifndef CIRCUIT_BUILDER_H__
#define CIRCUIT_BUILDER_H__

#pragma once

#include <string>
#include <vector>
#include <functional>
#include <exception>
#include "stringhelper.h"
#include "Gate.h"
#include "BooleanCircuit.h"


namespace Circuit
{
	enum class BinaryGates
	{
		AND,
		XOR,
		NOT
	};

	/// <summary>
	/// Class to construct a circuit from a file.
	/// </summary>
	class BooleanCircuitBuilder
	{
	private:
		int numGates = 0;
		int numWires = 0;
		std::vector<IGate *> gates;
		int inputSize = 0;
		int outputSize = 0;

		/// <summary>
		/// Creates a circuit from a file.
		/// </summary>
		/// <param name="filePath">  Filepath for file containing circuit. </param>
		/// <returns> </returns>
	public:
		BooleanCircuit *CreateCircuit(const std::string &filePath);

		/// <summary>
		/// Converts String to int.
		/// </summary>
	private:
		int ConvertToInt(const std::string &str);

		/// <summary>
		/// Adds number of gate and number of wires to circuit.
		/// </summary>
		void AddCircuitSize(const std::string &str);

		/// <summary>
		/// 
		/// </summary>
		void AddInputOutputSizes(const std::string &str);

		/// <summary>
		/// Adds a gate to circuit description.
		/// </summary>
		void AddGate(const std::string &str, int index);

	private:
		class UnreachableCodeException : public std::exception
		{
		public:
			UnreachableCodeException(const std::string &str);
		};
	};
}

#endif
