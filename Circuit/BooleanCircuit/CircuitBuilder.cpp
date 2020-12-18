#include <cassert>
#include <fstream>
#include <iostream>
#include "CircuitBuilder.h"
#include "BooleanGate.h"
#include "BooleanCircuit.h"
#include "stringhelper.h"

namespace Circuit
{

	BooleanCircuit *BooleanCircuitBuilder::CreateCircuit(const std::string &filePath)
	{
		{
			std::ifstream sr(filePath);
			std::string line;
			
			std::getline(sr, line);
			AddCircuitSize(line);
			
			std::getline(sr, line);
			AddInputOutputSizes(line);
			
			std::getline(sr, line); // empty line
			
			gates = std::vector<IGate *>(numGates);

			for (int i = 0; i < numGates; i++)
			{
				std::getline(sr, line);
				AddGate(line, i);
			}
		}

		return new BooleanCircuit(numGates, numWires, inputSize,outputSize,gates);
	}

	int BooleanCircuitBuilder::ConvertToInt(const std::string &str)
	{
		return std::stoi(str);
	}

	void BooleanCircuitBuilder::AddCircuitSize(const std::string &str)
	{
		auto sizes = StringHelper::split(str, ' ');
		assert(sizes.size() == 2);
		numGates = std::stoi(sizes[0]);
		numWires = std::stoi(sizes[1]);
	}

	void BooleanCircuitBuilder::AddInputOutputSizes(const std::string &str)
	{
		auto variables = StringHelper::split(str, ' ');
		assert(variables.size() == 2);

		std::function<int(int)> f = [variables] (int x)
		{
			return std::stoi(variables[x]);
		};

		inputSize = f(0);
		outputSize = f(1);

	}

	void BooleanCircuitBuilder::AddGate(const std::string &str, int index)
	{
		auto variables = StringHelper::split(str, ' ');

		std::function<int(int)> f = [variables] (int x)
		{
			return std::stoi(variables[x]);
		};

		int inputSize = f(0);

		if (inputSize == 1)
		{
			assert(variables.size() == 5);
			int inputWire = f(2);
			int outputWire = f(3);
			gates[index] = new NotGate(inputWire, outputWire);
		}
		else
		{
			assert(variables.size() == 6);
			int leftWire = f(2);
			int rightWire = f(3);
			int outputWire = f(4);
			
			std::string gate = variables[5];

			if (gate[0] == 'A')
			{
				gates[index] = new AndGate(leftWire, rightWire, outputWire);
			}
			else if (gate[0] == 'X')
			{
				gates[index] = new XorGate(leftWire, rightWire, outputWire);
			}
			else
			{
				throw UnreachableCodeException("Gate " + gate + " does not exist.");
			}
		}
	}

	BooleanCircuitBuilder::UnreachableCodeException::UnreachableCodeException(const std::string &str) : std::exception(/*str*/)
	{
	}
}
