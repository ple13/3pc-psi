#include <cassert>
#include <fstream>
#include <iostream>
#include <functional>
#include "stringhelper.h"
#include "ArithmeticGate.h"
#include "ArithmeticCircuit.h"
#include "ArithmeticCircuitBuilder.h"

namespace Circuit
{
	ArithmeticCircuit *ArithmeticCircuitBuilder::CreateCircuit(const std::string &filePath)
	{
		{
			std::ifstream sr(filePath);
			std::string line;
			
			std::getline(sr, line);
			
			std::getline(sr, line);
			AddCircuitSize(line);
			
			std::getline(sr, line);
			
			std::getline(sr, line);
			AddInputOutputSizes(line);
			
			std::getline(sr, line); // empty line
			
			gates = std::vector<IGate *>(numGates);

			int i = 0;
			while(std::getline(sr, line))
			{
				if(line[0] != '#')
				{
					AddGate(line, i);
					i++;
				}
			}
		}

		return new ArithmeticCircuit(numGates, numWires, inputSize, outputSize, gates);
	}

	int ArithmeticCircuitBuilder::ConvertToInt(const std::string &str)
	{
		return std::stoi(str);
	}

	void ArithmeticCircuitBuilder::AddCircuitSize(const std::string &str)
	{
		std::cout << "AddCircuitSize: " << str << std::endl;
		auto sizes = StringHelper::split(str, ' ');
		assert(sizes.size() == 2);
		numGates = std::stoi(sizes[0]);
		numWires = std::stoi(sizes[1]);
	}

	void ArithmeticCircuitBuilder::AddInputOutputSizes(const std::string &str)
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

	void ArithmeticCircuitBuilder::AddGate(const std::string &str, int index)
	{
// 		std::cout << "AddGate: " << str << std::endl;
		auto variables = StringHelper::split(str, ' ');

		std::function<int(int)> f = [variables] (int x)
		{
			return std::stoi(variables[x]);
		};

		int inputSize = f(0);
		
// 		std::cout << "inputSize: " << inputSize << std::endl;
		
		if(2 == inputSize)
		{
			assert(variables.size() == 6);
			
			int inputWireOrConstant = f(2);
			int inputWire = f(3);
			int outputWire = f(4);
			
			std::string gate = variables[5];
// 			std::cout << "gate: " << gate << std::endl;
			
			/// Gate's names: ADD, MUL, CADD, CMUL
			if (gate[0] == 'A')
			{
				gates[index] = new AdditionGate(inputWireOrConstant, inputWire, outputWire);
			}
			else if (gate[0] == 'S')
			{
				gates[index] = new SubtractionGate(inputWireOrConstant, inputWire, outputWire);
			}
			else if (gate[0] == 'M')
			{
				gates[index] = new MultiplicationGate(inputWireOrConstant, inputWire, outputWire);
			}
			else if (gate[0] == 'N')
			{
				if(gate[1] == 'A')
				{
					std::vector<int> inputWires;
					for(int idx = 0; idx < inputSize; idx++)
					{
						inputWires.push_back(f(2+idx));
					}
					int outputWire = f(2 + inputSize);
					gates[index] = new NaryAdditionGate(inputWires, outputWire);
				}
				else if (gate[1] == 'D')
				{
					std::vector<int> inputWires;
					for(int idx = 0; idx < inputSize; idx++)
					{
						inputWires.push_back(f(2+idx));
					}
					int outputWire = f(2 + inputSize);
					gates[index] = new NaryDotGate(inputWires, outputWire);
				}
			}
			else if (gate[0] == 'D' && gate[1] == 'O')
			{
				gates[index] = new DotProductGate(inputWireOrConstant, inputWire, outputWire);
			}
			else if (gate[0] == 'D' && gate[1] == 'I')
			{
// 				std::cout << "Div gate" << std::endl;
				gates[index] = new DivisionGate(inputWireOrConstant, inputWire, outputWire);
			}
			else if (gate[1] == 'A')
			{
				assert(gate[0] == 'C');
				gates[index] = new ConstantAdditionGate(inputWireOrConstant, inputWire, outputWire);
			}
			else if (gate[1] == 'M')
			{
				assert(gate[0] == 'C');
				float inputWireOrConstant = std::stof(variables[2]);
				gates[index] = new ConstantMultiplicationGate(inputWireOrConstant, inputWire, outputWire);
			}
			else
			{
				std::cout << "Here : gate[0] = " << gate[0] << "\tgate[1] = " << gate[1] << std::endl;
				throw UnreachableCodeException("Gate " + gate + " does not exist.");
			}
		}
		else if(inputSize != 2)
		{
			std::string gate = variables[inputSize + 3];
// 			std::cout << gate << std::endl;
			if(gate[1] == 'A')
			{
				std::vector<int> inputWires;
				for(int idx = 0; idx < inputSize; idx++)
				{
					inputWires.push_back(f(2+idx));
				}
				int outputWire = f(2 + inputSize);
				gates[index] = new NaryAdditionGate(inputWires, outputWire);
			}
			else if (gate[1] == 'D')
			{
				std::vector<int> inputWires;
				for(int idx = 0; idx < inputSize; idx++)
				{
					inputWires.push_back(f(2+idx));
				}
				int outputWire = f(2 + inputSize);
				gates[index] = new NaryDotGate(inputWires, outputWire);
			}
		}
	}

	ArithmeticCircuitBuilder::UnreachableCodeException::UnreachableCodeException(const std::string &str) : std::exception(/*str*/)
	{
	}
}
