#ifndef GRADIENT_DESCEND_CIRCUIT_BUILDER_H__
#define GRADIENT_DESCEND_CIRCUIT_BUILDER_H__

#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include "Gate.h"
#include "stringhelper.h"
#include "ArithmeticGate.h"
#include "ArithmeticCircuit.h"
#include "ArithmeticCircuitBuilder.h"
#include "InvalidOperationException.h"
#include "RandomArithmeticCircuitBuilder.h"
#include "../Utility/CryptoUtility.h"

#define DIM 10

using namespace Circuit;

class GradientDescendCircuitBuilder
	{
	public:
		static void BuildSingleUser(int numItems,float gamma, float muy, const std::string& outfile)
		{
			std::stringstream ss;
			ss << "#gates #wires\n";
			ss << std::to_string(5 + 4*numItems) << " " << std::to_string(6 + 7*numItems) << "\n";
			
			ss << "InputSize OutputSize\n";
			ss << std::to_string(1 + 3*numItems) << " " << std::to_string(1) << "\n";
			ss << "\n";
			
			int count = 1 + 3*numItems;
			
			for(int idx = 0; idx < numItems; idx++)
			{
				ss << "#####---item " << idx << "---#####\n";
				ss << "2 1 0 " << (1 + 3*idx) << " " << count << " DOT\n"; count++;				// <u,v>
				ss << "2 1 " << (2 + 3*idx) << " " << (count - 1) << " " << count << " SUB\n"; count++;		// r - <u,v>
				ss << "2 1 " << (count - 1) << " " << (1 + 3*idx) << " " << count << " MUL\n"; count++;		// v(r - <u,v>)
				ss << "2 1 " << (3 + 3*idx) << " " << (count - 1) << " " << count << " MUL\n"; count++;		// isReal*v(r - <u,v>)
			}
			
			// Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)
			ss << numItems << " 1 ";
			for(int idx = 1; idx <= numItems; idx++)
			{
				 ss << (3*numItems + 4*idx) << " ";
			}
			ss << count << " NADD\n"; 
			count++;
			
			int gamma_int = (int)(gamma*1048576.0);
			int muy_int = (int)(muy*1048576.0);
			
			ss << "2 1 " << gamma << " " << (count - 1) << " " << count << " CMUL\n"; count++;			// gamma*[Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)]
			ss << "2 1 " << muy << " " << 0 << " " << count << " CMUL\n"; count++;					// muy*u_i
			ss << "2 1 " << (count - 2) << " " << (count - 1) << " " << count << " ADD\n"; count++;			// gamma*[Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)] + muy*u_i
			ss << "2 1 " << 0 << " " << (count - 1) << " " << count << " ADD\n"; count++;				// u_i + gamma*[Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)] + muy*u_i
			
			std::ofstream myfile;
			myfile.open(outfile);
			myfile << ss.str();
			myfile.close();
		}
		
		static void BuildMultipleUsers(std::vector<int> numItems,float gamma, float muy, const std::string& outfile)
		{
			std::stringstream ss;
			
			int numGates = 0;
			int numWires = 0;
			int numInputWires = 0;
			int numOutputWires = 0;
			
			for(int idx = 0; idx < numItems.size(); idx++)
			{
				numGates += (5 + 4*numItems[idx]);
				numWires += (6 + 7*numItems[idx]);
				numInputWires += (1 + 3*numItems[idx]);
			}
			numOutputWires = numItems.size();
			
			ss << "#gates #wires\n";
			ss << std::to_string(numGates) << " " << std::to_string(numWires) << "\n";
			
			ss << "InputSize OutputSize\n";
			ss << std::to_string(numInputWires) << " " << std::to_string(numOutputWires) << "\n";
			ss << "\n";
			
			int shift = 0;
			
			for(int kdx = 0; kdx < numItems.size(); kdx++)
			{
				ss << "#####---user " << kdx << "---#####\n";
				int count = 1 + 3*numItems[kdx];
				
				for(int idx = 0; idx < numItems[kdx]; idx++)
				{
					ss << "#####---item " << idx << "---#####\n";
					ss << "2 1 0 " << (shift + 1 + 3*idx) << " " << (shift + count) << " DOT\n"; count++;				// <u,v>
					ss << "2 1 " << (shift + 2 + 3*idx) << " " << (shift + count - 1) << " " << (shift + count) << " SUB\n"; count++;		// r - <u,v>
					ss << "2 1 " << (shift + count - 1) << " " << (shift + 1 + 3*idx) << " " << (shift + count) << " MUL\n"; count++;		// v(r - <u,v>)
					ss << "2 1 " << (shift + 3 + 3*idx) << " " << (shift + count - 1) << " " << (shift + count) << " MUL\n"; count++;		// isReal*v(r - <u,v>)
				}
				
				ss << "#####---Adding---#####\n";
				// Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)
				ss << numItems[kdx] << " 1 ";
				for(int idx = 1; idx <= numItems[kdx]; idx++)
				{
					ss << (shift + 3*numItems[kdx] + 4*idx) << " ";
				}
				ss << (shift + count) << " NADD\n"; 
				count++;
				
				int gamma_int = (int)(gamma*1048576.0);
				int muy_int = (int)(muy*1048576.0);
				
				ss << "2 1 " << gamma << " " << (shift + count - 1) << " " << shift + count << " CMUL\n"; count++;			// gamma*[Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)]
				ss << "2 1 " << muy << " " << shift << " " << (shift + count) << " CMUL\n"; count++;					// muy*u_i
				ss << "2 1 " << (shift + count - 2) << " " << (shift + count - 1) << " " << (shift + count) << " ADD\n"; count++;			// gamma*[Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)] + muy*u_i
				ss << "2 1 " << shift << " " << (shift + count - 1) << " " << (shift + count) << " ADD\n"; count++;				// u_i + gamma*[Sum isReal_ij*v_j*(r_ij - <u_i, v_j>)] + muy*u_i
				
				shift += count;
			}
			std::ofstream myfile;
			myfile.open(outfile);
			myfile << ss.str();
			myfile.close();
		}
		
		static ArithmeticCircuit * BuildMultipleUsersByLayers(std::vector<int> numItems, float gamma, float muy)
		{
			std::stringstream ss;
			
			int numGates = 0;
			int numWires = 0;
			int numInputWires = 0;
			int numOutputWires = 0;
			      
			for(int idx = 0; idx < numItems.size(); idx++)
			{
				numGates += (5 + 3*numItems[idx]);
				numWires += (6 + 6*numItems[idx]);
				numInputWires += (1 + 3*numItems[idx]);
			}
			numOutputWires = numItems.size();
			
			std::vector<IGate *> gates(numGates);
			
// 			ss << "#gates #wires\n";
// 			ss << std::to_string(numGates) << " " << std::to_string(numWires) << "\n";
// 			
// 			ss << "InputSize OutputSize\n";
// 			ss << std::to_string(numInputWires) << " " << std::to_string(numOutputWires) << "\n";
// 			ss << "\n";
			
			int numUsers = numItems.size(); // 0 --> n-1
			int sumNumItems = 0;		// t_0 + ... + t_{n-1}
			for(int idx = 0; idx < numUsers; idx++)
			{
				sumNumItems += numItems[idx];
			}
			
			int gateCount = 0;
			
			// input wires has values from 0 --> (numUsers - 1) + 3*sumNumItems
			
			// Generate Dot Product Gates for layer 0
			
			// shift for user 0
			int shift = 0;
			int count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << "2 1 " << shift << " " << (shift + 1 + 3*kdx) << " " << (numUsers + 3*sumNumItems + count) << " DOT\n";
					gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
					count++;
				}
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate SUB gates for layer 1
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << "2 1 " << (shift + 2 + 3*kdx) << " " << (numUsers + 3*sumNumItems + count) << " " << (numUsers + 4*sumNumItems + count) << " SUB\n";
					gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
					count++;
				}
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate MUL gates for layer 2
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << "2 1 " << (numUsers + 4*sumNumItems + count) << " " << (shift + 1 + 3*kdx) << " " << (numUsers + 5*sumNumItems + count) << " MUL\n";
					gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
					count++;
				}
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate NDOT gates for layer 3
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << 2*numItems[idx] << " 1 ";
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << (shift + 3 + 3*kdx) << " ";
					
				}
				
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << (numUsers + 5*sumNumItems + count) << " ";
					count++;
				}
				
				ss << (numUsers + 6*sumNumItems + idx) << " NDOT\n";
				gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
					
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate CMUL gates for layer 2 and 4
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << "2 1 " << muy << " " << shift << " " << (2*numUsers + 6*sumNumItems + idx) << " CMUL\n";
				gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
				
				ss << "2 1 " << gamma << " " << (numUsers + 6*sumNumItems + idx) << " " << (3*numUsers + 6*sumNumItems + idx) << " CMUL\n";
				gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
				
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate ADD gates for layer 5
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << "2 1 " << (2*numUsers + 6*sumNumItems + idx) << " " << (3*numUsers + 6*sumNumItems + idx) << " " << (4*numUsers + 6*sumNumItems + idx) << " ADD\n";
				gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
			}
			
			// Generate ADD gates for layer 6
			// shift for user 0
			shift = 0;
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << "2 1 " << shift << " " << (4*numUsers + 6*sumNumItems + idx) << " " << (5*numUsers + 6*sumNumItems + idx) << " ADD\n";
				gates[gateCount] = CreateGate(ss.str()); ss.str(""); gateCount++;
				
				shift += (1 + 3*numItems[idx]);
			}
			
			return new ArithmeticCircuit(numGates, numWires, numInputWires, numOutputWires, gates);
			// First layer: Dot product gates
			// user 0, item v_0: n + 3Sum t_1  
// 			std::ofstream myfile;
// 			myfile.open(outfile);
// 			myfile << ss.str();
// 			myfile.close();
		}
		
		static ArithmeticCircuit * BuildHistogramCircuit(std::vector<int> numItems, float gamma, float muy)
		{
			std::stringstream ss;
			
			int numGates = 0;
			int numWires = 0;
			int numInputWires = 0;
			int numOutputWires = 0;
			      
			for(int idx = 0; idx < numItems.size(); idx++)
			{
				numGates += 1;
				numWires += numItems[idx] + 1;
				numInputWires += numItems[idx];
			}
			
			numOutputWires = numItems.size();
			
			std::vector<IGate *> gates(numGates);
			
// 			ss << "#gates #wires\n";
// 			ss << std::to_string(numGates) << " " << std::to_string(numWires) << "\n";
// 			
// 			ss << "InputSize OutputSize\n";
// 			ss << std::to_string(numInputWires) << " " << std::to_string(numOutputWires) << "\n";
// 			ss << "\n";
			
			int numUsers = numItems.size(); // 0 --> n-1
			int sumNumItems = 0;		// t_0 + ... + t_{n-1}
			for(int idx = 0; idx < numUsers; idx++)
			{
				sumNumItems += numItems[idx];
			}
			
			int gateCount = 0;
			
			// input wires has values from 0 --> sumNumItems - 1
			
			// Generate NADD Gates for layer 0
			
			// shift for user 0
			int inputShift = 0;
			int outputShift = sumNumItems; // increment for output wire of each gate
			
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				// Number of input wires
				ss << numItems[idx] << " 1 ";
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << inputShift << " "; inputShift++;
				}
				
				ss << outputShift << " "; outputShift++;
				
				ss << "NADD\n";
				
				gates[idx] = CreateGate(ss.str()); 
				ss.str("");
			}
			
			
			
			return new ArithmeticCircuit(numGates, numWires, numInputWires, numOutputWires, gates);
		}
		
		static IGate * CreateGate(const std::string &str)
		{
			IGate *gates;
			
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
					gates = new AdditionGate(inputWireOrConstant, inputWire, outputWire);
				}
				else if (gate[0] == 'S')
				{
					gates = new SubtractionGate(inputWireOrConstant, inputWire, outputWire);
				}
				else if (gate[0] == 'M')
				{
					gates = new MultiplicationGate(inputWireOrConstant, inputWire, outputWire);
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
						gates = new NaryAdditionGate(inputWires, outputWire);
					}
					else if (gate[1] == 'D')
					{
						std::vector<int> inputWires;
						for(int idx = 0; idx < inputSize; idx++)
						{
							inputWires.push_back(f(2+idx));
						}
						int outputWire = f(2 + inputSize);
						gates = new NaryDotGate(inputWires, outputWire);
					}
				}
				else if (gate[0] == 'D' && gate[1] == 'O')
				{
					gates = new DotProductGate(inputWireOrConstant, inputWire, outputWire);
				}
				else if (gate[0] == 'D' && gate[1] == 'I')
				{
	// 				std::cout << "Div gate" << std::endl;
					gates = new DivisionGate(inputWireOrConstant, inputWire, outputWire);
				}
				else if (gate[1] == 'A')
				{
					assert(gate[0] == 'C');
					gates = new ConstantAdditionGate(inputWireOrConstant, inputWire, outputWire);
				}
				else if (gate[1] == 'M')
				{
					assert(gate[0] == 'C');
					float inputWireOrConstant = std::stof(variables[2]);
					gates = new ConstantMultiplicationGate(inputWireOrConstant, inputWire, outputWire);
				}
				else
				{
					std::cout << "Here : gate[0] = " << gate[0] << "\tgate[1] = " << gate[1] << std::endl;
					std::cout << "Gate " + gate + " does not exist." << std::endl;
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
					gates = new NaryAdditionGate(inputWires, outputWire);
				}
				else if (gate[1] == 'D')
				{
					std::vector<int> inputWires;
					for(int idx = 0; idx < inputSize; idx++)
					{
						inputWires.push_back(f(2+idx));
					}
					int outputWire = f(2 + inputSize);
					gates = new NaryDotGate(inputWires, outputWire);
				}
			}
			
			return gates;
		}
		
		static void BuildMultipleUsersByLayers(std::vector<int> numItems, float gamma, float muy, const std::string& outfile)
		{
			std::stringstream ss;
			
			int numGates = 0;
			int numWires = 0;
			int numInputWires = 0;
			int numOutputWires = 0;
			      
			for(int idx = 0; idx < numItems.size(); idx++)
			{
				numGates += (5 + 3*numItems[idx]);
				numWires += (6 + 6*numItems[idx]);
				numInputWires += (1 + 3*numItems[idx]);
			}
			numOutputWires = numItems.size();
			
			ss << "#gates #wires\n";
			ss << std::to_string(numGates) << " " << std::to_string(numWires) << "\n";
			
			ss << "InputSize OutputSize\n";
			ss << std::to_string(numInputWires) << " " << std::to_string(numOutputWires) << "\n";
			ss << "\n";
			
			int numUsers = numItems.size(); // 0 --> n-1
			int sumNumItems = 0;		// t_0 + ... + t_{n-1}
			for(int idx = 0; idx < numUsers; idx++)
			{
				sumNumItems += numItems[idx];
			}
			
			// input wires has values from 0 --> (numUsers - 1) + 3*sumNumItems
			
			// Generate Dot Product Gates for layer 0
			
			// shift for user 0
			int shift = 0;
			int count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << "2 1 " << shift << " " << (shift + 1 + 3*kdx) << " " << (numUsers + 3*sumNumItems + count) << " DOT\n";
					count++;
				}
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate SUB gates for layer 1
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << "2 1 " << (shift + 2 + 3*kdx) << " " << (numUsers + 3*sumNumItems + count) << " " << (numUsers + 4*sumNumItems + count) << " SUB\n";
					count++;
				}
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate MUL gates for layer 2
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << "2 1 " << (numUsers + 4*sumNumItems + count) << " " << (shift + 1 + 3*kdx) << " " << (numUsers + 5*sumNumItems + count) << " MUL\n";
					count++;
				}
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate NDOT gates for layer 3
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << 2*numItems[idx] << " 1 ";
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << (shift + 3 + 3*kdx) << " ";
					
				}
				
				for(int kdx = 0; kdx < numItems[idx]; kdx++)
				{
					ss << (numUsers + 5*sumNumItems + count) << " ";
					count++;
				}
				
				ss << (numUsers + 6*sumNumItems + idx) << " NDOT\n";
					
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate CMUL gates for layer 2 and 4
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << "2 1 " << muy << " " << shift << " " << (2*numUsers + 6*sumNumItems + idx) << " CMUL\n";
				ss << "2 1 " << gamma << " " << (numUsers + 6*sumNumItems + idx) << " " << (3*numUsers + 6*sumNumItems + idx) << " CMUL\n";
				shift += (1 + 3*numItems[idx]);
			}
			
			// Generate ADD gates for layer 5
			// shift for user 0
			shift = 0;
			count = 0; // increment for output wire of each gate
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << "2 1 " << (2*numUsers + 6*sumNumItems + idx) << " " << (3*numUsers + 6*sumNumItems + idx) << " " << (4*numUsers + 6*sumNumItems + idx) << " ADD\n";
			}
			
			// Generate ADD gates for layer 6
			// shift for user 0
			shift = 0;
			// Go through each user
			for(int idx = 0; idx < numUsers; idx++)
			{
				ss << "2 1 " << shift << " " << (4*numUsers + 6*sumNumItems + idx) << " " << (5*numUsers + 6*sumNumItems + idx) << " ADD\n";
				shift += (1 + 3*numItems[idx]);
			}
			
			// First layer: Dot product gates
			// user 0, item v_0: n + 3Sum t_1  
			std::ofstream myfile;
			myfile.open(outfile);
			myfile << ss.str();
			myfile.close();
		}
		
		
// 		static void BuildMultipleUsersByLayers(std::vector<int> numItems, float gamma, float muy, const std::string& outfile)
// 		{
// 			std::stringstream ss;
// 			
// 			int numGates = 0;
// 			int numWires = 0;
// 			int numInputWires = 0;
// 			int numOutputWires = 0;
// 			      
// 			for(int idx = 0; idx < numItems.size(); idx++)
// 			{
// 				numGates += (5 + 4*numItems[idx]);
// 				numWires += (6 + 7*numItems[idx]);
// 				numInputWires += (1 + 3*numItems[idx]);
// 			}
// 			numOutputWires = numItems.size();
// 			
// 			ss << "#gates #wires\n";
// 			ss << std::to_string(numGates) << " " << std::to_string(numWires) << "\n";
// 			
// 			ss << "InputSize OutputSize\n";
// 			ss << std::to_string(numInputWires) << " " << std::to_string(numOutputWires) << "\n";
// 			ss << "\n";
// 			
// 			int numUsers = numItems.size(); // 0 --> n-1
// 			int sumNumItems = 0;		// t_0 + ... + t_{n-1}
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				sumNumItems += numItems[idx];
// 			}
// 			
// 			// input wires has values from 0 --> (numUsers - 1) + 3*sumNumItems
// 			
// 			// Generate Dot Product Gates for layer 0
// 			
// 			// shift for user 0
// 			int shift = 0;
// 			int count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				for(int kdx = 0; kdx < numItems[idx]; kdx++)
// 				{
// 					ss << "2 1 " << shift << " " << (shift + 1 + 3*kdx) << " " << (numUsers + 3*sumNumItems + count) << " DOT\n";
// 					count++;
// 				}
// 				shift += (1 + 3*numItems[idx]);
// 			}
// 			
// 			// Generate SUB gates for layer 1
// 			// shift for user 0
// 			shift = 0;
// 			count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				for(int kdx = 0; kdx < numItems[idx]; kdx++)
// 				{
// 					ss << "2 1 " << (shift + 2 + 3*kdx) << " " << (numUsers + 3*sumNumItems + count) << " " << (numUsers + 4*sumNumItems + count) << " SUB\n";
// 					count++;
// 				}
// 				shift += (1 + 3*numItems[idx]);
// 			}
// 			
// 			// Generate MUL gates for layer 2
// 			// shift for user 0
// 			shift = 0;
// 			count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				for(int kdx = 0; kdx < numItems[idx]; kdx++)
// 				{
// 					ss << "2 1 " << (numUsers + 4*sumNumItems + count) << " " << (shift + 1 + 3*kdx) << " " << (numUsers + 5*sumNumItems + count) << " MUL\n";
// 					count++;
// 				}
// 				shift += (1 + 3*numItems[idx]);
// 			}
// 			
// 			// Generate MUL gates for layer 3
// 			// shift for user 0
// 			shift = 0;
// 			count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				for(int kdx = 0; kdx < numItems[idx]; kdx++)
// 				{
// 					ss << "2 1 " << (shift + 3 + 3*kdx) << " " << (numUsers + 5*sumNumItems + count) << " " << (numUsers + 6*sumNumItems + count) << " MUL\n";
// 					count++;
// 				}
// 				shift += (1 + 3*numItems[idx]);
// 			}
// 			
// 			// Generate Nary ADD gates for layer 4
// 			// shift for user 0
// 			shift = 0;
// 			count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				ss << numItems[idx] << " 1 ";
// 				for(int kdx = 0; kdx < numItems[idx]; kdx++)
// 				{
// 					ss << (shift + numUsers + 6*sumNumItems + kdx) << " ";
// 				}
// 				ss << (numUsers + 7*sumNumItems + idx) << " NADD\n";
// 				shift += numItems[idx];
// 			}
// 			
// 			// Generate CMUL gates for layer 2 and 5
// 			// shift for user 0
// 			shift = 0;
// 			count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				ss << "2 1 " << muy << " " << shift << " " << (2*numUsers + 7*sumNumItems + idx) << " CMUL\n";
// 				ss << "2 1 " << gamma << " " << (numUsers + 7*sumNumItems + idx) << " " << (3*numUsers + 7*sumNumItems + idx) << " CMUL\n";
// 				shift += (1 + 3*numItems[idx]);
// 			}
// 			
// 			// Generate ADD gates for layer 6
// 			// shift for user 0
// 			shift = 0;
// 			count = 0; // increment for output wire of each gate
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				ss << "2 1 " << (2*numUsers + 7*sumNumItems + idx) << " " << (3*numUsers + 7*sumNumItems + idx) << " " << (4*numUsers + 7*sumNumItems + idx) << " ADD\n";
// 			}
// 			
// 			// Generate ADD gates for layer 7
// 			// shift for user 0
// 			shift = 0;
// 			// Go through each user
// 			for(int idx = 0; idx < numUsers; idx++)
// 			{
// 				ss << "2 1 " << shift << " " << (4*numUsers + 7*sumNumItems + idx) << " " << (5*numUsers + 7*sumNumItems + idx) << " ADD\n";
// 				shift += (1 + 3*numItems[idx]);
// 			}
// 			
// 			// First layer: Dot product gates
// 			// user 0, item v_0: n + 3Sum t_1  
// 			std::ofstream myfile;
// 			myfile.open(outfile);
// 			myfile << ss.str();
// 			myfile.close();
// 		}
	};
#endif
