#include <cassert>
#include <vector>
#include <iostream>
#include "Gate.h"
#include "ArithmeticGate.h"
#include "ArithmeticCircuit.h"
#include "InvalidOperationException.h"
#include "LayeredArithmeticCircuit.h"
#include "PlainArithmeticCircuitEvaluator.h"
#include "../Utility/CryptoUtility.h"

// #define PRECISION_BIT_LENGTH 20

using namespace Utility;

namespace Circuit
{
	void PlainArithmeticCircuitEvaluator::EvaluateArithmeticCircuit(ArithmeticCircuit *circuit, LayeredArithmeticCircuit *lc, const std::vector<int64_t>& input, std::vector<int64_t>& evaluation)
	{
		assert(circuit->InputLength == input.size());
		
		/// Value of input gates are known
		for(int idx = 0; idx < input.size(); idx++)
		{
			evaluation[idx] = input[idx];
		}
		
		/// Reset the rest of the wires to 0
		for(int idx = input.size(); idx < circuit->NumWires; idx++)
		{
			evaluation[idx] = 0;
		}
		
		for(int ldx = 0; ldx < lc->Depth; ldx++)
		{
			std::vector<IGate *> lgates;
			
			/// Get all the gates of the ldx-th layer
			for (auto &g : (lc->operator[](ldx))->AddGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->SubGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->MulGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->DivGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->CAddGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->CMulGates) {lgates.push_back(g);}
			
			/// Evaluate each gate
			for (int i = 0; i < lgates.size();i++)
			{
				IGate *g = lgates[i];
				
				if (dynamic_cast<BinaryGate*>(g) != nullptr)
				{
					BinaryGate *bg = (BinaryGate *)g;
					
					int64_t leftValue  = evaluation[bg->LeftWire];
					int64_t rightValue = evaluation[bg->RightWire];
					
					int outputWire = bg->OutputWire;
					
					if (dynamic_cast<AdditionGate*>(bg) != nullptr)
					{
						evaluation[outputWire] = ArithmeticOperation::add(leftValue, rightValue);
					}
					else if (dynamic_cast<SubtractionGate*>(bg) != nullptr)
					{
						evaluation[outputWire] = ArithmeticOperation::sub(leftValue, rightValue);
					}
					else if (dynamic_cast<MultiplicationGate*>(bg) != nullptr)
					{
						evaluation[outputWire] = ArithmeticOperation::mul(leftValue, rightValue);
					}
					else if (dynamic_cast<DivisionGate*>(bg) != nullptr)
					{
						assert(rightValue != 0);
						evaluation[outputWire] = ArithmeticOperation::div(leftValue, rightValue);
					}
					else
					{
						throw InvalidOperationException();
					}
				}
				else if (dynamic_cast<ConstantGate *>(g) != nullptr)
				{
					ConstantGate *ug = (ConstantGate *)g;
					
					int64_t inputValue = evaluation[ug->InputWire];
					int64_t constant = ug->constant;
					
					int outputWire = ug->OutputWire;

					if (dynamic_cast<ConstantAdditionGate*>(ug) != nullptr)
					{
						evaluation[outputWire] = ArithmeticOperation::add(constant, inputValue);
					}
					else if (dynamic_cast<ConstantMultiplicationGate*>(ug) != nullptr)
					{
						evaluation[outputWire] = ArithmeticOperation::mul(constant, inputValue);
					}
					else
					{
						throw InvalidOperationException();
					}
				}
			}
		}
	}
	
	void PlainArithmeticCircuitEvaluator::EvaluateArithmeticCircuit(ArithmeticCircuit *circuit, LayeredArithmeticCircuit *lc, const std::vector<std::vector<int64_t> >& input, const std::vector<int>& itemsPerUser, std::vector<std::vector<int64_t> >& evaluation)
	{
		assert(circuit->InputLength == input.size());
		
		evaluation.resize(circuit->NumWires);
		/// Value of input gates are known
		
		for(int idx = 0; idx < input.size(); idx++)
		{
			evaluation[idx] = input[idx];
		}
		
		for(int ldx = 0; ldx < lc->Depth; ldx++)
		{
			std::vector<IGate *> lgates;
			
			/// Get all the gates of the ldx-th layer
			for (auto &g : (lc->operator[](ldx))->AddGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->SubGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->MulGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->DotGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->NDotGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->NAddGates) {lgates.push_back(g);}
			for (auto &g : (lc->operator[](ldx))->CMulGates) {lgates.push_back(g);}
			
			/// Evaluate each gate
			for (int i = 0; i < lgates.size();i++)
			{
				IGate *g = lgates[i];
				
				if (dynamic_cast<BinaryGate*>(g) != nullptr)
				{
					BinaryGate *bg = (BinaryGate *)g;
					
					std::vector<int64_t> leftValue  = evaluation[bg->LeftWire];
					std::vector<int64_t> rightValue = evaluation[bg->RightWire];
					
					int outputWire = bg->OutputWire;
					
					if (dynamic_cast<AdditionGate*>(bg) != nullptr)
					{
// 						std::cout << "Add: " << leftValue.size() << " " << rightValue.size() << std::endl;
						evaluation[outputWire].resize(leftValue.size());
						for(int idx = 0; idx < leftValue.size(); idx++)
						{
							evaluation[outputWire][idx] = ArithmeticOperation::add(leftValue[idx], rightValue[idx]);
						}
					}
					else if (dynamic_cast<SubtractionGate*>(bg) != nullptr)
					{
						evaluation[outputWire].resize(rightValue.size());
						for(int idx = 0; idx < leftValue.size(); idx++)
						{
							evaluation[outputWire][idx] = ArithmeticOperation::sub(leftValue[idx], rightValue[idx]);
						}
					}
					else if (dynamic_cast<MultiplicationGate*>(bg) != nullptr)
					{
// 						std::cout << "Mult: " << leftValue.size() << " " << rightValue.size() << std::endl;
						if(leftValue.size() == 1){
							evaluation[outputWire].resize(rightValue.size());
							for(int idx = 0; idx < rightValue.size(); idx++)
							{
								evaluation[outputWire][idx] = ArithmeticOperation::mul(leftValue[0], rightValue[idx]);
							}
						}
						else if(rightValue.size() == 1){
							evaluation[outputWire].resize(leftValue.size());
							for(int idx = 0; idx < leftValue.size(); idx++)
							{
								evaluation[outputWire][idx] = ArithmeticOperation::mul(leftValue[idx], rightValue[0]);
							}
						}
					}
					else if (dynamic_cast<DotProductGate *>(bg) != nullptr)
					{
// 						std::cout << "Dot Product: " << bg->LeftWire << " " << bg->RightWire << " " << leftValue.size() << " " << rightValue.size() << std::endl;
						evaluation[outputWire].resize(1);
						evaluation[outputWire][0] = 0;
						for(int idx = 0; idx < rightValue.size(); idx++)
						{
							evaluation[outputWire][0] += leftValue[idx]*rightValue[idx];
						}
						evaluation[outputWire][0] = (evaluation[outputWire][0] >> PRECISION_BIT_LENGTH);
					}
					else
					{
						throw InvalidOperationException();
					}
				}
				else if (dynamic_cast<ConstantGate *>(g) != nullptr)
				{
					ConstantGate *ug = (ConstantGate *)g;
					
					std::vector<int64_t> inputValue = evaluation[ug->InputWire];
					int64_t constant = ug->constant;
					
					int outputWire = ug->OutputWire;
					evaluation[outputWire].resize(inputValue.size());
					
					if (dynamic_cast<ConstantAdditionGate*>(ug) != nullptr)
					{
						for(int idx = 0; idx < inputValue.size(); idx++)
						{
							evaluation[outputWire][idx] = ArithmeticOperation::add(constant, inputValue[idx]);
						}
					}
					else if (dynamic_cast<ConstantMultiplicationGate*>(ug) != nullptr)
					{
// 						std::cout << "CMUL " << (int)constant << std::endl;
						for(int idx = 0; idx < inputValue.size(); idx++)
						{
							evaluation[outputWire][idx] = ArithmeticOperation::mul(constant, inputValue[idx]);
						}
					}
					else
					{
						throw InvalidOperationException();
					}
				}
				else if (dynamic_cast<NaryAdditionGate *>(g) != nullptr)
				{
					NaryAdditionGate *ng  = (NaryAdditionGate *)g;
					int outputWire = ng->OutputWire;
					evaluation[outputWire].resize(dimension);
					
					std::fill(evaluation[outputWire].begin(), evaluation[outputWire].end(), 0);
					
					for(int idx = 0; idx < dimension; idx++)
					{
						for(int kdx = 0; kdx < ng->InputWires.size(); kdx++)
						{
							evaluation[outputWire][idx] += evaluation[ng->InputWires[kdx]][idx];
						}
					}
				}
				else if (dynamic_cast<NaryDotGate *>(g) != nullptr)
				{
					NaryDotGate *ng  = (NaryDotGate *)g;
					int outputWire = ng->OutputWire;
					evaluation[outputWire].resize(dimension);
					
					std::fill(evaluation[outputWire].begin(), evaluation[outputWire].end(), 0);
					
					int step = ng->InputWires.size()/2;
					
					for(int kdx = 0; kdx < step; kdx++)
					{
						evaluation[outputWire] = VectorOperation::Add(VectorOperation::Mul(evaluation[ng->InputWires[kdx]][0], evaluation[ng->InputWires[kdx + step]]), evaluation[outputWire]);
					}
				}
			}
		}
	}
}
