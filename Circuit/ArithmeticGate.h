#ifndef _ARITHMETIC_GATE_H__
#define _ARITHMETIC_GATE_H__

#include <stdio.h>
#include "Gate.h"
#include "math.h"
#include "../Utility/CryptoUtility.h"

namespace Circuit
{  
	class AdditionGate : public BinaryGate
	{
	public:
		/// Constructor
		AdditionGate(int leftwire, int rightwire, int outwire) : BinaryGate(leftwire, rightwire, outwire){}
		
		/// Destructor 
		~AdditionGate(){}
	};

	class SubtractionGate : public BinaryGate
	{
	public:
		/// Constructor
		SubtractionGate(int leftwire, int rightwire, int outwire) : BinaryGate(leftwire, rightwire, outwire){}
		
		/// Destructor 
		~SubtractionGate(){}
	};
	
	class MultiplicationGate : public BinaryGate
	{
	public:
		/// Constructor
		MultiplicationGate(int leftwire, int rightwire, int outwire) : BinaryGate(leftwire, rightwire, outwire){} 
		
		/// Destructor 
		~MultiplicationGate(){}
	};
	
	class DotProductGate : public BinaryGate
	{
	public:
		DotProductGate(int leftwire, int rightwire, int outwire) : BinaryGate(leftwire, rightwire, outwire){} 
		
		~DotProductGate(){};
	};
	
	class DivisionGate : public BinaryGate
	{
	public:
		/// Constructor
		DivisionGate(int leftwire, int rightwire, int outwire) : BinaryGate(leftwire, rightwire, outwire){} 
		
		/// Destructor 
		~DivisionGate(){}
	};
	
	class ConstantGate : public UnitaryGate
	{
	public:
		ConstantGate(float constant, int inputWire, int outputWire) : UnitaryGate(inputWire, outputWire)
		{
			this->constant = (int)(constant*pow(2.0, PRECISION_BIT_LENGTH));
		}
		
		~ConstantGate(){}
		
		std::string ToString()
		{
			return std::to_string(InputWire) + "|(" + std::to_string(constant/pow(2.0, PRECISION_BIT_LENGTH)) + ")|" + std::to_string(OutputWire);
		}
		
		int constant;
	};
	
	class ConstantAdditionGate : public ConstantGate
	{
	public:
		/// Constructor
		ConstantAdditionGate(float constant, int inputwire, int outputwire) : ConstantGate(constant, inputwire, outputwire){}
		
		/// Destructor
		~ConstantAdditionGate(){};
	};

	class ConstantMultiplicationGate : public ConstantGate
	{
	public:
		/// Constructor
		ConstantMultiplicationGate(float constant, int inputwire, int outputwire) : ConstantGate(constant, inputwire, outputwire){}
		
		/// Destructor
		~ConstantMultiplicationGate(){};
	};
	
	class NaryAdditionGate : public NaryGate
	{
	public:
		NaryAdditionGate(const std::vector<int>& inputWires, int outputWire) : NaryGate(inputWires, outputWire){} 
	    
		~NaryAdditionGate(){}
	};
	
	class NaryDotGate : public NaryGate
	{
	public:
		NaryDotGate(const std::vector<int>& inputWires, int outputWire) : NaryGate(inputWires, outputWire){} 
		~NaryDotGate(){}
	};
}

#endif