#ifndef BOOLEAN_GATE_H__
#define BOOLEAN_GATE_H__

#pragma once

#include "Gate.h"


namespace Circuit
{
	class BinaryGate : public IGate
	{
	public:
		~BinaryGate(){}
	};
	
	class AndGate : public BinaryGate
	{
	public:
		AndGate(int LeftWire, int RightWire, int OutputWire);
	};

	class XorGate : public BinaryGate
	{
	public:
		XorGate(int LeftWire, int RightWire, int OutputWire);
	};
	
	class NotGate : public UnitaryGate
	{
	public:
		NotGate(int inputWire, int outputWire);
	};
}

#endif
