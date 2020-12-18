#ifndef GATE_H__
#define GATE_H__

#pragma once

#include <vector>
#include <string>

namespace Circuit
{
	class IGate
	{
	public:
	    virtual ~IGate(){};
	    std::string ToString(){};
	};

	class UnitaryGate : public IGate
	{
	public:
		std::vector<uint64_t> InputWire;
		std::vector<uint64_t> OutputWire;

		UnitaryGate(int inputWire, int outputWire);
		UnitaryGate(const UnitaryGate &obj);
		UnitaryGate(const IGate &obj);
		
		std::string ToString();
		~UnitaryGate(){};
	};

	class BinaryGate : public IGate
	{
	public:
		int LeftWire;
		int RightWire;
		int OutputWire;

		/// <summary>
		/// Abstract Instantiation an arithmetic gate.
		/// </summary>
		BinaryGate(int LeftWire, int RightWire, int OutputWire);
		BinaryGate(const BinaryGate &obj);
		BinaryGate(const IGate &obj);
		
		std::string ToString();
		~BinaryGate(){};
	};
	
	class NaryGate : public IGate
	{
	public:
		std::vector<int> InputWires;
		int OutputWire;
		
		NaryGate(std::vector<int> inputWires, int outputWire);
		std::string ToString();
		~NaryGate(){};
	};
}

#endif
