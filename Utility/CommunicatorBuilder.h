#ifndef COMMUNICATOR_BUILDER_H__
#define COMMUNICATOR_BUILDER_H__

#pragma once

#include "../Network.hpp"
#include "Communicator.h"

namespace Utility
{
	class CommunicatorBuilder
	{
	public:
		static void BuildCommunicator(int party, Machine *machine, Communicator * &c)
		{
			Network *cPartner, *cAlice, *cBob, *cCharlie;
			
// 			bool isAlice = true;
			
			if(Alice == party)
			{
				cBob = machine->partiesDown[Bob - Alice - 1];
				cCharlie   = machine->partiesDown[Charlie - Alice - 1];
			  
				std::vector<Network *> input = {cBob, cCharlie};
				std::vector<Network *> output = {cBob, cCharlie};
				
				c = new Communicator(party, input, output);
			}
			else if(Bob == party)
			{
				cAlice = machine->partiesUp[Alice];
				cCharlie   = machine->partiesDown[Charlie - Bob - 1];
			  
				std::vector<Network *> input = {cAlice, cCharlie};
				std::vector<Network *> output = {cAlice, cCharlie};
				
				c = new Communicator(party, input, output);
			}
			else if(Charlie == party)
			{
				cAlice   = machine->partiesUp[Alice];
				cBob     = machine->partiesUp[Bob];
			  
				std::vector<Network *> input = {cAlice, cBob};
				std::vector<Network *> output = {cAlice, cBob};
				
				c = new Communicator(party, input, output);
			}
		}
	};
}

#endif

