#ifndef COMMUNICATOR_H__
#define COMMUNICATOR_H__

#pragma once

#include "../Network.hpp"
#include <vector>

namespace Utility
{
	class Communicator
	{
	public:
		int party;
		uint64_t sendingSize;
		uint64_t receivingSize;
		
		std::vector<Network *> const inputChannels;
		std::vector<Network *> const outputChannels;

		static constexpr int PARTNER = 0;
		
		/// <summary>
		/// Order of input channels is Partner, Alice, Bob.
		/// </summary>
		/// <param name="isAlice"></param>
		/// <param name="inputChannels"></param>
		/// <param name="outputChannels"></param>
		
		Communicator(const int party, std::vector<Network *> &inputChannels, std::vector<Network *> &outputChannels);

		void SendMessage(int player, unsigned char *msg, unsigned int size);

		unsigned int AwaitMessage(int player, unsigned char *msg, unsigned int bufferSize);

		void SendEvaluationPartner(unsigned char *msg, unsigned int size);
		
		unsigned int AwaitEvaluationPartner(unsigned char *msg, unsigned int size);

		void SendAlice(unsigned char *msg, unsigned int size);

		unsigned int AwaitAlice(unsigned char *msg, unsigned int size);

		void SendBob(unsigned char *msg, unsigned int size);
		
		unsigned int AwaitBob(unsigned char *msg, unsigned int size);
		
		void SendCharlie(unsigned char *msg, unsigned int size);
		
		unsigned int AwaitCharlie(unsigned char *msg, unsigned int size);
		
		void SendNext(unsigned char *msg, unsigned int size);
		void SendPrevious(unsigned char *msg, unsigned int size);
		
		unsigned int AwaitNext(unsigned char *msg, unsigned int size);
		unsigned int AwaitPrevious(unsigned char *msg, unsigned int size);
		
		std::vector<unsigned char> getCommonSeed();
		std::vector<unsigned char> getSharedRandomSeed(int party, int partner);
	private:
		int ALICE;
		int BOB;
		int CHARLIE;
	};
}

#endif
