#include <cassert>
#include "Communicator.h"
#include "CryptoUtility.h"
#include <emp-tool/utils/block.h>

#define MAX_PACKAGE_SIZE 1073741824

namespace Utility
{
	Communicator::Communicator(const int party, std::vector<Network *> &inputChannels, std::vector<Network *> &outputChannels) : party(party), inputChannels(inputChannels), outputChannels(outputChannels){
		if(Alice == party)
		{
			ALICE = -1; BOB = 0; CHARLIE = 1;
		}
		else if(Bob == party)
		{
			ALICE = 0; BOB = -1; CHARLIE = 1;
		}
		else if(Charlie == party)
		{
			ALICE = 0; BOB = 1; CHARLIE = -1;
		}
		
		sendingSize = 0;
		receivingSize = 0;
	}

	void Communicator::SendMessage(int player, unsigned char *msg, unsigned int size)
	{
		sendingSize += size;
		
		if(size > 100*1024*1024) std::cout << "sending " << size << " bytes" << std::endl;
		
		unsigned int count = 0;
		unsigned int bufferSize;
		
		while(count < size)
		{
			if (count + MAX_PACKAGE_SIZE < size)
			{
				bufferSize = MAX_PACKAGE_SIZE;
			}
			else
			{
				bufferSize = size - count;
			}
			
			unsigned int bytes_sent = 0;

			while(bytes_sent < bufferSize)
			{
				unsigned int bs = send(inputChannels[player]->sock_fd, msg + count, bufferSize - bytes_sent, 0);

				if(bs < 0)
				{
					perror("error: send");
					exit(1);
				}

				bytes_sent += bs;
				count += bs;
			}
		}
		
		assert(count == size);
	}

	unsigned int Communicator::AwaitMessage(int player, unsigned char *msg, unsigned int size)
	{
		receivingSize += size;
		if(size > 100*1024*1024) std::cout << "awaiting " << size << " bytes" << std::endl;
		
		
		int sd = outputChannels[player]->sock_fd;
		fd_set reading_set;
		
		while (true)
		{
			FD_ZERO(&reading_set); //clear the socket set 
			FD_SET(sd, &reading_set); //add master socket to set
			int activity = select(sd+1, &reading_set , NULL , NULL , NULL);   //&timeout
			if ((activity < 0) && (errno!=EINTR))  { perror("select failed"); exit(1); }
			if (activity > 0)
			{
				if (FD_ISSET(sd, &reading_set)){
					unsigned int count = 0;
					unsigned int bufferSize;
					unsigned int bytes_received = 0;
					
					while(count < size)
					{
						if (count + MAX_PACKAGE_SIZE < size)
						{
							bufferSize = MAX_PACKAGE_SIZE;
						}
						else
						{
							bufferSize = size - count;
						}
						
						unsigned int arrived = 0;

						while (arrived < bufferSize){
							unsigned int br = recv(sd, msg + count, bufferSize - arrived, 0);
							arrived += br;
							count += br;
						}

						bytes_received += arrived;

						// count += MAX_PACKAGE_SIZE;
					}
					assert(bytes_received == size);
					return bytes_received;
				}
				
			}
			else if (activity == 0){
				return 0;
			}
		}
	}

	void Communicator::SendAlice(unsigned char *msg, unsigned int size)
	{
		SendMessage(ALICE, msg, size);
	}

	unsigned int Communicator::AwaitAlice(unsigned char *msg, unsigned int size)
	{
		return AwaitMessage(ALICE, msg, size);
	}

	void Communicator::SendBob(unsigned char *msg, unsigned int size)
	{
		SendMessage(BOB, msg, size);
	}
	
	unsigned int Communicator::AwaitBob(unsigned char *msg, unsigned int size)
	{
		return AwaitMessage(BOB, msg, size);
	}
	
	void Communicator::SendCharlie(unsigned char *msg, unsigned int size)
	{
		SendMessage(CHARLIE, msg, size);
	}

	unsigned int Communicator::AwaitCharlie(unsigned char *msg, unsigned int size)
	{
		return AwaitMessage(CHARLIE, msg, size);
	}
	
	void Communicator::SendNext(unsigned char *msg, unsigned int size)
	{
		if(party == Alice)
		{
			SendMessage(BOB, msg, size);
		}
		else if(party == Bob)
		{
			SendMessage(CHARLIE, msg, size);
		}
		else
		{
			SendMessage(ALICE, msg, size);
		}
	}
	
	void Communicator::SendPrevious(unsigned char *msg, unsigned int size)
	{
		if(party == Alice)
		{
			SendMessage(CHARLIE, msg, size);
		}
		else if(party == Bob)
		{
			SendMessage(ALICE, msg, size);
		}
		else
		{
			SendMessage(BOB, msg, size);
		}
	}

	unsigned int Communicator::AwaitPrevious(unsigned char *msg, unsigned int size)
	{
		if(party == Alice)
		{
			return AwaitMessage(CHARLIE, msg, size);
		}
		else if(party == Bob)
		{
			return AwaitMessage(ALICE, msg, size);
		}
		else
		{
			return AwaitMessage(BOB, msg, size);
		}
	}

	unsigned int Communicator::AwaitNext(unsigned char *msg, unsigned int size)
	{
		if(party == Alice)
		{
			return AwaitMessage(BOB, msg, size);
		}
		else if(party == Bob)
		{
			return AwaitMessage(CHARLIE, msg, size);
		}
		else
		{
			return AwaitMessage(ALICE, msg, size);
		}
	}
	
	void Communicator::SendEvaluationPartner(unsigned char *msg, unsigned int size)
	{
		assert(party != Charlie);
		SendMessage(PARTNER, msg, size);
	}

	unsigned int Communicator::AwaitEvaluationPartner(unsigned char *msg, unsigned int size)
	{
		assert(party != Charlie);
		return AwaitMessage(PARTNER, msg, size);
	}
	
	std::vector<unsigned char> Communicator::getCommonSeed()
	{
		// Generate random seed for shuffling
		std::vector<unsigned char> seed = CryptoUtility::SampleByteArray(32);
		std::vector<unsigned char> seed_prev(32), seed_next(32);
		
		for(int idx = 0; idx < 16; idx++) seed[idx] = emp::fix_key[idx];
		for(int idx = 16; idx < 32; idx++) seed[idx] = 0;
		return seed;
		
		if(party == Alice)
		{
			// Charlie - Alice - Bob
			SendNext((unsigned char *)(seed.data()), seed.size());
			SendPrevious((unsigned char *)(seed.data()), seed.size());
			AwaitNext((unsigned char *)(seed_next.data()), seed.size());
			AwaitPrevious((unsigned char *)(seed_prev.data()), seed.size());
		}
		else if(party == Bob)
		{
			// Alice - Bob - Charlie
			AwaitPrevious((unsigned char *)(seed_prev.data()), seed.size());
			SendNext((unsigned char *)(seed.data()), seed.size());
			SendPrevious((unsigned char *)(seed.data()), seed.size());
			AwaitNext((unsigned char *)(seed_next.data()), seed.size());
		}
		else if(party == Charlie)
		{
			// Bob - Charlie - Alice
			AwaitNext((unsigned char *)(seed_next.data()), seed.size());
			AwaitPrevious((unsigned char *)(seed_prev.data()), seed.size());
			SendNext((unsigned char *)(seed.data()), seed.size());
			SendPrevious((unsigned char *)(seed.data()), seed.size());
		}
		
		for(int idx = 0; idx < 32; idx++)
		{
			seed[idx] ^= (seed_prev[idx] ^ seed_next[idx]);
		}
		
		return seed;
	}
	
	std::vector<unsigned char> Communicator::getSharedRandomSeed(int party, int partner)
	{
		std::vector<unsigned char> seed(32), seed1, seed2;
		for(int idx = 0; idx < 16; idx++) seed[idx] = emp::fix_key[idx];
		for(int idx = 16; idx < 32; idx++) seed[idx] = 0;
		return seed;
		
		seed1 = CryptoUtility::SampleByteArray(32);
		seed2.resize(32);
		
		if(party > partner)
		{
			if(partner == Alice)
			{
				SendAlice(seed1.data(), seed1.size());
				AwaitAlice(seed2.data(), seed2.size());
			}
			else if(partner == Bob)
			{
				SendBob(seed1.data(), seed1.size());
				AwaitBob(seed2.data(), seed2.size());
			}
		}
		else
		{
			if(partner == Bob)
			{
				AwaitBob(seed2.data(), seed2.size());
				SendBob(seed1.data(), seed1.size());
			}
			else if(partner == Charlie)
			{
				AwaitCharlie(seed2.data(), seed2.size());
				SendCharlie(seed1.data(), seed1.size());
			}
		}
		
		for(int idx = 0; idx < 32; idx++)
		{
			seed[idx] = seed1[idx] + seed2[idx];
		}
		
		return seed;
	}
}
