#ifndef GRAPH_H__
#define GRAPH_H__

#pragma once
#include <cstddef>
#include <iterator>
#include <sstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include <cstring>
#include <openssl/rand.h>
#include <openssl/sha.h>
#include <bitset>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>

// #ifdef UNIX_PLATFORM

#include <unistd.h>
#include <arpa/inet.h>
#include <sys/types.h>
#include <netinet/tcp.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <errno.h>
#include <sys/ioctl.h>
#include <sys/poll.h>
#include <sys/time.h>

#include "Circuit/Gate.h"

#include "Utility/Range.h"
#include "Utility/Timer.h"
#include "Utility/Communicator.h"
#include "Utility/Commitment.h"
#include "Utility/CommunicatorBuilder.h"
#include "Utility/ISecureRNG.h"
#include "Utility/CryptoUtility.h"
#include "Utility/Stringhelper.h"
#include "Shares.h"
#include "Wires.hpp"
#include "ReplicatedArithmeticShares.h"
#include "ReplicatedRingShares.h"
#include "ReplicatedNewRingShares.h"

#include "Triplets.h"
#include "Machine.hpp"
#include "GraphData.hpp"
#include "Polynomial.h"

// #include "emp-ot/emp-ot/emp-ot.h"
// #include <emp-tool/emp-tool.h>

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_limbs.h>
#include <NTL/BasicThreadPool.h>

#include <emmintrin.h>

// using namespace emp;
using namespace std;
using namespace NTL;
using namespace Utility;
using namespace Circuit;

#define PSI 1

class Graph 
{ 
public:
	int machineId;
	int party;
	int partner;
	int totalItems;
	int totalParties;
	float rate;
	Machine *machine;
	Communicator *com, *comRecv;
	
	ZZ p;
	
	Graph (int machineId, int party, int totalItems, Machine *machine, float rate) {
		this->totalParties = 3;
		this->machineId = machineId;
		this->party = party;
		this->totalItems = totalItems;
		this->machine = machine;
		this->rate = rate;
		
		if(party == Alice) partner = Bob;
		if(party == Bob)   partner = Alice;
		
		CommunicatorBuilder::BuildCommunicator(party, machine, com);
		CommunicatorBuilder::BuildCommunicator(party, machine, comRecv);
		
		p = (ZZ(1) << 80) - 65; // n = 2^20
		ZZ_p::init(p);
	}

	~Graph(){
		delete machine;
	}
	
	void GraphComputation(int test) {
		int nItems = totalItems;
		Timer t;
// 		TestRSAMultiEval(1024, 2048);
		
		if(test == 0)
		{
			PSI_CA_Polynomial(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("PSI_CA_Polynomial************************************");
		}
		else if(test == 1)
		{
			PSI_CA_Circuit(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("PSI_CA_Circuit***************************************");
		}
		else if(test == 2)
		{
			PSI_CA_Hybrid(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("PSI_CA_Hybrid****************************************");
		}
		else if(test == 3)
		{
			fPSI(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("fPSI*************************************************");
		}
		else if(test == 4)
		{
			SA_PSI(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("SA_PSI***********************************************");
		}
		else if(test == 5)
		{
			PSI_Circuit(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("PSI_Circuit******************************************");
		}
		else if(test == 6)
		{
			fPSI_MCS(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("fPSI*************************************************");
		}
		else if(test == 61)
		{
			fPSI_MCA(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("fPSI*************************************************");
		}
		else if(test == 7)
		{
			int pl = rate;
			rate = 0.5;
			
			PSI_CA_Circuit_Payload(party, nItems, pl, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("PSI_CA_Circuit***************************************");
		}
		else if(test == 8)
		{
			PSI_CA_Circuit_Arithmetic(party, nItems, com, machine, rate);
			std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
			com->sendingSize = 0; com->receivingSize = 0;
			t.Tick("PSI_CA_Circuit***************************************");
		}
		else
		{
			std::cout << "0 <= test no <= 5" << std::endl;
		}
		
// 		fPSI_soPRF(party, nItems, com, machine, rate);
		
// 		t.Tick("fPSI_soPRF*******************************************");
		
// 		PSI_CA_Cristefaro(party, nItems, com, machine, rate);
// 		std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
	}
	
	void PSI_3PC_Circuit(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		
	}
	
	void SA_PSI(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		
		int fieldSize = 80; int nBytes = 16;
// 		ZZ pp = (ZZ(1) << 88) - 299; // n = 2^20
// 		ZZ_p::init(pp);
		
		std::vector<uint64_t> data;
		std::vector<__uint128_t> dataFx, dataPSI;
		
		std::cout << "Generate input" << std::endl;
		GenerateInput(party, data, nItems, rate);
		
		t.Tick("Generate input + dummies");
		
		int psiSize = 0;
		
		// Alice & Bob Compute f(k, x) and send to Charlie
		if(party == Alice || party == Bob)
		{
			std::vector<unsigned char> rand  = machine->getSharedRandomSeed(party, partner, com);
			AESRNG *rng = new AESRNG(rand.data());
			rng->ComputeRingFx(dataFx, data);
			
			com->SendCharlie((unsigned char *)(dataFx.data()), nBytes*nItems);
// 			t.Tick("Encode data");
			
			com->AwaitCharlie((unsigned char *)(&psiSize), sizeof(int));
			
			std::vector<__uint128_t> PSISet(psiSize), out(psiSize);
			
			com->AwaitCharlie((unsigned char *)(PSISet.data()), nBytes*psiSize);
			
			std::vector<unsigned char> hash = CryptoUtility::ComputeHash((unsigned char *)(PSISet.data()), nBytes*psiSize);
			std::vector<unsigned char> theirHash(hash.size());
			
			if(party == Alice)
			{
				com->SendBob(hash.data(), hash.size());
				com->AwaitBob(theirHash.data(), hash.size());
			}
			else
			{
				com->AwaitAlice(theirHash.data(), hash.size());
				com->SendAlice(hash.data(), hash.size());
			}
			
			assert(hash == theirHash);
			
			rng->ComputeRingFx(PSISet, out);
			
			std::vector<std::string> str;
			char *ptr2 = (char *)(out.data());
			
			for(int idx = 0; idx < psiSize; idx++)
			{
				std::string ss((const char *)ptr2, nBytes);
				str.push_back(ss);
				ptr2 += nBytes;
			}
			
			std::unordered_map<std::string, int> freq; 
			for (int idx = 0; idx < psiSize; idx++) 
			    freq[str[idx]]++; 
		      
			int count = 0; 
			for (auto it = freq.begin(); it != freq.end(); it++)
			{
				if(it->second == 1 || it->second == 3) count++;
			}
			
			std::cout << count << std::endl;
		}
		else if(party == Charlie)
		{
			// Charlie count the size of the union and send it to Alice & Bob 
			std::vector<__uint128_t> F_k1_xy(2*nItems);
			
			com->AwaitAlice((unsigned char *)(F_k1_xy.data()), nBytes*nItems);
			com->AwaitBob((unsigned char *)(F_k1_xy.data() + nItems), nBytes*nItems);
			
			dataFx.resize(2*nItems);

// 			t.Tick("Encode data");
			std::vector<__uint128_t> PSISet;
			
			std::vector<std::string> str;
			char *ptr2 = (char *)(F_k1_xy.data());
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				std::string ss((const char *)ptr2, nBytes);
				str.push_back(ss);
				ptr2 += nBytes;
			}
			
			std::unordered_map<std::string, int> freq; 
			for (int idx = 0; idx < 2*nItems; idx++) 
			    freq[str[idx]]++; 
		      
			std::vector<std::string> psiStr;
			__uint128_t temp;
			
			unsigned char * cstr = new unsigned char [nBytes + 1];
			int count = 0;
			
			for (auto itr=freq.begin(); itr!=freq.end(); itr++) 
			{ 
			    // if frequency is more than 1 
			    // print the element 
			    if (itr->second > 1) 
			    { 
				psiStr.push_back(itr->first);   
				std::strcpy ((char *)cstr, (itr->first).c_str());
				temp = *((__uint128_t *)cstr);
				PSISet.push_back(temp);
			    } 
			}
			
			psiSize = PSISet.size();
			com->SendAlice((unsigned char *)(&psiSize), sizeof(int));
			com->SendBob((unsigned char *)(&psiSize), sizeof(int));
			
			std::cout << "PSI Size: " << psiSize << std::endl;
			
			com->SendAlice((unsigned char *)(PSISet.data()), nBytes*psiSize);
			com->SendBob((unsigned char *)(PSISet.data()), nBytes*psiSize);
			
			t.Tick("PSI");
			std::cout << "Done sending PSI" << std::endl;
		}
		
		t.Tick("PSI");
	}
	
	
	void SA_PSI2(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		
		int fieldSize = 88; int nBytes = 11;
		ZZ pp = (ZZ(1) << 88) - 299; // n = 2^20
		ZZ_p::init(pp);
		
		std::cout << "nItems: " << nItems << "\tnBytes: " << nBytes << "\t" << ZZ_p::ModulusSize() << std::endl;
		
		std::vector<uint64_t> data;
		std::vector<ZZ_p> dataFx, dataPSI;
		
		std::cout << "Generate input" << std::endl;
		GenerateInput(party, data, nItems, rate);
		
		t.Tick("Generate input + dummies");
		
		int psiSize = 0;
		unsigned char *bytes = new unsigned char[nBytes*2*nItems];
		
		// Alice & Bob Compute f(k, x) and send to Charlie
		if(party == Alice || party == Bob)
		{
			std::vector<unsigned char> rand  = machine->getSharedRandomSeed(party, partner, com);
			AESRNG *rng = new AESRNG(rand.data());
			rng->ComputeFx(dataFx, data);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, dataFx.data(), nBytes, nItems);
			com->SendCharlie(bytes, nBytes*nItems);
			t.Tick("Encode data");
			
			com->AwaitCharlie((unsigned char *)(&psiSize), sizeof(int));
			
			std::cout << "PSI Size: " << psiSize << std::endl;
			
			unsigned char *psi = new unsigned char[nBytes*psiSize];
			
			com->AwaitCharlie(psi, nBytes*psiSize);
			
			std::vector<unsigned char> hash = CryptoUtility::ComputeHash(psi, nBytes*psiSize);
			std::vector<unsigned char> theirHash(hash.size());
			
			if(party == Alice)
			{
				com->SendBob(hash.data(), hash.size());
				com->AwaitBob(theirHash.data(), hash.size());
			}
			else
			{
				com->AwaitAlice(theirHash.data(), hash.size());
				com->SendAlice(hash.data(), hash.size());
			}
			
			assert(hash == theirHash);
			
			std::vector<ZZ_p> PSISet(psiSize), out(psiSize);
			ArrayEncoder::ByteArray2ZZArray(PSISet.data(), psi, nBytes, psiSize);
			
			rng->ComputeFx(PSISet, out, nBytes);
			
			delete [] psi;
		}
		else if(party == Charlie)
		{
			// Charlie count the size of the union and send it to Alice & Bob 
			unsigned char *ptr = bytes;
			
			com->AwaitAlice(bytes, nBytes*nItems);
			com->AwaitBob((bytes + nBytes*nItems), nBytes*nItems);
			
			dataFx.resize(2*nItems);

			std::vector<ZZ> F_k1_xy(2*nItems);
			
			ArrayEncoder::ByteArray2ZZArray(F_k1_xy.data(), bytes, nBytes, 2*nItems);
			t.Tick("Encode data");
			// Sort array according to f(k,x) and remember the permutation
			std::cout << "Sort data based on f(k,x)" << std::endl;
			
			// Sort
// 			std::sort(F_k1_xy.begin(), F_k1_xy.end());
// 			
// 			t.Tick("Sort");
			
			std::vector<ZZ> PSISet;
			
// 			for(int idx = 0; idx < F_k1_xy.size(); idx++)
// 			{
// 				if(F_k1_xy[idx] == F_k1_xy[idx+1]) PSISet.push_back(F_k1_xy[idx]);
// 			}
// 			t.Tick("Push into vector");
			
// 			PSISet.resize(0);
			
			std::vector<std::string> str;
			char *ptr2 = (char *)bytes;
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				std::string ss((const char *)ptr2, nBytes);
				str.push_back(ss);
				ptr2 += nBytes;
			}
			
			std::unordered_map<std::string, int> freq; 
			for (int idx = 0; idx < 2*nItems; idx++) 
			    freq[str[idx]]++; 
		      
			std::vector<std::string> psiStr;
			ZZ temp;
			
			unsigned char * cstr = new unsigned char [nBytes + 1];
			int count = 0;
			
			for (auto itr=freq.begin(); itr!=freq.end(); itr++) 
			{ 
			    // if frequency is more than 1 
			    // print the element 
			    if (itr->second > 1) 
			    { 
				psiStr.push_back(itr->first);   std::strcpy ((char *)cstr, (itr->first).c_str());
				ZZFromBytes(temp, cstr, nBytes);
				PSISet.push_back(temp);
			    } 
			} 
			t.Tick("Duplicates with hashmap");
			
			
			
			psiSize = PSISet.size();
			com->SendAlice((unsigned char *)(&psiSize), sizeof(int));
			com->SendBob((unsigned char *)(&psiSize), sizeof(int));
			
			std::cout << "PSI Size: " << psiSize << std::endl;
			
			ArrayEncoder::ZZArray2ByteArray(bytes, PSISet.data(), nBytes, psiSize);
			com->SendAlice(bytes, nBytes*psiSize);
			com->SendBob(bytes, nBytes*psiSize);
		}
		
		t.Tick("PSI");
	}
	
	void ComputeSOPRF2(int party, std::vector<uint64_t>& data, RepArithShares& soPRF, int nItems, int bitlength, Communicator *com)
	{
		Timer t;
		
		int nBytes;
		int length = data.size();
		
		MPCOperations mpc;
		
		ZZ_pContext context;
		context.save();
		// Distribute (shared/whole) prf keys
		// 1024 bit prime number
		ZZ p = conv<ZZ>("124325339146889384540494091085456630009856882741872806181731279018491820800119460022367403769795008250021191767583423221479185609066059226301250167164084041279837566626881119772675984258163062926954046545485368458404445166682380071370274810671501916789361956272226105723317679562001235501455748016154805420913");
		
		ZZ q = conv<ZZ>("1399252811935680595399801714158014275474696840019");
		ZZ g = conv<ZZ>("115740200527109164239523414760926155534485715860090261532154107313946218459149402375178179458041461723723231563839316251515439564315555249353831328479173170684416728715378198172203100328308536292821245983596065287318698169565702979765910089654821728828592422299160041156491980943427556153020487552135890973413");
		
		ZZ_p::init(q); nBytes = 20;
		
		RepArithShares r0; 
		
		std::vector<RepArithShares> sharedPrfKeys(bitlength);
		std::vector<RepArithShares> dataShares(bitlength);
		std::vector<RepArithShares> dataSharesBob(bitlength);
		
		std::cout << "Generate random keys" << std::endl;
		// Generate shared keys
		mpc.GenerateRandomShares(r0, nBytes, 1, com);
		
		for(int bdx = 0; bdx < bitlength; bdx++)
		{
			mpc.GenerateRandomShares(sharedPrfKeys[bdx], nBytes, 1, com);
		}
		
		std::cout << "Secret share data" << std::endl;
		// Secret share data
		std::vector<ZZ_p> temp;
		
		if(party == Alice || party == Bob)
		{
			std::vector<std::vector<ZZ_p> > dataInBits(bitlength);
			for(int bdx = 0; bdx < bitlength; bdx++) dataInBits[bdx].resize(nItems);
			
			for(int idx = 0; idx < nItems; idx++)
			{
				for(int bdx = 0; bdx < bitlength; bdx++)
				{
					dataInBits[bdx][idx] = ((data[idx] & (1ULL << bdx)) >> bdx);
				}
			}
			
			for(int bdx = 0; bdx < bitlength; bdx++)
			{
				if(party == Alice)
				{
					mpc.MPC_ShareData(dataShares[bdx], Alice, dataInBits[bdx], nBytes, nItems, com);
					mpc.MPC_ShareData(dataSharesBob[bdx], Bob, temp, nBytes, nItems, com);
				}
				else if(party == Bob)
				{
					mpc.MPC_ShareData(dataShares[bdx], Alice, temp, nBytes, nItems, com);
					mpc.MPC_ShareData(dataSharesBob[bdx], Bob, dataInBits[bdx], nBytes, nItems, com);
				}
			}
		}
		else
		{
			for(int bdx = 0; bdx < bitlength; bdx++)
			{
				mpc.MPC_ShareData(dataShares[bdx], Alice, temp, nBytes, nItems, com);
				mpc.MPC_ShareData(dataSharesBob[bdx], Bob, temp, nBytes, nItems, com);
			}
		}
		
		t.Tick("Secret Share data");
		
		std::cout << "Compute soPRF" << std::endl;
		// Verify inputs (may not need)
		
		// Compute (r_i+1)^x_i = r_i*x_i + 1;
		// We assume that the sampled keys are r_i + 1, thus no need to subtract 1 anymore
		
		for(int bdx = 0; bdx < bitlength; bdx++)
		{
			dataShares[bdx] = mpc.repMul(dataShares[bdx], sharedPrfKeys[bdx], nBytes, com);
			dataSharesBob[bdx] = mpc.repMul(dataSharesBob[bdx], sharedPrfKeys[bdx], nBytes, com);
			
			if(party == Alice)
			{
				for(int idx = 0; idx < nItems; idx++)
				{
					dataShares[bdx].x[idx] += 1;
					dataSharesBob[bdx].x[idx] += 1;
				}
			}
			else if(party == Bob)
			{
				for(int idx = 0; idx < nItems; idx++)
				{
					dataShares[bdx].y[idx] += 1;
					dataSharesBob[bdx].y[idx] += 1;
				}
			}
		}
		
		for(int bdx = bitlength - 2; bdx >= 0; bdx--)
		{
			dataShares[bdx] = mpc.repMul(dataShares[bdx], dataShares[bdx + 1], nBytes, com);
			dataSharesBob[bdx] = mpc.repMul(dataSharesBob[bdx], dataSharesBob[bdx + 1], nBytes, com);
		}
		
		dataShares[0] = mpc.repMul(dataShares[0], r0, nBytes, com);
		dataSharesBob[0] = mpc.repMul(dataSharesBob[0], r0, nBytes, com);
		
		dataShares[0].x.insert(dataShares[0].x.end(), /*std::make_move_iterator*/(dataSharesBob[0].x.begin()), /*std::make_move_iterator*/(dataSharesBob[0].x.end()));
		dataShares[0].y.insert(dataShares[0].y.end(), /*std::make_move_iterator*/(dataSharesBob[0].y.begin()), /*std::make_move_iterator*/(dataSharesBob[0].y.end()));
		
		context.save(); 
		ZZ_p::init(p); nBytes = 128;
		
		soPRF.x.resize(dataShares[0].x.size());
		soPRF.y.resize(dataShares[0].x.size());
		t.Tick("reset");
		std::vector<ZZ_p> table(160);
		table[0] = conv<ZZ_p>(g);
		for(int idx = 1; idx < 160; idx++)
		{
			table[idx] = table[idx-1]*table[idx-1];
		}
		
		for(int idx = 0; idx < dataShares[0].x.size(); idx++)
		{
			soPRF.x[idx] = 1;
			soPRF.y[idx] = 1;
			for(int kdx = 0; kdx < 160; kdx++)
			{
				if(IsOdd(rep(dataShares[0].x[idx]) >> kdx)) mul(soPRF.x[idx], soPRF.x[idx], table[kdx]);
				if(IsOdd(rep(dataShares[0].y[idx]) >> kdx)) mul(soPRF.y[idx], soPRF.y[idx], table[kdx]);
			}
		}
		t.Tick("exp");
		
		std::cout << "Reconstruct PRF" << std::endl;
		std::vector<ZZ_p> prf = mpc.openAllParties(soPRF, nBytes, com);
		for(int idx = 0; idx < 10; idx++) std::cout << idx << ": " << prf[idx] << std::endl;
	}
	
	void fPSI_soPRF(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		int bitlength = 64;
		std::vector<uint64_t> data;
		RepArithShares soPRF;
		
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
		
		ComputeSOPRF2(party, data, soPRF, nItems, bitlength, com);
	}
	
	// Compute a function over the intersection
	void fPSI_MCS(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		
		int nBytes = 10; int bitlength = rate; nBytes = bitlength/8; rate = 0.5;
		
		unsigned char *bytes = new unsigned char[nBytes*nItems];
		
		
		std::vector<Shares> sx(bitlength), sy(bitlength);
		
		std::vector<uint64_t> data(nItems*nBytes/8);
		
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems*nBytes/8, rate);
			
			if(party == Bob)
			{
				std::reverse(data.begin(), data.end());
			}
			
			// Transpose bits
			sse_trans((uint8_t *)(data.data()), bytes, nItems, bitlength);
			
			int numBlks = nItems/64;
			
			uint64_t *ptr = (uint64_t *)bytes;
			
			std::cout << "Secret share" << std::endl;
			for(int idx = 0; idx < bitlength; idx++)
			{
				std::vector<uint64_t> vec(ptr, ptr + numBlks);
				
				sx[idx].SecretShare(vec, nItems, com, Alice);
				sy[idx].SecretShare(vec, nItems, com, Bob);
				
				ptr += numBlks;
			}
		}
		else
		{
			std::cout << "Secret share" << std::endl;
			for(int idx = 0; idx < bitlength; idx++)
			{
				std::vector<uint64_t> vec(nItems/64);
				sx[idx].SecretShare(vec, nItems, com, Alice);
				sy[idx].SecretShare(vec, nItems, com, Bob);
			}
		}
		
		// Verify inputs
		Triplets tps;
		std::vector<Shares> sx2(bitlength), sy2(bitlength);
		for(int idx = 0; idx < bitlength; idx++)
		{
			sx2[idx] = sx[idx].LeftShiftOne();
			sy2[idx] = sy[idx].RightShiftOne();
		}
		
		if(tps.IsIncreasing(sx2, sx, com)) std::cout << "Passing input 1 check" << std::endl;
		if(tps.IsIncreasing(sy2, sy, com)) std::cout << "Passing input 2 check" << std::endl;
		
// 		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
// 		tpsInput.VerifyTriples(rngseed, com);
// 		
		t.Tick("Secret share inputs");
		
// 		Triplets tps;
		
		std::cout << "size: " << (1 + log2(nItems)) << std::endl;
		// CompareAndSwap
		for(int idx = 0; idx < (1 + log2(nItems)); idx++)
		{
			tps.CompareAndSwap(sx, sy, com);
		}
		
		std::vector<Shares> sx1 = sx;
		
// 		tps.VerifyTriples(rngseed, com);
		
		t.Tick("CompareAndSwap");
		
		// Pairwise Comparison
		// Merge sx and sy
		for(int idx = 0; idx < sx.size(); idx++)
		{
			sx[idx].s.insert(sx[idx].s.end(), sy[idx].s.begin(), sy[idx].s.end());
			sx[idx].t.insert(sx[idx].t.end(), sy[idx].t.begin(), sy[idx].t.end());
		}
		
		
		for(int idx = 0; idx < bitlength; idx++)
		{
			sx2[idx] = sx[idx].LeftShiftOne();
// 			sx3[idx] = sx2[idx].LeftShiftOne();
		}
		
		tps.PairwiseComparison(sx, sx2, com);
		t.Tick("PairwiseComparison");
		
		std::cout << sx1[0].s.size() << " " << sy[0].s.size() << std::endl;
		for(int idx = 0; idx < 2*log2(2*nItems); idx++)
		{
			tps.Swapping(sx1, sy, com);
		}
		
		t.Tick("Shuffle");
		std::vector<unsigned char> seed = com->getCommonSeed();
		tps.VerifyTriples(seed, com);
		
		t.Tick("Verify Triplets");
		
		std::cout << "Comm fpsi: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		std::vector<std::vector<uint64_t> > sout(sx.size());
// 		
// 		for(int idx = 0; idx < sx.size(); idx++)
// 		{
// 			if(party == Alice)
// 			{
// 				sout[idx] = sx[idx].OpenOneParty(com, Alice);
// 				sx[idx].OpenOneParty(com, Bob);
// 			}
// 			else if(party == Bob)
// 			{
// 				sx[idx].OpenOneParty(com, Alice);
// 				sout[idx] = sx[idx].OpenOneParty(com, Bob);
// 			}
// 			else
// 			{
// 				sx[idx].OpenOneParty(com, Alice);
// 				sx[idx].OpenOneParty(com, Bob);
// 			}
// 		}
// 		
// 		if(party == Alice || party == Bob)
// 		{
// 			std::vector<uint64_t> out;
// 			for(int idx = 0; idx < sx.size(); idx++)
// 			{
// 				out.insert(out.end(), sout[idx].begin(), sout[idx].end());
// 			}
// 			
// 			
// 			unsigned char *outarray = new unsigned char[out.size()*sizeof(uint64_t)];
// 			
// 			sse_trans((uint8_t *)(out.data()), outarray, bitlength, nItems);
// 		}
// 		t.Tick("Open psi");
	}
	
	void fPSI_MCA(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		
		int nBytes = 10; int bitlength = rate; nBytes = bitlength/8; rate = 0.5;
		
		unsigned char *bytes = new unsigned char[nBytes*nItems];
		
		
		std::vector<Shares> sx(bitlength), sy(bitlength);
		
		std::vector<uint64_t> data(nItems*nBytes/8);
		
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems*nBytes/8, rate);
			
			if(party == Bob)
			{
				std::reverse(data.begin(), data.end());
			}
			
			// Transpose bits
			sse_trans((uint8_t *)(data.data()), bytes, nItems, bitlength);
			
			int numBlks = nItems/64;
			
			uint64_t *ptr = (uint64_t *)bytes;
			
			std::cout << "Secret share" << std::endl;
			for(int idx = 0; idx < bitlength; idx++)
			{
				std::vector<uint64_t> vec(ptr, ptr + numBlks);
				
				sx[idx].SecretShare(vec, nItems, com, Alice);
				sy[idx].SecretShare(vec, nItems, com, Bob);
				
				ptr += numBlks;
			}
		}
		else
		{
			std::cout << "Secret share" << std::endl;
			for(int idx = 0; idx < bitlength; idx++)
			{
				std::vector<uint64_t> vec(nItems/64);
				sx[idx].SecretShare(vec, nItems, com, Alice);
				sy[idx].SecretShare(vec, nItems, com, Bob);
			}
		}
		
		// Verify inputs
		Triplets tps;
		std::vector<Shares> sx2(bitlength), sy2(bitlength);
		for(int idx = 0; idx < bitlength; idx++)
		{
			sx2[idx] = sx[idx].LeftShiftOne();
			sy2[idx] = sy[idx].RightShiftOne();
		}
		
		if(tps.IsIncreasing(sx2, sx, com)) std::cout << "Passing input 1 check" << std::endl;
		if(tps.IsIncreasing(sy2, sy, com)) std::cout << "Passing input 2 check" << std::endl;
		
// 		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
// 		tpsInput.VerifyTriples(rngseed, com);
// 		
		t.Tick("Secret share inputs");
		
// 		Triplets tps;
		
		std::cout << "size: " << (1 + log2(nItems)) << std::endl;
		// CompareAndSwap
		for(int idx = 0; idx < (1 + log2(nItems)); idx++)
		{
			tps.CompareAndSwap(sx, sy, com);
		}
		
		std::vector<Shares> sx1 = sx;
		
// 		tps.VerifyTriples(rngseed, com);
		
		t.Tick("CompareAndSwap");
		
		// Pairwise Comparison
		// Merge sx and sy
		for(int idx = 0; idx < sx.size(); idx++)
		{
			sx[idx].s.insert(sx[idx].s.end(), sy[idx].s.begin(), sy[idx].s.end());
			sx[idx].t.insert(sx[idx].t.end(), sy[idx].t.begin(), sy[idx].t.end());
		}
		
		
		for(int idx = 0; idx < bitlength; idx++)
		{
			sx2[idx] = sx[idx].LeftShiftOne();
// 			sx3[idx] = sx2[idx].LeftShiftOne();
		}
		
		Shares eq, eq2;
		
		tps.PairwiseComparison2(eq, sx, sx2, com);
		std::cout << "eq len: " << eq.s.size() << std::endl;
		t.Tick("PairwiseComparison");
		
		int len = eq.s.size();
		eq2.s.insert(eq2.s.end(), eq.s.begin() + len/2, eq.s.end());
		eq2.t.insert(eq2.t.end(), eq.t.begin() + len/2, eq.t.end());
		
		eq.s.resize(len/2);
		eq.t.resize(len/2);
		
		std::vector<Shares> in1(1), in2(1), out1, out2;
		in1[0] = eq; in2[0] = eq2;
			
		
		while(len > 2)
		{
// 			std::cout << "len: " << len << std::endl;
			out1 = tps.RippleCarryAdder(in1, in2, com);
			out2.resize(out1.size());
			
			len = len/2;
			for(int idx = 0; idx < out1.size(); idx++)
			{
				out2[idx].s.resize(0); out2[idx].t.resize(0);
				out2[idx].s.insert(out2[idx].s.end(), out1[idx].s.begin() + len/2, out1[idx].s.end());
				out2[idx].t.insert(out2[idx].t.end(), out1[idx].t.begin() + len/2, out1[idx].t.end());
				out1[idx].s.resize(len/2);
				out1[idx].t.resize(len/2);
			}
			
			in1 = out1; in2 = out2;
		}
		
		std::cout << "Final addition" << in1.size() << " " << in1[0].s.size() << std::endl;
		out1 = tps.RippleCarryAdder(in1, in2, com);
// 		std::cout << "Final addition" << std::endl;
		out2 = tps.RippleCarryAdder(in1, in2, com);
		
		t.Tick("Addition Circuit");
		std::vector<unsigned char> seed = com->getCommonSeed();
		tps.VerifyTriples(seed, com);
		
		t.Tick("Verify Triplets");
		
		std::cout << "Comm fpsi: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		std::vector<std::vector<uint64_t> > sout(sx.size());
// 		
// 		for(int idx = 0; idx < sx.size(); idx++)
// 		{
// 			if(party == Alice)
// 			{
// 				sout[idx] = sx[idx].OpenOneParty(com, Alice);
// 				sx[idx].OpenOneParty(com, Bob);
// 			}
// 			else if(party == Bob)
// 			{
// 				sx[idx].OpenOneParty(com, Alice);
// 				sout[idx] = sx[idx].OpenOneParty(com, Bob);
// 			}
// 			else
// 			{
// 				sx[idx].OpenOneParty(com, Alice);
// 				sx[idx].OpenOneParty(com, Bob);
// 			}
// 		}
// 		
// 		if(party == Alice || party == Bob)
// 		{
// 			std::vector<uint64_t> out;
// 			for(int idx = 0; idx < sx.size(); idx++)
// 			{
// 				out.insert(out.end(), sout[idx].begin(), sout[idx].end());
// 			}
// 			
// 			
// 			unsigned char *outarray = new unsigned char[out.size()*sizeof(uint64_t)];
// 			
// 			sse_trans((uint8_t *)(out.data()), outarray, bitlength, nItems);
// 		}
// 		t.Tick("Open psi");
	}
	
// 	void ComputeMCSCircuit(int party, Triplets& tps, const std::vector<Shares>& sx, const std::vector<Shares>& sy, int nItems, const std::string circuitFileName, Communicator *com)
// 	{
// 		// Load AES circuit 
// 		int nWires = 0;
// 		int nGates = 0;
// 		
// 		// Read file
// 		std::ifstream sr(circuitFileName);
// 		std::string line;
// 		
// 		std::getline(sr, line);
// 		
// 		auto s1 = StringHelper::split(line, ' ');
// 		nGates = std::stoi(s1[0]);
// 		nWires = std::stoi(s1[1]);
// 		
// 		// Read input/output size (we already know: 128)
// 		std::getline(sr, line);
// 		auto s2 = StringHelper::split(line, ' ');
// 		int inputSize  = std::stoi(s2[0]);
// 		int outputSize = std::stoi(s2[2]);
// 		
// 		std::getline(sr, line); // empty line
// 	  
// 		// Prepare input
// 		std::vector<Wires> wires(nWires);
// 		
// 		for(int idx = 0; idx < inputSize; idx++)
// 		{
// 			wires[idx].val = sx[idx];
// 			wires[idx + inputSize] = sy[idx];
// 		}
// 		
// 		std::cout << "nWires: " << nWires << " nGates: " << nGates << " input length: " << inputSize << std::endl;
// 		Shares zero;
// 		zero.s.resize(input[0].s.size());
// 		zero.t.resize(input[0].s.size());
// 		
// 		for(int idx = 0; idx < zero.s.size(); idx++)
// 		{
// 			zero.s[idx] = 0; zero.t[idx] = 0;
// 		}
// 		
// 		for(int idx = 0; idx < nWires; idx++) wires[idx].layer = 0;
// 		
// 		int length = wires[0].val.s.size();
// 		
// 		// Layering the circuit
// 		std::vector<std::vector<IGate *> > Layers;
// 		
// 		// Init the first layer, which is empty.
// 		Layers.resize(1); Layers[0].resize(0);
// 		
// 		std::cout << "Layering the circuit" << std::endl;
// 		// Add gates to current layers and append new layers if needed
// 		for (uint64_t idx = 0; idx < nGates; idx++)
// 		{
// 			std::getline(sr, line);
// 			
// 			auto variables = StringHelper::split(line, ' ');
// 			
// 			uint64_t nInputWires = std::stoi(variables[0]);
// 			
// 			if(nInputWires == 1)
// 			{
// 				uint64_t in = std::stoi(variables[2]);
// 				uint64_t out = std::stoi(variables[3]);
// 				
// 				wires[out].layer = wires[in].layer;
// 				IGate *gate = new UnitaryGate(in, out);
// 				Layers[wires[out].layer].push_back(gate);
// 			}
// 			else if(nInputWires == 2)
// 			{
// 				uint64_t in1 = std::stoi(variables[2]);
// 				uint64_t in2 = std::stoi(variables[3]);
// 				uint64_t out = std::stoi(variables[4]);
// 				
// 				IGate *gate = new BinaryGate(in1, in2, out);
// 				
// 				uint64_t maxLayer = (wires[in1].layer >= wires[in2].layer) ? (wires[in1].layer) : (wires[in2].layer);
// 				
// 				if(variables[5][0] == 'A')
// 				{
// 					wires[out].layer = maxLayer + 1;
// 					if(wires[out].layer >= Layers.size())
// 					{
// 						std::vector<IGate *> temp;
// 						Layers.push_back(temp);
// 					}
// 				}
// 				else
// 				{
// 					wires[out].layer = maxLayer;
// 				}
// 				
// 				Layers[maxLayer].push_back(gate);
// 			}
// 			else
// 			{
// 				assert(0);
// 			}
// 		}
// 		
// 		std::cout << "Executing the circuit: " << Layers.size() << " layers" << std::endl;
// 		
// 		// Execute circuit by layers
// 		for(uint64_t idx = 0; idx < Layers.size(); idx++)
// 		{
// 			std::cout << "Execute layer: " << idx << " " << Layers[idx].size() << " gates" << std::endl;
// 			std::vector<Shares> left, right;
// 			std::vector<uint64_t> outIdx;
// 			
// 			for(uint64_t gdx = 0; gdx < Layers[idx].size(); gdx++)
// 			{
// 				if (dynamic_cast<UnitaryGate *>(Layers[idx][gdx]) != nullptr)
// 				{
// // 					std::cout << gdx << ": NOT" << std::endl;
// 					UnitaryGate *gate = (UnitaryGate *)(Layers[idx][gdx]);
// 				
// 					wires[gate->OutputWire].val = Shares::ComputeNOT(wires[gate->InputWire].val);
// 				}
// 				else 
// 				{
// 					BinaryGate *gate = (BinaryGate *)(Layers[idx][gdx]);
// 					if((wires[gate->OutputWire].layer == wires[gate->LeftWire].layer) || (wires[gate->OutputWire].layer == wires[gate->RightWire].layer))
// 					{
// // 						std::cout << gdx << ": XOR" << std::endl;
// 						wires[gate->OutputWire].val = Shares::ComputeXOR(wires[gate->LeftWire].val, wires[gate->RightWire].val);
// 					}
// 					else
// 					{
// // 						std::cout << gdx << ": AND" << std::endl;
// 						left.push_back(wires[gate->LeftWire].val);
// 						right.push_back(wires[gate->RightWire].val);
// 						outIdx.push_back(gate->OutputWire);
// 					}
// 				}
// 			}
// 				
// 			if(outIdx.size() > 0)
// 			{
// 				std::cout << "Execute AND" << std::endl;
// 				std::vector<Shares> out = tps.ComputeAND(left, right, com);
// 				for(uint64_t idx = 0; idx < out.size(); idx++) wires[outIdx[idx]].val = out[idx];
// 			}
// 		}
// 		
// 		output.resize(outputSize);
// 		
// 		uint64_t start = nWires - outputSize - 1;
// 		for(uint64_t idx = 0; idx < outputSize; idx++)
// 		{
// 			output[idx] = wires[start + idx].val;
// 		}
// 	}
// 	
	// Compute a function over the intersection
	void fPSI(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		
		int fieldSize = 64; int nBytes = 10; int bitlength = 64;
		
		unsigned char *bytes = new unsigned char[nBytes*nItems];
		
		ZZ_pContext context;
		context.save();
		
		ZZ prime64bit = (ZZ(1) << 64) - 59;
		ZZ_p::init(prime64bit);
		nBytes = 8;
		
		// Secret share the data
		MPCOperations mpc;
		RepArithShares MacKeyData, MacKeyPRF, sx, sy, sz1, sz2; 
		
		std::vector<ZZ_p> Fx(nItems);
		
		if(party != Charlie)
		{
			std::vector<uint64_t> data;
			GenerateInput(party, data, nItems, rate);
			
			for(int idx = 0; idx < nItems; idx++)
			{
				Fx[idx] = (data[idx]);
			}
		}
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Sharing input");
		
		mpc.MPC_Shuffle(sx, nBytes, com);
		
		t.Tick("3pc Shuffling");
		
		// Convert shares in Z2 to binary shares
		int size = sx.x.size();
		
		unsigned char *inX = new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size];
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
// 		// Compact the bytes
// 		unsigned char *inX = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.x.data()); //new unsigned char[nBytes*size];
// 		unsigned char *inY = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.y.data()); //new unsigned char[nBytes*size];
// 		
// 		unsigned char *ptr = 0;
// 		int count = 0;
// 		
// 		ptr = (unsigned char *)(sx.x.data());
// 		for(int idx = 0; idx < sx.x.size(); idx++)
// 		{
// 			for(int kdx = 0; kdx < nBytes; kdx++)
// 			{
// 				inX[count] = *ptr;
// 				count++;
// 				ptr++;
// 			}
// 			ptr += (16 - nBytes);
// 		}
// 		
// 		count = 0;
// 		ptr = (unsigned char *)(sx.y.data());
// 		for(int idx = 0; idx < sx.y.size(); idx++)
// 		{
// 			for(int kdx = 0; kdx < nBytes; kdx++)
// 			{
// 				inY[count] = *ptr;
// 				count++;
// 				ptr++;
// 			}
// 			ptr += (16 - nBytes);
// 		}
// 		
// 		t.Tick("Compact bytes");
		
		std::vector<Shares> bx1(bitlength);
		std::vector<Shares> bx2(bitlength);
		std::vector<Shares> bx3(bitlength);
		
		ArrayEncoder::ZZArray2ByteArray(inX, sx.x.data(), nBytes, size);
		ArrayEncoder::ZZArray2ByteArray(inY, sx.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, bitlength);
		sse_trans((uint8_t *)inY, outY, size, bitlength);
		
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < bitlength; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
		sum.resize(fieldSize);
		// sum = sum mod p
// 		tps.modp_fast(sum, prime64bit, com);
		
		std::vector<unsigned char> seed = com->getCommonSeed();
		tps.VerifyTriples(seed, com);
		
		context.restore();
		
		delete [] bytes;
		
		t.Tick("Convert [x]_A to [x]_B");
		
		Triplets tpsAES;
		
		std::vector<Shares> output;
		ComputeAES_fast(party, tpsAES, output, sum, nItems, nBytes, com);
		t.Tick("AES");
		
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tpsAES.VerifyTriples(rngseed, com);
		
		std::vector< std::vector<uint64_t> > rOut(output.size());
		
		for(int idx = 0; idx < output.size(); idx++)
		{
			rOut[idx] = output[idx].OpenAllParty(com);
		}
		
		t.Tick("Verification");
	}
	
	
	
	// Execute multiple circuits at once
	void ComputeAES_fast(int party, Triplets& tps, std::vector<Shares>& output, std::vector<Shares>& input, int nItems, int nBytes, Communicator *com)
	{
		// Load AES circuit 
		int nWires = 0;
		int nGates = 0;
		
		const int inputSize = 128;
		const int outputSize = 128;
		const int keySize = 128;
		
		std::string filePath("./AES-non-expanded.txt");
		
		// Read file
		std::ifstream sr(filePath);
		std::string line;
		
		std::getline(sr, line);
		
		auto sizes = StringHelper::split(line, ' ');
		nGates = std::stoi(sizes[0]);
		nWires = std::stoi(sizes[1]);
		
		// Read input/output size (we already know: 128)
		std::getline(sr, line);
		
		std::getline(sr, line); // empty line
	  
		// Prepare input
		std::vector<Wires> wires(nWires);
		
		for(int idx = 0; idx < input.size(); idx++)
		{
			wires[idx].val = input[idx];
		}
		
		Shares zero;
		zero.s.resize(input[0].s.size());
		zero.t.resize(input[0].s.size());
		
		for(int idx = 0; idx < zero.s.size(); idx++)
		{
			zero.s[idx] = 0; zero.t[idx] = 0;
		}
		
		for(int idx = input.size(); idx < inputSize; idx++)
		{
			wires[idx].val = zero;
		}
		
		for(int idx = 0; idx < nWires; idx++) wires[idx].layer = 0;
		
		// Prepare key 
		Shares key;
		key.GenerateRandomShares(keySize, com);
		
		std::string s = std::bitset<64>(key.s[0]).to_string() + std::bitset<64>(key.s[1]).to_string();
		std::string t = std::bitset<64>(key.t[0]).to_string() + std::bitset<64>(key.t[1]).to_string();
		
		std::vector<Shares> keys(keySize);
		
		int length = wires[0].val.s.size();
		
		for(int idx = 0; idx < keySize; idx++)
		{
			keys[idx].s.resize(length);
			keys[idx].t.resize(length);
			
			if(s[idx] == '0')
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].s[kdx] = 0;
				}
			}
			else
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].s[kdx] = 0xFFFFFFFFFFFFFFFFULL;
				}
			}
			
			if(t[idx] == '0')
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].t[kdx] = 0;
				}
			}
			else
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].t[kdx] = 0xFFFFFFFFFFFFFFFFULL;
				}
			}
		}
		
		for(int idx = inputSize; idx < inputSize + keySize; idx++)
		{
			wires[idx].isConstant = false;
			wires[idx].val = keys[idx - inputSize];
		}
		
		// Layering the circuit
		std::vector<std::vector<IGate *> > Layers;
		
		// Init the first layer, which is empty.
		Layers.resize(1); Layers[0].resize(0);
		
		std::cout << "Layering the circuit" << std::endl;
		// Add gates to current layers and append new layers if needed
		for (uint64_t idx = 0; idx < nGates; idx++)
		{
			std::getline(sr, line);
			
			auto variables = StringHelper::split(line, ' ');
			
			uint64_t nInputWires = std::stoi(variables[0]);
			
			if(nInputWires == 1)
			{
				uint64_t in = std::stoi(variables[2]);
				uint64_t out = std::stoi(variables[3]);
				
				wires[out].layer = wires[in].layer;
				IGate *gate = new UnitaryGate(in, out);
				Layers[wires[out].layer].push_back(gate);
			}
			else if(nInputWires == 2)
			{
				uint64_t in1 = std::stoi(variables[2]);
				uint64_t in2 = std::stoi(variables[3]);
				uint64_t out = std::stoi(variables[4]);
				
				IGate *gate = new BinaryGate(in1, in2, out);
				
				uint64_t maxLayer = (wires[in1].layer >= wires[in2].layer) ? (wires[in1].layer) : (wires[in2].layer);
				
				if(variables[5][0] == 'A')
				{
					wires[out].layer = maxLayer + 1;
					if(wires[out].layer >= Layers.size())
					{
						std::vector<IGate *> temp;
						Layers.push_back(temp);
					}
				}
				else
				{
					wires[out].layer = maxLayer;
				}
				
				Layers[maxLayer].push_back(gate);
			}
			else
			{
				assert(0);
			}
		}
		
		std::cout << "Executing the circuit" << std::endl;
		
		// Execute circuit by layers
		for(uint64_t idx = 0; idx < Layers.size(); idx++)
		{
// 			std::cout << "Execute layer: " << idx << " " << Layers[idx].size() << " gates" << std::endl;
			std::vector<Shares> left, right;
			std::vector<uint64_t> outIdx;
			
			for(uint64_t gdx = 0; gdx < Layers[idx].size(); gdx++)
			{
				if (dynamic_cast<UnitaryGate *>(Layers[idx][gdx]) != nullptr)
				{
					UnitaryGate *gate = (UnitaryGate *)(Layers[idx][gdx]);
				
					wires[gate->OutputWire].val = Shares::ComputeNOT(wires[gate->InputWire].val);
				}
				else 
				{
					BinaryGate *gate = (BinaryGate *)(Layers[idx][gdx]);
					if((wires[gate->OutputWire].layer == wires[gate->LeftWire].layer) || (wires[gate->OutputWire].layer == wires[gate->RightWire].layer))
					{
						wires[gate->OutputWire].val = Shares::ComputeXOR(wires[gate->LeftWire].val, wires[gate->RightWire].val);
					}
					else
					{
						left.push_back(wires[gate->LeftWire].val);
						right.push_back(wires[gate->RightWire].val);
						outIdx.push_back(gate->OutputWire);
					}
				}
			}
				
			if(outIdx.size() > 0)
			{
				std::vector<Shares> out = tps.ComputeAND(left, right, com);
				for(uint64_t idx = 0; idx < out.size(); idx++) wires[outIdx[idx]].val = out[idx];
			}
		}
		
		output.resize(outputSize);
		
		uint64_t start = nWires - outputSize - 1;
		for(uint64_t idx = 0; idx < outputSize; idx++)
		{
			output[idx] = wires[start + idx].val;
		}
	}
	
	void ComputeAES(int party, Triplets& tps, std::vector<Shares>& output, std::vector<Shares>& input, int nItems, int nBytes, Communicator *com)
	{
		// Load AES circuit 
		int nWires = 0;
		int nGates = 0;
		
		const int inputSize = 128;
		const int outputSize = 128;
		const int keySize = 128;
		
		std::string filePath("./AES-non-expanded.txt");
		
		// Read file
		std::ifstream sr(filePath);
		std::string line;
		
		std::getline(sr, line);
		
		auto sizes = StringHelper::split(line, ' ');
		nGates = std::stoi(sizes[0]);
		nWires = std::stoi(sizes[1]);
		
		// Read input/output size (we already know: 128)
		std::getline(sr, line);
		
		std::getline(sr, line); // empty line
	  
		// Prepare input
		std::vector<Wires> wires(nWires);
		
		for(int idx = 0; idx < input.size(); idx++)
		{
			wires[idx].val = input[idx];
		}
		
		Shares zero;
		zero.s.resize(input[0].s.size());
		zero.t.resize(input[0].s.size());
		
		for(int idx = 0; idx < zero.s.size(); idx++)
		{
			zero.s[idx] = 0; zero.t[idx] = 0;
		}
		
		for(int idx = input.size(); idx < inputSize; idx++)
		{
			wires[idx].val = zero;
		}
		
		// Prepare key 
		Shares key;
		key.GenerateRandomShares(keySize, com);
		
		std::string s = std::bitset<64>(key.s[0]).to_string() + std::bitset<64>(key.s[1]).to_string();
		std::string t = std::bitset<64>(key.t[0]).to_string() + std::bitset<64>(key.t[1]).to_string();
		
		std::vector<Shares> keys(keySize);
		
		int length = wires[0].val.s.size();
		
		for(int idx = 0; idx < keySize; idx++)
		{
			keys[idx].s.resize(length);
			keys[idx].t.resize(length);
			
			if(s[idx] == '0')
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].s[kdx] = 0;
				}
			}
			else
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].s[kdx] = 0xFFFFFFFFFFFFFFFFULL;
				}
			}
			
			if(t[idx] == '0')
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].t[kdx] = 0;
				}
			}
			else
			{
				for(int kdx = 0; kdx < length; kdx++)
				{
					keys[idx].t[kdx] = 0xFFFFFFFFFFFFFFFFULL;
				}
			}
		}
		
		for(int idx = inputSize; idx < inputSize + keySize; idx++)
		{
			wires[idx].isConstant = false;
			wires[idx].val = keys[idx - inputSize];
		}
		
		// Execute the circuit
		for (int idx = 0; idx < nGates; idx++)
		{
			std::getline(sr, line);
			
			auto variables = StringHelper::split(line, ' ');
			
			int nInputWires = std::stoi(variables[0]);
			
			if(nInputWires == 1)
			{
				int in = std::stoi(variables[2]);
				int out = std::stoi(variables[3]);
				
				wires[out].val = Shares::ComputeNOT(wires[in].val);
			}
			else if(nInputWires == 2)
			{
				int in1 = std::stoi(variables[2]);
				int in2 = std::stoi(variables[3]);
				int out = std::stoi(variables[4]);
				
				if(variables[5][0] == 'A')
				{
					wires[out].val = tps.ComputeAND(wires[in1].val, wires[in2].val, com);
				}
				else
				{
					wires[out].val = Shares::ComputeXOR(wires[in1].val, wires[in2].val);
				}
			}
			else
			{
				assert(0);
			}
		}
		
		output.resize(outputSize);
		
		int start = nWires - outputSize - 1;
		for(int idx = 0; idx < outputSize; idx++)
		{
			output[idx] = wires[start + idx].val;
		}
	}
	
	void ComputeSOPRF(int party, std::vector<ZZ_p>& Fx, std::vector<Shares>& sum, int nItems, int nBytes, int bitlength, Communicator *com)
	{
		unsigned char *bytes = new unsigned char[nBytes*nItems];
		
		ZZ_pContext context;
		context.save();
		// Distribute (shared/whole) prf keys
		// 1024 bit prime number
		ZZ p = conv<ZZ>("124325339146889384540494091085456630009856882741872806181731279018491820800119460022367403769795008250021191767583423221479185609066059226301250167164084041279837566626881119772675984258163062926954046545485368458404445166682380071370274810671501916789361956272226105723317679562001235501455748016154805420913");
		
		ZZ q = conv<ZZ>("1399252811935680595399801714158014275474696840019");
		ZZ g = conv<ZZ>("115740200527109164239523414760926155534485715860090261532154107313946218459149402375178179458041461723723231563839316251515439564315555249353831328479173170684416728715378198172203100328308536292821245983596065287318698169565702979765910089654821728828592422299160041156491980943427556153020487552135890973413");
		
		
		ZZ_p::init(q); nBytes = 20;
		
		// Distribute keys
		int nKeys = bitlength + 1;
		std::vector<ZZ_p> prfKeys(nKeys);
		std::vector<ZZ_p> sharedPrfKeys(nKeys);
		
		if(party == Alice)
		{
			for(int idx = 0; idx < nKeys; idx++)
			{
				random(prfKeys[idx]);
				random(sharedPrfKeys[idx]);
			}
			
			prfKeys[0] = conv<ZZ_p>(1);
			sharedPrfKeys[0] = conv<ZZ_p>(1);
		
			ArrayEncoder::ZZArray2ByteArray(bytes, prfKeys.data(), nBytes, nKeys);
			com->SendBob(bytes, nBytes*nKeys);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, sharedPrfKeys.data(), nBytes, nKeys);
			com->SendCharlie(bytes, nBytes*nKeys);
			
			for(int idx = 0; idx < nKeys; idx++)
			{
				mul(prfKeys[idx], prfKeys[idx], sharedPrfKeys[idx]);
			}
		}
		else if(party == Bob)
		{
			com->AwaitAlice(bytes, nBytes*nKeys);
			ArrayEncoder::ByteArray2ZZArray(sharedPrfKeys.data(), bytes, nBytes, nKeys);
		}
		else
		{
			com->AwaitAlice(bytes, nBytes*nKeys);
			ArrayEncoder::ByteArray2ZZArray(sharedPrfKeys.data(), bytes, nBytes, nKeys);
		}
		
		// Compute PRF and soPRF
		std::vector<ZZ_p> prfs(nItems);
		std::vector<ZZ_p> sprfs(2*nItems);
		
		if(party == Alice)
		{
			for(int idx = 0; idx < nItems; idx++)
			{
				ZZ_p exp = conv<ZZ_p>(1); 
				for(int kdx = 0; kdx < bitlength; kdx++)
				{
					if(IsOne((conv<ZZ>(Fx[idx]) & (ZZ(1) << kdx)) >> kdx))
					{
						mul(exp, exp, prfKeys[kdx + 1]);
					}
				}
				if(idx == 0) std::cout << "exp: " << exp << std::endl;
				
				mul(exp, exp, prfKeys[0]);
				
				if(idx == 0) std::cout << "g/beta: " << g << std::endl;
				
				context.save(); 
				ZZ_p::init(p);
				
				power(prfs[idx], conv<ZZ_p>(g), conv<ZZ>(exp));
				
				if((idx % 10000) == 0) std::cout << "prf[" << idx << "]: " << prfs[idx] << std::endl;
				context.restore();
			}
		}
		else
		{
// 			NetIO * io = new NetIO(party==Alice ? nullptr:"127.0.0.1", 40000);
			std::vector<ZZ_p> ab(bitlength);
			
			std::vector<ZZ_p> left(bitlength), right(bitlength);
			std::vector<ZZ_p> selectedOT(bitlength);
			
			unsigned char *larray = new unsigned char[bitlength*nBytes];
			unsigned char *rarray = new unsigned char[bitlength*nBytes];
			unsigned char *bytes = new unsigned char[128]; // To hold beta
			
			for(int idx = 0; idx < nItems; idx++)
			{
				for(int kdx = 0; kdx < bitlength; kdx++)
				{
					random(ab[kdx]);
				}
				
				if(party == Charlie)
				{
					for(int kdx = 0; kdx < bitlength; kdx++)
					{
						// Bob prepares two array
						// Get the idx-th bit of sum[kdx]
						int r = idx / bitlength;
						int c = idx % bitlength;
						
						bool b = ((sum[kdx].s[r] >> c) & 1);
// 						if(idx == 0) std::cout << kdx << " " << b << std::endl;
						if(b)
						{
							left[kdx] = ab[kdx]*sharedPrfKeys[kdx + 1];
							right[kdx] = ab[kdx];
						}
						else
						{
							left[kdx] = ab[kdx];
							right[kdx] = ab[kdx]*sharedPrfKeys[kdx + 1];
						}
					}
					
					// Performing OT with Charlie
					ArrayEncoder::ZZArray2ByteArray(larray, left.data(), nBytes, bitlength);
					ArrayEncoder::ZZArray2ByteArray(rarray, right.data(), nBytes, bitlength);
					
					com->SendBob(larray, nBytes*bitlength);
					com->SendBob(rarray, nBytes*bitlength);
					
					// Second round OT
					com->AwaitBob(larray, nBytes*bitlength);
					com->AwaitBob(rarray, nBytes*bitlength);
					
					ArrayEncoder::ByteArray2ZZArray(left.data(), larray, nBytes, bitlength);
					ArrayEncoder::ByteArray2ZZArray(right.data(), rarray, nBytes, bitlength);
					
					for(int kdx = 0; kdx < bitlength; kdx++)
					{
						// Bob prepares two array
						// Get the idx-th bit of sum[kdx]
						int r = idx / bitlength;
						int c = idx % bitlength;
						
						bool b = ((sum[kdx].s[r] >> c) & 1);
						if(b)
						{
							selectedOT[kdx] = right[kdx];
						}
						else
						{
							selectedOT[kdx] = left[kdx];
						}
					}
					
					// Compute alpha
					ZZ_p v = conv<ZZ_p>(1), a = conv<ZZ_p>(1);
					for(int idx = 0; idx < bitlength; idx++)
					{
						mul(v, v, selectedOT[idx]);
						mul(a, a, ab[idx]);
					}
					
					ZZ_p alpha = sharedPrfKeys[0]*v*inv(a);
					
					context.save();
					ZZ_p::init(p); nBytes = 128;
					
					ZZ temp;
					
					com->AwaitBob(bytes, nBytes);
					ZZFromBytes(temp, bytes, nBytes);
					
					power(prfs[idx], conv<ZZ_p>(temp), conv<ZZ>(alpha));
					
					if((idx % 10000) == 0) std::cout << "prf[" << idx << "]: " << prfs[idx] << std::endl;
					context.restore(); nBytes = 20;
				}
				else
				{
					com->AwaitCharlie(larray, nBytes*bitlength);
					com->AwaitCharlie(rarray, nBytes*bitlength);
					
					ArrayEncoder::ByteArray2ZZArray(left.data(), larray, nBytes, bitlength);
					ArrayEncoder::ByteArray2ZZArray(right.data(), rarray, nBytes, bitlength);
					
					for(int kdx = 0; kdx < bitlength; kdx++)
					{
						// Bob prepares two array
						// Get the idx-th bit of sum[kdx]
						int r = idx / 64;
						int c = idx % 64;
						
						bool b = ((sum[kdx].t[r] >> c) & 1);
// 						if(idx == 0) std::cout << kdx << " " << b << std::endl;
						if(b)
						{
							selectedOT[kdx] = right[kdx];
						}
						else
						{
							selectedOT[kdx] = left[kdx];
						}
					}
					
					// Second round OT
					for(int kdx = 0; kdx < bitlength; kdx++)
					{
						// Bob prepares two array
						// Get the idx-th bit of sum[kdx]
						int r = idx / 64;
						int c = idx % 64;
						
						bool b = ((sum[kdx].t[r] >> c) & 1);
						if(b)
						{
							left[kdx] = ab[kdx]*sharedPrfKeys[kdx + 1]*selectedOT[kdx];
							right[kdx] = ab[kdx]*selectedOT[kdx];
						}
						else
						{
							left[kdx] = ab[kdx]*selectedOT[kdx];
							right[kdx] = ab[kdx]*sharedPrfKeys[kdx + 1]*selectedOT[kdx];
						}
					}
					
					// Performing OT with Charlie
					ArrayEncoder::ZZArray2ByteArray(larray, left.data(), nBytes, bitlength);
					ArrayEncoder::ZZArray2ByteArray(rarray, right.data(), nBytes, bitlength);
					
					com->SendCharlie(larray, nBytes*bitlength);
					com->SendCharlie(rarray, nBytes*bitlength);
					
					// Compute beta
					ZZ_p b = conv<ZZ_p>(1);
					for(int idx = 0; idx < bitlength; idx++)
					{
						mul(b, b, ab[idx]);
					}
					
					ZZ_p beta = sharedPrfKeys[0]*inv(b);
					
					context.save();
					
					ZZ_p::init(p); nBytes = 128;
					
					ZZ_p temp;
					
					power(temp, conv<ZZ_p>(g), rep(beta));
					
					BytesFromZZ(bytes, rep(temp), nBytes);
					
					com->SendCharlie(bytes, nBytes);
					
					context.restore(); nBytes = 20;
				}
			}
		}
		
		delete [] bytes;
	}
	
	void PSI_Circuit(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
		
		MPCRingOperations mpc;
		
		RepRingShares sx, sy, sz1, sz2;
		
		std::vector<__uint128_t> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
		
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
			
			rng->ComputeRingFx(Fx, data);
		}
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
		std::cout << "Shuffling shares" << std::endl;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
		t.Tick("Sorting");
		
		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepRingShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = (mpc.modulo + sz1.x[2*idx] - sz1.x[2*idx + 1]);
			psi.y[idx] = (mpc.modulo + sz1.y[2*idx] - sz1.y[2*idx + 1]);
			mpc.modp(psi.x[idx]);
			mpc.modp(psi.y[idx]);
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		nBytes = 10;
		if(sz1.x.size()/2 == nItems)
		{
			std::vector<__uint128_t> psiOutput;
		
			if(party == Alice)
			{
				psiOutput = mpc.openOneParty(psi, nBytes, com, Alice);
				mpc.openOneParty(psi, nBytes, com, Bob);
			}
			else if(party == Bob)
			{
				mpc.openOneParty(psi, nBytes, com, Alice);
				psiOutput = mpc.openOneParty(psi, nBytes, com, Bob);
			}
			else
			{
				mpc.openOneParty(psi, nBytes, com, Alice);
				mpc.openOneParty(psi, nBytes, com, Bob);
			}
			
			return;
		}
		
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		nBytes = 10;
		
		// Compact the bytes
		unsigned char *inX = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.x.data()); //new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.y.data()); //new unsigned char[nBytes*size];
		
		unsigned char *ptr = 0;
		int count = 0;
		
		ptr = (unsigned char *)(sz2.x.data());
		for(int idx = 0; idx < sz2.x.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inX[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		count = 0;
		ptr = (unsigned char *)(sz2.y.data());
		for(int idx = 0; idx < sz2.y.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inY[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		t.Tick("Compact bytes");
		
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
		std::vector<Shares> bx1(fieldSize);
		std::vector<Shares> bx2(fieldSize);
		std::vector<Shares> bx3(fieldSize);
		
// 		ArrayEncoder::ZZArray2ByteArray(inX, sz2.x.data(), nBytes, size);
// 		ArrayEncoder::ZZArray2ByteArray(inY, sz2.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, fieldSize);
		sse_trans((uint8_t *)inY, outY, size, fieldSize);
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < fieldSize; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
// 		sum.resize(fieldSize);
// 		sum.resize(sum.size() - 2);
		// sum = sum mod p
		tps.modp_fast(sum, p, com);
		
		std::vector<Shares> sum2(fieldSize);
		for(int idx = 0; idx < fieldSize; idx++)
		{
			sum2[idx] = sum[idx].LeftShiftOne();
		}
		
		if(tps.IsIncreasing(sum2, sum, com))
		{
			std::cout << "Passed verification" << std::endl;
		}
		else
		{
			std::cout << "Failed verification" << std::endl;
		}
		
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tps.VerifyTriples(rngseed, com);
		t.Tick("PSU proof");
		
		std::vector<__uint128_t> psiOutput;
		
		if(party == Alice)
		{
			psiOutput = mpc.openOneParty(psi, nBytes, com, Alice);
			mpc.openOneParty(psi, nBytes, com, Bob);
		}
		else if(party == Bob)
		{
			mpc.openOneParty(psi, nBytes, com, Alice);
			psiOutput = mpc.openOneParty(psi, nBytes, com, Bob);
		}
		else
		{
			mpc.openOneParty(psi, nBytes, com, Alice);
			mpc.openOneParty(psi, nBytes, com, Bob);
		}
		
		t.Tick("Open PSI");
	}
	
	
	void PSI_Circuit2(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
// 		if(nItems > 1000000)
// 		{
// 			ZZ pp = (ZZ(1) << 96) - 17;
// 			ZZ_p::init(pp);
// 			fieldSize = 96; nBytes = 12;
// 		}
		
		MPCOperations mpc;
		
		RepArithShares sx, sy, sz1, sz2;
		
		std::vector<ZZ_p> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
		
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
			
			rng->ComputeFx(Fx, data);
		}
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
		std::cout << "Shuffling shares" << std::endl;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
		t.Tick("Sorting");
		
		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepArithShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = sz1.x[2*idx] - sz1.x[2*idx + 1];
			psi.y[idx] = sz1.y[2*idx] - sz1.y[2*idx + 1];
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		unsigned char *inX = new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size];
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
		std::vector<Shares> bx1(fieldSize);
		std::vector<Shares> bx2(fieldSize);
		std::vector<Shares> bx3(fieldSize);
		
		ArrayEncoder::ZZArray2ByteArray(inX, sz2.x.data(), nBytes, size);
		ArrayEncoder::ZZArray2ByteArray(inY, sz2.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, fieldSize);
		sse_trans((uint8_t *)inY, outY, size, fieldSize);
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < fieldSize; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
// 		sum.resize(sum.size() - 2);
		// sum = sum mod p
		tps.modp_fast(sum, p, com);
		
		std::vector<Shares> sum2(fieldSize);
		for(int idx = 0; idx < fieldSize; idx++)
		{
			sum2[idx] = sum[idx].LeftShiftOne();
		}
		
		if(tps.IsIncreasing(sum2, sum, com))
		{
			std::cout << "Passed verification" << std::endl;
		}
		else
		{
			std::cout << "Failed verification" << std::endl;
		}
		
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tps.VerifyTriples(rngseed, com);
		t.Tick("PSU proof");
		
		std::vector<ZZ_p> psiOutput;
		
		if(party == Alice)
		{
			psiOutput = mpc.openOneParty(psi, nBytes, com, Alice);
			mpc.openOneParty(psi, nBytes, com, Bob);
		}
		else if(party == Bob)
		{
			mpc.openOneParty(psi, nBytes, com, Alice);
			psiOutput = mpc.openOneParty(psi, nBytes, com, Bob);
		}
		else
		{
			mpc.openOneParty(psi, nBytes, com, Alice);
			mpc.openOneParty(psi, nBytes, com, Bob);
		}
		
		t.Tick("Open PSI");
	}
	
	void PSI_CA_Circuit_Payload(int party, int nItems, int pl, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
		
		MPCRingOperations mpc;
		
		RepRingShares sx, sy, sz1, sz2;
		
		std::vector<__uint128_t> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
// 		t.Tick("GenerateInput");
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
// 			t.Tick("rng");
			rng->ComputeRingFx(Fx, data);
		}
		
// 		t.Tick("rng->ComputeRingFx");
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
		RepRingShares sx2 = sx;
		std::vector<__uint128_t> Fx2 = Fx;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
// 		std::cout << "Shuffling shares" << std::endl;
		for(int idx = 0; idx < pl; idx++)
		{
			RepRingShares sx3 = sx2;
			std::vector<__uint128_t> Fx3 = Fx2;
			RepRingShares sz12, sz22;
			mpc.MPC_Sort(sz12, sz22, Fx3, sx3, nBytes, com);
		}
		
		
		t.Tick("Sorting");
		
// 		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepRingShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = (mpc.modulo + sz1.x[2*idx] - sz1.x[2*idx + 1]) % mpc.modulo;
			psi.y[idx] = (mpc.modulo + sz1.y[2*idx] - sz1.y[2*idx + 1]) % mpc.modulo;
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		std::cout << "****************Sending: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		com->sendingSize = 10*com->sendingSize/16;
		
		nBytes = 10;
		
		// Compact the bytes
		unsigned char *inX = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.x.data()); //new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.y.data()); //new unsigned char[nBytes*size];
		
		unsigned char *ptr = 0;
		int count = 0;
		
		ptr = (unsigned char *)(sz2.x.data());
		for(int idx = 0; idx < sz2.x.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inX[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		count = 0;
		ptr = (unsigned char *)(sz2.y.data());
		for(int idx = 0; idx < sz2.y.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inY[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		t.Tick("Compact bytes");
		
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
		std::vector<Shares> bx1(fieldSize);
		std::vector<Shares> bx2(fieldSize);
		std::vector<Shares> bx3(fieldSize);
		
// 		ArrayEncoder::ZZArray2ByteArray(inX, sz2.x.data(), nBytes, size);
// 		ArrayEncoder::ZZArray2ByteArray(inY, sz2.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, fieldSize);
		sse_trans((uint8_t *)inY, outY, size, fieldSize);
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < fieldSize; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		t.Tick("Share conversion");
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
		t.Tick("RippleCarryAdder");
// 		sum.resize(fieldSize);
		// sum = sum mod p
		tps.modp_fast(sum, p, com);
		
		t.Tick("modp_fast");
		
		std::vector<Shares> sum2(fieldSize);
		for(int idx = 0; idx < fieldSize; idx++)
		{
			sum2[idx] = sum[idx].LeftShiftOne();
		}
		t.Tick("Left Shift");
		if(tps.IsIncreasing(sum2, sum, com))
		{
			std::cout << "Passed verification" << std::endl;
		}
		else
		{
			std::cout << "Failed verification" << std::endl;
		}
		t.Tick("Is Increasing");
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tps.VerifyTriples(rngseed, com);
		t.Tick("PSU proof");
	}
	
	
	void PSI_CA_Circuit(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
		
		MPCRingOperations mpc;
		
		RepRingShares sx, sy, sz1, sz2;
		
		std::vector<__uint128_t> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
// 		t.Tick("GenerateInput");
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
// 			t.Tick("rng");
			rng->ComputeRingFx(Fx, data);
		}
		
// 		t.Tick("rng->ComputeRingFx");
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
// 		std::cout << "Shuffling shares" << std::endl;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
		t.Tick("Sorting");
		
// 		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepRingShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = (mpc.modulo + sz1.x[2*idx] - sz1.x[2*idx + 1]) % mpc.modulo;
			psi.y[idx] = (mpc.modulo + sz1.y[2*idx] - sz1.y[2*idx + 1]) % mpc.modulo;
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		std::cout << "****************Sending: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		com->sendingSize = 10*com->sendingSize/16;
		
		nBytes = 10;
		
		// Compact the bytes
		unsigned char *inX = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.x.data()); //new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.y.data()); //new unsigned char[nBytes*size];
		
		unsigned char *ptr = 0;
		int count = 0;
		
		ptr = (unsigned char *)(sz2.x.data());
		for(int idx = 0; idx < sz2.x.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inX[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		count = 0;
		ptr = (unsigned char *)(sz2.y.data());
		for(int idx = 0; idx < sz2.y.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inY[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		t.Tick("Compact bytes");
		
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
		std::vector<Shares> bx1(fieldSize);
		std::vector<Shares> bx2(fieldSize);
		std::vector<Shares> bx3(fieldSize);
		
// 		ArrayEncoder::ZZArray2ByteArray(inX, sz2.x.data(), nBytes, size);
// 		ArrayEncoder::ZZArray2ByteArray(inY, sz2.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, fieldSize);
		sse_trans((uint8_t *)inY, outY, size, fieldSize);
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < fieldSize; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		t.Tick("Share conversion");
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
		t.Tick("RippleCarryAdder");
// 		sum.resize(fieldSize);
		// sum = sum mod p
		tps.modp_fast(sum, p, com);
		
		t.Tick("modp_fast");
		
		std::vector<Shares> sum2(fieldSize);
		for(int idx = 0; idx < fieldSize; idx++)
		{
			sum2[idx] = sum[idx].LeftShiftOne();
		}
		t.Tick("Left Shift");
		if(tps.IsIncreasing(sum2, sum, com))
		{
			std::cout << "Passed verification" << std::endl;
		}
		else
		{
			std::cout << "Failed verification" << std::endl;
		}
		t.Tick("Is Increasing");
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tps.VerifyTriples(rngseed, com);
		t.Tick("PSU proof");
	}
	
	void PSI_CA_Circuit_Arithmetic(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
		
		MPCRingOperations mpc;
		
		RepRingShares sx, sy, sz1, sz2;
		
		std::vector<__uint128_t> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
// 		t.Tick("GenerateInput");
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
// 			t.Tick("rng");
			rng->ComputeRingFx(Fx, data);
		}
		
// 		t.Tick("rng->ComputeRingFx");
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
// 		std::cout << "Shuffling shares" << std::endl;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
		t.Tick("Sorting");
		
// 		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepRingShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = (mpc.modulo + sz1.x[2*idx] - sz1.x[2*idx + 1]) % mpc.modulo;
			psi.y[idx] = (mpc.modulo + sz1.y[2*idx] - sz1.y[2*idx + 1]) % mpc.modulo;
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		std::cout << "****************Sending: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		com->sendingSize = 10*com->sendingSize/16;
		
		nBytes = 10;
		
		RepRingShares b1, b2, b3, mb1, mb2, mb3;
		
		b1.x.resize(size);
		b1.y.resize(size);
		b2.x.resize(size);
		b2.y.resize(size);
		b3.x.resize(size);
		b3.y.resize(size);
		
		for(int idx = 0; idx < size; idx++)
		{
			sz2.x[idx] = (2*sz2.x[idx] % mpc.modulo) % 2;
			sz2.y[idx] = (2*sz2.y[idx] % mpc.modulo) % 2;
			
			
			if(party == Alice)
			{
				b1.x[idx] = sz2.x[idx];
				b1.y[idx] = 0;
				
				b2.x[idx] = 0;
				b2.y[idx] = sz2.y[idx];
				
				b3.x[idx] = 0;
				b3.y[idx] = 0;
			}
			else if(party == Bob)
			{
				b1.x[idx] = 0;
				b1.y[idx] = 0;
				
				b2.x[idx] = sz2.x[idx];
				b2.y[idx] = 0;
				
				b3.x[idx] = 0;
				b3.y[idx] = sz2.y[idx];
			}
			else
			{
				b1.x[idx] = 0;
				b1.y[idx] = sz2.y[idx];
				
				b2.x[idx] = 0;
				b2.y[idx] = 0;
				
				b3.x[idx] = sz2.x[idx];
				b3.y[idx] = 0;
			}
		}
		
		RepRingShares MacKey; mpc.GenerateRandomShares(MacKey, nBytes, 1, com);
		RepRingShares MACs;
		
		mb1 = mpc.repMul(b1, MacKey, nBytes, com);
		mb2 = mpc.repMul(b2, MacKey, nBytes, com);
		mb3 = mpc.repMul(b3, MacKey, nBytes, com);
		
		RepRingShares b1xb2 = mpc.repMul(b1, b2, nBytes, com);
		b1xb2 = mpc.repConstMul(b1xb2, 2, nBytes, com);
		
		RepRingShares mb1xb2 = mpc.repMul(mb1, b2, nBytes, com);
		mb1xb2 = mpc.repConstMul(mb1xb2, 2, nBytes, com);
		
		for(int idx = 0; idx < size; idx++)
		{
			b1xb2.x[idx] -= (b1.x[idx] + b2.x[idx]);
			b1xb2.y[idx] -= (b1.y[idx] + b2.y[idx]);
			
			mb1xb2.x[idx] -= (mb1.x[idx] + mb2.x[idx]);
			mb1xb2.y[idx] -= (mb1.y[idx] + mb2.y[idx]);
		}
		
		RepRingShares b1xb2xb3 = mpc.repMul(b1xb2, b3, nBytes, com);
		b1xb2xb3 = mpc.repConstMul(b1xb2xb3, 2, nBytes, com);
		
		RepRingShares mb1xb2xb3 = mpc.repMul(mb1xb2, b3, nBytes, com);
		mb1xb2xb3 = mpc.repConstMul(mb1xb2xb3, 2, nBytes, com);
		
		for(int idx = 0; idx < size; idx++)
		{
			b1xb2xb3.x[idx] -= (b1xb2.x[idx] + b3.x[idx]);
			b1xb2xb3.y[idx] -= (b1xb2.y[idx] + b3.y[idx]);
			
			mb1xb2xb3.x[idx] -= (mb1xb2.x[idx] + mb3.x[idx]);
			mb1xb2xb3.y[idx] -= (mb1xb2.y[idx] + mb3.y[idx]);
		}
		
		
		
		t.Tick("PSU proof");
	}
	
	
	void PSI_CA_Circuit_Ring(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 16;
		
		MPCNewRingOperations mpc;
		
		RepNewRingShares sx, sy, sz1, sz2;
		
		std::vector<__uint128_t> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
// 		t.Tick("GenerateInput");
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
// 			t.Tick("rng");
			rng->ComputeRingFx(Fx, data);
		}
		
// 		t.Tick("rng->ComputeRingFx");
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
// 		std::cout << "Shuffling shares" << std::endl;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
		t.Tick("Sorting");
		
// 		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepNewRingShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = (mpc.modulo + sz1.x[2*idx] - sz1.x[2*idx + 1]) % mpc.modulo;
			psi.y[idx] = (mpc.modulo + sz1.y[2*idx] - sz1.y[2*idx + 1]) % mpc.modulo;
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		std::cout << "****************Sending: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		com->sendingSize = 10*com->sendingSize/16;
		
		nBytes = 10;
		
		// Compact the bytes
		unsigned char *inX = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.x.data()); //new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size]; //(unsigned char *)(sz2.y.data()); //new unsigned char[nBytes*size];
		
		unsigned char *ptr = 0;
		int count = 0;
		
		ptr = (unsigned char *)(sz2.x.data());
		for(int idx = 0; idx < sz2.x.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inX[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		count = 0;
		ptr = (unsigned char *)(sz2.y.data());
		for(int idx = 0; idx < sz2.y.size(); idx++)
		{
			for(int kdx = 0; kdx < nBytes; kdx++)
			{
				inY[count] = *ptr;
				count++;
				ptr++;
			}
			ptr += (16 - nBytes);
		}
		
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
		std::vector<Shares> bx1(fieldSize);
		std::vector<Shares> bx2(fieldSize);
		std::vector<Shares> bx3(fieldSize);
		
// 		ArrayEncoder::ZZArray2ByteArray(inX, sz2.x.data(), nBytes, size);
// 		ArrayEncoder::ZZArray2ByteArray(inY, sz2.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, fieldSize);
		sse_trans((uint8_t *)inY, outY, size, fieldSize);
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < fieldSize; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
// 		sum.resize(fieldSize);
		// sum = sum mod p
		tps.modp_fast(sum, p, com);
		
		t.Tick("Share conversion");
		
		std::vector<Shares> sum2(fieldSize);
		for(int idx = 0; idx < fieldSize; idx++)
		{
			sum2[idx] = sum[idx].LeftShiftOne();
		}
		
		if(tps.IsIncreasing(sum2, sum, com))
		{
			std::cout << "Passed verification" << std::endl;
		}
		else
		{
			std::cout << "Failed verification" << std::endl;
		}
		
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tps.VerifyTriples(rngseed, com);
		t.Tick("PSU proof");
	}
	
	
	void PSI_CA_Circuit2(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
// 		if(nItems > 1000000)
// 		{
// 			ZZ pp = (ZZ(1) << 96) - 17;
// 			ZZ_p::init(pp);
// 			fieldSize = 96; nBytes = 12;
// 		}
		
		MPCOperations mpc;
		
		RepArithShares sx, sy, sz1, sz2;
		
		std::vector<ZZ_p> Fx;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
		}
		t.Tick("GenerateInput");
		
		if(party != Charlie)
		{
			std::vector<unsigned char> seed  = com->getSharedRandomSeed(party, partner);
			
			AESRNG *rng = new AESRNG(seed.data());
			
			rng->ComputeFx(Fx, data);
		}
		t.Tick("ComputeFx");
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Data commitment");
		
		std::cout << "Shuffling shares" << std::endl;
		mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		
		t.Tick("Sorting");
		
		std::cout << "PSI size: " << sz1.x.size()/2 << std::endl;
		RepArithShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = sz1.x[2*idx] - sz1.x[2*idx + 1];
			psi.y[idx] = sz1.y[2*idx] - sz1.y[2*idx + 1];
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		unsigned char *inX = new unsigned char[nBytes*size];
		unsigned char *inY = new unsigned char[nBytes*size];
		unsigned char *outX = new unsigned char[nBytes*size];
		unsigned char *outY = new unsigned char[nBytes*size];
		
		std::vector<Shares> bx1(fieldSize);
		std::vector<Shares> bx2(fieldSize);
		std::vector<Shares> bx3(fieldSize);
		
		ArrayEncoder::ZZArray2ByteArray(inX, sz2.x.data(), nBytes, size);
		ArrayEncoder::ZZArray2ByteArray(inY, sz2.y.data(), nBytes, size);
		
		sse_trans((uint8_t *)inX, outX, size, fieldSize);
		sse_trans((uint8_t *)inY, outY, size, fieldSize);
		
		int numBlks = size/64;
		
		uint64_t *out64X = (uint64_t *)outX;
		uint64_t *out64Y = (uint64_t *)outY;
		
		for(int idx = 0; idx < fieldSize; idx++)
		{
			std::vector<uint64_t> tempX(out64X, out64X + numBlks);
			std::vector<uint64_t> tempY(out64Y, out64Y + numBlks);
			
			if(party == Alice)
			{
				bx1[idx].s = tempX; 
				bx1[idx].t = tempX;
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t.resize(numBlks, 0);
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t = tempY;
			}
			else if(party == Bob)
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t = tempY;
				bx2[idx].s = tempX;
				bx2[idx].t = tempX;
				bx3[idx].s.resize(numBlks, 0);
				bx3[idx].t.resize(numBlks, 0);
			}
			else
			{
				bx1[idx].s.resize(numBlks, 0);
				bx1[idx].t.resize(numBlks, 0);
				bx2[idx].s.resize(numBlks, 0);
				bx2[idx].t = tempY;
				bx3[idx].s = tempX;
				bx3[idx].t = tempX;
			}
			
			out64X += numBlks; out64Y += numBlks;
		}
		
		delete [] inX;
		delete [] inY;
		delete [] outX;
		delete [] outY;
		
		Triplets tps;
		
		// sum = bx1 + bx2
		std::vector<Shares> sum = tps.RippleCarryAdder(bx1, bx2, bx3, com);
		
// 		sum.resize(sum.size() - 2);
		// sum = sum mod p
		tps.modp_fast(sum, p, com);
		
		t.Tick("Share conversion");
		
		std::vector<Shares> sum2(fieldSize);
		for(int idx = 0; idx < fieldSize; idx++)
		{
			sum2[idx] = sum[idx].LeftShiftOne();
		}
		
		if(tps.IsIncreasing(sum2, sum, com))
		{
			std::cout << "Passed verification" << std::endl;
		}
		else
		{
			std::cout << "Failed verification" << std::endl;
		}
		
		std::vector<unsigned char> rngseed = machine->getCommonSeed(com);
		
		tps.VerifyTriples(rngseed, com);
		t.Tick("PSU proof");
	}
	
	void PSI_CA_Cristefaro(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		ZZ_pContext context;
		
		// Distribute (shared/whole) prf keys
		// 1024 bit prime number
		ZZ p = conv<ZZ>("124325339146889384540494091085456630009856882741872806181731279018491820800119460022367403769795008250021191767583423221479185609066059226301250167164084041279837566626881119772675984258163062926954046545485368458404445166682380071370274810671501916789361956272226105723317679562001235501455748016154805420913");
		
		// 160 bit prime
		ZZ q = conv<ZZ>("1399252811935680595399801714158014275474696840019");
		
		// Generator
		ZZ g = conv<ZZ>("115740200527109164239523414760926155534485715860090261532154107313946218459149402375178179458041461723723231563839316251515439564315555249353831328479173170684416728715378198172203100328308536292821245983596065287318698169565702979765910089654821728828592422299160041156491980943427556153020487552135890973413");
		
		
		if(party != Charlie)
		{
			ZZ_p::init(p);
			
			context.save();
			ZZ_p::init(q);
			ZZ_p R_cs;
			random(R_cs);
			context.restore();
			// Bytes required for each f(x)
			int nBytes = 128; //ceil(40 + 2*20)/8;
			
			std::cout << "nItems: " << nItems << "\tnBytes: " << nBytes << "\t" << ZZ_p::ModulusSize() << std::endl;
			unsigned char *bytes = new unsigned char[nBytes*nItems];
			
			std::vector<uint64_t> data;
			std::vector<ZZ_p> Fx, A, A1, A2, Fx1, Fx2;
			
			std::cout << "Generate input" << std::endl;
			GenerateInput(party, data, nItems, rate);
		
			std::vector<unsigned char> seed0  = machine->getSharedRandomSeed(party, partner, com);
			std::vector<unsigned char> seed1  = machine->getSharedRandomSeed(party, partner, com);
			
			AESRNG *rngH = new AESRNG(seed0.data());
			AESRNG *rngH1 = new AESRNG(seed1.data());
			
			// hc_i
			rngH->ComputeFx(Fx, data);
			
			
			if(party == Alice)
			{
				context.save();
				ZZ_p::init(q);
				ZZ_p Rc, iRc;
				random(Rc); inv(iRc, Rc);
				context.restore();
				A.resize(Fx.size());
				
				context.save(); 
				ZZ_p::init(p);
				
				// a_i
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					power(A[idx], Fx[idx], rep(Rc));
				}
				
				context.restore();
				
				ZZ_p PoKAlice = conv<ZZ_p>(1); 
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					PoKAlice = PoKAlice*A[idx];
				}
				
				ArrayEncoder::ZZArray2ByteArray(bytes, A.data(), nBytes, nItems);
				com->SendEvaluationPartner(bytes, nBytes*nItems);
				
				// Receive A1 from Bob
				A1.resize(nItems);
				com->AwaitEvaluationPartner(bytes, nBytes*nItems);
				ArrayEncoder::ByteArray2ZZArray(A1.data(), bytes, nBytes, nItems);
				
				A2.resize(nItems);
				// a_i''
				context.save();
				ZZ_p::init(q);
				ZZ_p Rc1; random(Rc1);
				context.restore();
				
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					power(A2[idx], A1[idx], rep(Rc1));
				}
				
				std::random_shuffle(A2.begin(), A2.end());
				
				ArrayEncoder::ZZArray2ByteArray(bytes, A2.data(), nBytes, nItems);
				com->SendEvaluationPartner(bytes, nBytes*nItems);
				
				// Receive A3 from Bob 
				std::vector<ZZ_p> A3(nItems);
				com->AwaitEvaluationPartner(bytes, nBytes*nItems);
				ArrayEncoder::ByteArray2ZZArray(A3.data(), bytes, nBytes, nItems);
				
				context.save();
				ZZ_p::init(q);
				ZZ_p iRc1; inv(iRc1, Rc1);
				context.restore();
				
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					power(A3[idx], A3[idx], rep(iRc1));
				}
				
				std::vector<ZZ> temp1(nItems), temp2(nItems);
				
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					temp1[idx] = rep(A3[idx]);
					temp2[idx] = rep(A[idx]);
				}
				
				std::sort(temp1.begin(), temp1.end());
				std::sort(temp2.begin(), temp2.end());
				
				for(int idx = 0; idx < 5; idx++)
				{
					std::cout << temp1[idx] << " " << temp2[idx] << std::endl;
				}
				
// 				assert(temp1 == temp2);
				
				for(int idx = 0; idx < nItems; idx++)
				{
					power(A1[idx], A1[idx], rep(iRc));
				}
				
// 				rngH1->ComputeFx(Fx, A1, nBytes);
				unsigned char *a = new unsigned char[nBytes];
				
				for(int idx = 0; idx < nItems; idx++)
				{
					BytesFromZZ(a, rep(A1[idx]), nBytes);
					std::string str(reinterpret_cast<const char *>(a));
					unsigned char *val = rngH1->AES(str);
					ZZ temp;
					ZZFromBytes(temp, val, nBytes);
					conv(Fx[idx], temp);
				}
				
				com->AwaitEvaluationPartner(bytes, nBytes*nItems);
				ArrayEncoder::ByteArray2ZZArray(A.data(), bytes, nBytes, nItems);
				
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					temp1[idx] = rep(A1[idx]);
					temp2[idx] = rep(A[idx]);
				}
				
				temp1.insert(temp1.end(), temp2.begin(), temp2.end());
				
				std::sort(temp1.begin(), temp1.end());
				
				temp1.erase(std::unique(temp1.begin(), temp1.end()), temp1.end());
				
				std::cout << "Union Cardinality: " << temp1.size() << std::endl;
			}
			else
			{
				context.save();
				ZZ_p::init(q);
				ZZ_p Rs;
				random(Rs);
				context.restore();
				A.resize(Fx.size());
				
				for(int idx = 0; idx < Fx.size(); idx++)
				{
					power(A[idx], Fx[idx], rep(Rs));
				}
				
				std::cout << "Compute ts_j" << std::endl;
// 				rngH1->ComputeFx(Fx, A, nBytes);
				
				unsigned char *a = new unsigned char[nBytes];
				
				for(int idx = 0; idx < nItems; idx++)
				{
					BytesFromZZ(a, rep(A[idx]), nBytes);
					std::string str(reinterpret_cast<const char *>(a));
					unsigned char *val = rngH1->AES(str);
					ZZ temp;
					ZZFromBytes(temp, val, nBytes);
					conv(Fx[idx], temp);
				}
				std::cout << "Waiting for data from Alice" << std::endl;
				
				com->AwaitEvaluationPartner(bytes, nBytes*nItems);
				A1.resize(nItems);
				
				ArrayEncoder::ByteArray2ZZArray(A1.data(), bytes, nBytes, nItems);
				for(int idx = 0; idx < A1.size(); idx++)
				{
					power(A1[idx], A1[idx], rep(R_cs));
				}
				
				ZZ_p PoKBob = conv<ZZ_p>(1);
				for(int idx = 0; idx < A1.size(); idx++)
				{
					PoKBob = PoKBob*A1[idx];
				}
				
				// Shuffle the set 
				std::random_shuffle(A1.begin(), A1.end());
				
				ArrayEncoder::ZZArray2ByteArray(bytes, A1.data(), nBytes, nItems);
				com->SendEvaluationPartner(bytes, nBytes*nItems);
				
				// Send proof of knowledge here
				
				
				// Receive A2 from Alice
				std::vector<ZZ_p> & A2 = A1;
				com->AwaitEvaluationPartner(bytes, nBytes*nItems);
				ArrayEncoder::ByteArray2ZZArray(A2.data(), bytes, nBytes, nItems);
				
				context.save();
				ZZ_p::init(q);
				ZZ_p iRs; 
				inv(iRs, Rs);
				context.restore();
				
				std::vector<ZZ_p> & A3 = A2;
				
				for(int idx = 0; idx < A2.size(); idx++)
				{
					power(A3[idx], A2[idx], rep(iRs));
				}
				
				std::random_shuffle(A3.begin(), A3.end());
				
				ArrayEncoder::ZZArray2ByteArray(bytes, A3.data(), nBytes, nItems);
				com->SendEvaluationPartner(bytes, nBytes*nItems);
				
				ArrayEncoder::ZZArray2ByteArray(bytes, A.data(), nBytes, nItems);
				com->SendEvaluationPartner(bytes, nBytes*nItems);
			}
		}
		
		std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
	}
	
// 	int PSI_Poly(int party, int nItems, Communication *com, Machine *machine, float rate)
// 	{
// 		ZZ pp = (ZZ(1) << 80) - 65; // n = 2^20
// 		ZZ_p::init(pp);
// 		
// 		// Bytes required for each f(x)
// 		int nBytes = ceil(40 + 2*20)/8;
// 		
// 		std::cout << "nItems: " << nItems << "\tnBytes: " << nBytes << "\t" << ZZ_p::ModulusSize() << std::endl;
// 		
// 		std::vector<uint64_t> data;
// 		std::vector<ZZ_p> dataFx, dataPSI;
// 		
// 		std::cout << "Generate input" << std::endl;
// 		GenerateInput(party, data, nItems, rate);
// 		
// 		Timer t;
// 		
// 		int psiSize = 0, unionSize = 0;
// 		unsigned char *bytes = new unsigned char[nBytes*2*nItems];
// 		
// 		if(party == Alice || party == Bob)
// 		{
// 			std::vector<unsigned char> rand  = machine->getSharedRandomSeed(party, partner, com);
// 			AESRNG *rng = new AESRNG(rand.data());
// 			rng->ComputeFx(dataFx, data);
// 		}
// 		
// 		std::vector<unsigned __int128> coeffs = rng->GetUInt128Array(2*nItems);
// 			
// 			ZZ_pX f, g;
// 			f.SetLength(unionSize);
// 			g.SetLength(psiSize);
// 			
// 			for(int idx = 0; idx < unionSize; idx++)
// 			{
// 				f[idx] = coeffs[idx];
// 			}
// 			
// 		return 0;
// 	}
// 	
	int PSI_CA_Polynomial(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
// 		SetNumThreads(8);
		ZZ pp = (ZZ(1) << 80) - 65; // n = 2^20
		ZZ_p::init(pp);
		
		// Bytes required for each f(x)
		int nBytes = ceil(40 + 2*20)/8;
		
		std::cout << "nItems: " << nItems << "\tnBytes: " << nBytes << "\t" << ZZ_p::ModulusSize() << std::endl;
		
		std::vector<uint64_t> data;
		std::vector<ZZ_p> dataFx, dataPSI;
		
		std::cout << "Generate input" << std::endl;
		GenerateInput(party, data, nItems, rate);
		
		Timer t;
		
		int psiSize = 0, unionSize = 0;
		unsigned char *bytes = new unsigned char[nBytes*2*nItems];
		
		// Alice & Bob Compute f(k, x) and send to Charlie
		if(party == Alice || party == Bob)
		{
			std::vector<unsigned char> rand  = machine->getSharedRandomSeed(party, partner, com);
			AESRNG *rng = new AESRNG(rand.data());
			rng->ComputeFx(dataFx, data);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, dataFx.data(), nBytes, nItems);
			com->SendCharlie(bytes, nBytes*nItems);
			t.Tick("Encode data");
			
			com->AwaitCharlie((unsigned char *)(&unionSize), sizeof(int));
			psiSize = 2*nItems - unionSize;
			
			
			std::cout << "Degree f(union) & g(psi): " << (unionSize-1) << " " << (psiSize-1) << std::endl;
			
			std::vector<unsigned __int128> coeffs = rng->GetUInt128Array(2*nItems);
			ZZ_pX f, g;
			f.SetLength(unionSize);
			g.SetLength(psiSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				f[idx] = coeffs[idx];
			}
			
			for(int idx = 0; idx < psiSize; idx++)
			{
				g[idx] = coeffs[idx + unionSize];
			}
			
			coeffs.resize(0);
			coeffs.shrink_to_fit();
			
			printPolynomial(f, "Random polynomial F(x)", 1);
			printPolynomial(g, "Random polynomial G(x)", 1);
			
			std::cout << "Evaluate the polynomial over " << nItems << " points" << std::endl;
			// Sending data for the union
			ZZ_p *evals = new ZZ_p[nItems];
			
			PolynomialTree *tree = new PolynomialTree();
			buildTree(tree, dataFx.data(), nItems);
			t.Tick("Build tree");
			
			downtreeDivision(evals, tree, f, 0, nItems);
			t.Tick("Evaluate");
	
// 			evaluatePolynomial(evals, f, dataFx.data(), nItems);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, evals, nBytes, nItems);
			com->SendCharlie(bytes, nItems*nBytes);
			
			t.Tick("Sending f");
			
			// Sending data for the intersection
			rng->ComputeFx(dataPSI, dataFx, nBytes);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, dataPSI.data(), nBytes, nItems);
			
			if(party == Bob)
			{
// 				evaluatePolynomial(evals, g, dataFx.data(), nItems);
				downtreeDivision(evals, tree, g, 0, nItems);
				unsigned char *temp = new unsigned char[nBytes*nItems];
				ArrayEncoder::ZZArray2ByteArray(temp, evals, nBytes, nItems);
				for(int idx = 0; idx < nBytes*nItems; idx++)
				{
					bytes[idx] ^= temp[idx];
				}
				delete [] temp;
			}
			
			t.Tick("Evaluate g");
			
			com->SendCharlie(bytes, nItems*nBytes);
			
			t.Tick("Sending g");
			
			// Get commitment of f[0] from Charlie
			std::vector<unsigned char> commitment(COMMITMENT_LENGTH);
			com->AwaitCharlie(commitment.data(), commitment.size());
			
			// Send f(x) and g(x) to Charlie
			ArrayEncoder::ZZpX2ByteArray(bytes, f, nBytes, unionSize);
			com->SendCharlie(bytes, nBytes*unionSize);
			
			ArrayEncoder::ZZpX2ByteArray(bytes, g, nBytes, psiSize);
			com->SendCharlie(bytes, nBytes*psiSize);
			
			// Send rng to Charlie
			com->SendCharlie(rand.data(), rand.size());
			
			std::vector<unsigned char> seed(SEED_LENGTH), msg(nBytes);
			com->AwaitCharlie(seed.data(), seed.size());
			com->AwaitCharlie(msg.data(), msg.size());
			HashCommitment hc;
			assert(hc.Verification(commitment, msg, seed));
			
			t.Tick("Verification");
			
			delete rng;
			delete [] evals;
		}
		else if(party == Charlie)
		{
			// Charlie count the size of the union and send it to Alice & Bob 
			unsigned char *ptr = bytes;
			
			com->AwaitAlice(bytes, nBytes*nItems);
			com->AwaitBob((bytes + nBytes*nItems), nBytes*nItems);
			
			dataFx.resize(2*nItems);

			std::vector<ZZ> F_k1_xy(2*nItems);
			
			ArrayEncoder::ByteArray2ZZArray(F_k1_xy.data(), bytes, nBytes, 2*nItems);
			
			// Sort array according to f(k,x) and remember the permutation
			std::cout << "Sort data based on f(k,x)" << std::endl;
			
			std::vector<ZZ> F_k1_xy_copy = F_k1_xy;
			
			vector<int> index(2*nItems);
			for (int idx = 0 ; idx != index.size() ; idx++) {
				index[idx] = idx;
			}
			
			// Get the permutation
			std::sort(index.begin(), index.end(), [&](const int& a, const int& b) {return (F_k1_xy[a] < F_k1_xy[b]);});
			
			// Sort
			std::sort(F_k1_xy.begin(), F_k1_xy.end());
			
			// Remove duplicated items
			F_k1_xy.erase(std::unique(F_k1_xy.begin(), F_k1_xy.end()), F_k1_xy.end());
			
			unionSize = F_k1_xy.size();
			psiSize = 2*nItems - unionSize;
			
			dataFx.resize(unionSize);
			for(int idx = 0; idx < unionSize; idx++)
			{
				 conv(dataFx[idx], F_k1_xy[idx]);
			}
			
			t.Tick("Local Sort");
			
			com->SendAlice((unsigned char *)(&unionSize), sizeof(int));
			com->SendBob((unsigned char *)(&unionSize), sizeof(int));
			
			///-------------------------------------Union Interpolation-----------------------------------------------///
			
			///----------------Build Tree From F(k1,x)-----------------///
			std::cout << "Union Interpolation" << std::endl;
			// Build tree and evaluate derivative polynomial
			PolynomialTree *tree = new PolynomialTree();
			buildTree(tree, dataFx.data(), unionSize);
			
			t.Tick("Build tree");
			
			ZZ_pX derivativePolynomial;
			derivativePolynomial.SetLength(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				derivativePolynomial[idx] = (tree->val[idx + 1]*(idx + 1));
			}
			
			ZZ_p *derivativeEvals = new ZZ_p[unionSize];
			
			downtreeDivision(derivativeEvals, tree, derivativePolynomial, 0, unionSize);
			
			
			t.Tick("Downtree Division");
			
			ZZ_p *evals = new ZZ_p[2*nItems];
			
			///----------------Interpolate from p1(F(k1,x))-----------------///
			
			com->AwaitAlice(bytes, nItems*nBytes);
			com->AwaitBob((bytes + nItems*nBytes), nItems*nBytes);
			
			ptr = bytes;
			
			t.Tick("Wait time");
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				conv(evals[idx], ZZFromBytes(ptr, nBytes));
				ptr += nBytes;
			}
			
			std::cout << "Verify (x, p(x))" << std::endl;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(F_k1_xy_copy[index[idx - 1]] == F_k1_xy_copy[index[idx]])
				{
					assert(evals[index[idx - 1]] == evals[index[idx]]);
				}
			}
			
			ZZ_p *unionEvals = new ZZ_p[unionSize];
			
			int count = 0; 
			int count2 = 0;
			
			unionEvals[count] = evals[index[0]]; count++;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(F_k1_xy_copy[index[idx - 1]] != F_k1_xy_copy[index[idx]])
				{
					unionEvals[count] = evals[index[idx]];
					count++;
				}
			}
			
			assert(count == unionSize);
			
			t.Tick("Verify data");
			
			reconstruct(unionEvals, derivativeEvals, tree, 0, unionSize);
			
			ZZ_pX f = tree->interp;
			
			t.Tick("Reconstruction time for f");
			
			printPolynomial(f, "Interpolated polynomial", 1);
			
			///-------------------------------------PSI Interpolation-----------------------------------------------///
			std::cout << "PSI Interpolation" << std::endl;
			// Get the index of the matches
			dataFx.resize(psiSize);
			count = 0; count2 = 0;
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(F_k1_xy_copy[index[idx - 1]] == F_k1_xy_copy[index[idx]])
				{
					conv(dataFx[count], F_k1_xy_copy[index[idx-1]]);
					count++;
				}
			}
			
			assert(count == psiSize);
			
			///----------------Build Tree From F(k1,x)-----------------///
			
			PolynomialTree *psiTree = new PolynomialTree();
			buildTree(psiTree, dataFx.data(), psiSize);
			t.Tick("Build psi tree");
			ZZ_pX derivativePolynomialPSI;
			derivativePolynomialPSI.SetLength(psiSize);
			
			for(int idx = 0; idx < psiSize; idx++)
			{
				derivativePolynomialPSI[idx] = (psiTree->val[idx + 1]*(idx + 1));
			}
			
			ZZ_p *derivativeEvalsPSI = new ZZ_p[psiSize];
			
			downtreeDivision(derivativeEvalsPSI, psiTree, derivativePolynomialPSI, 0, psiSize);
			t.Tick("Downtree Division");
			///----------------Interpolate from p2(F(k1,x))-----------------///
			
			// Interpolate the intersection polynomial
			// Receiving F(k2, x) and F(k2, y) \xor g(y)
			com->AwaitAlice(bytes, nBytes*nItems);
			com->AwaitBob(bytes + nBytes*nItems, nBytes*nItems);
			
			// Hold F(k2, x~)
			std::vector<ZZ_p> AliceFxk2(nItems);
			std::vector<ZZ_p> BobFxPxk2(nItems);
			
			ArrayEncoder::ByteArray2ZZArray(AliceFxk2.data(), bytes, nBytes, nItems);
			ArrayEncoder::ByteArray2ZZArray(BobFxPxk2.data(), bytes + nBytes*nItems, nBytes, nItems);
			
			// Reconstruct g(y) for y in the intersection
			ZZ_p *psiEvals = new ZZ_p[psiSize];
			
			unsigned char *val = new unsigned char[nBytes];
			count = 0;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(F_k1_xy_copy[index[idx - 1]] == F_k1_xy_copy[index[idx]])
				{
					unsigned char *ptr1 = bytes + nBytes*index[idx-1];
					unsigned char *ptr2 = bytes + nBytes*index[idx];
					for(int idx = 0; idx < nBytes; idx++) val[idx] = ptr1[idx]^ptr2[idx];
					conv(psiEvals[count], ZZFromBytes(val, nBytes));
					count++;
				}
			}
			
			delete [] val;
			
			reconstruct(psiEvals, derivativeEvalsPSI, psiTree, 0, psiSize);
			
			ZZ_pX g = psiTree->interp;
			
			printPolynomial(g, "Interpolated polynomial G", 1);
			
			t.Tick("Reconstruction g");
			
			// Send Alice & Bob f[0]
			std::vector<unsigned char> seed;
			
			std::vector<unsigned char> msg(nBytes);
			BytesFromZZ(msg.data(), rep(f[0]), nBytes);
			
			HashCommitment hc;
			
			std::vector<unsigned char> commitment = hc.Commit(msg, seed);
			
			com->SendAlice(commitment.data(), commitment.size());
			com->SendBob(commitment.data(), commitment.size());
			
			// Receiver p(x) from Alice & Bob 
			unsigned char *pA = new unsigned char[nBytes*unionSize];
			unsigned char *pB = new unsigned char[nBytes*unionSize];
			unsigned char *pC = new unsigned char[nBytes*psiSize];
			unsigned char *pD = new unsigned char[nBytes*psiSize];
			com->AwaitAlice(pA, nBytes*unionSize);
			com->AwaitBob(pB, nBytes*unionSize);
			com->AwaitAlice(pC, nBytes*psiSize);
			com->AwaitBob(pD, nBytes*psiSize);
			
			assert(0 == std::memcmp(pA, pB, nBytes*unionSize));
			assert(0 == std::memcmp(pC, pD, nBytes*psiSize));
			
			// Verify that pA = f
			
			ptr = pA;
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				ZZ temp;
				ZZFromBytes(temp, ptr, nBytes);
				assert(temp == rep(f[idx]));
				ptr += nBytes;
			}
			
			ptr = pC;
			
			for(int idx = 0; idx < psiSize; idx++)
			{
				ZZ temp;
				ZZFromBytes(temp, ptr, nBytes);
				assert(temp == rep(g[idx]));
				ptr += nBytes;
			}
			
			// Receive k2 from PA and PB
			std::vector<unsigned char> randAlice(SEED_LENGTH/2), randBob(SEED_LENGTH/2);
			com->AwaitAlice(randAlice.data(), SEED_LENGTH/2);
			com->AwaitBob(randBob.data(), SEED_LENGTH/2);
			
			assert(randAlice == randBob);
			
			AESRNG *rng = new AESRNG(randAlice.data());
			
			ZZ_p *evalsPx = new ZZ_p[nItems];
			
			// Hold xhat = F(k1, x)
			std::vector<ZZ_p> AliceFxk1(nItems), BobFxk1(nItems);
			
			for(int idx = 0; idx < nItems; idx++)
			{
				conv(AliceFxk1[idx], F_k1_xy_copy[idx]);
				conv(BobFxk1[idx], F_k1_xy_copy[nItems + idx]);
			}
			
			// p1(yhat)
			evaluatePolynomial(evalsPx, g, BobFxk1.data(), nItems);
			
			// F(k_2, xhat)
			rng->ComputeFx(AliceFxk1, AliceFxk1, nBytes);
			rng->ComputeFx(BobFxk1, BobFxk1, nBytes);
			
			// Verify V1 is correct
			for(int idx = 0; idx < nItems; idx++)
			{
				assert(AliceFxk1[idx] == AliceFxk2[idx]);
			}
				
			ArrayEncoder::ZZArray2ByteArray(bytes, BobFxk1.data(), nBytes, nItems);
			ArrayEncoder::ZZArray2ByteArray(bytes + nBytes*nItems, evalsPx, nBytes, nItems);
			
			for(int idx = 0; idx < nBytes*nItems; idx++)
			{
				bytes[idx] ^= bytes[idx + nBytes*nItems];
			}
			
			ArrayEncoder::ByteArray2ZZArray(BobFxk1.data(), bytes, nBytes, nItems);
			
			//Verify W1 is correct
			for(int idx = 0; idx < nItems; idx++)
			{
				assert(BobFxk1[idx] == BobFxPxk2[idx]);
			}
			
			
			delete [] pA;
			delete [] pB;
			delete [] pC;
			delete [] pD;
			delete [] derivativeEvals;
			delete [] unionEvals;
			
			com->SendAlice(seed.data(), seed.size());
			com->SendBob(seed.data(), seed.size());
			
			com->SendAlice(msg.data(), msg.size());
			com->SendBob(msg.data(), msg.size());
			t.Tick("Verification g(x)");
		}
		
		delete [] bytes;
		
		std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
		
		com->sendingSize = 0;
		com->receivingSize = 0;
		
		return unionSize;
	}
	
	void PSI_CA_Hybrid(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
		
		MPCRingOperations mpc;
		
		RepRingShares sx, sy, sz1, sz2;
		
		std::vector<__uint128_t> Fx2;
		std::vector<ZZ_p> Fx;
		AESRNG *rng;
		std::vector<unsigned char> rand;
		
		std::vector<uint64_t> data;
		if(party != Charlie)
		{
			GenerateInput(party, data, nItems, rate);
			rand  = com->getSharedRandomSeed(party, partner);
			
			rng = new AESRNG(rand.data());
// 			t.Tick("rng");
			rng->ComputeRingFx(Fx2, data);
			rng->ComputeFx(Fx, data);
		}
		
		mpc.MPC_ShareData(sx, Alice, Fx2, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx2, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		std::cout << "Share data comm: " << com->sendingSize << std::endl;
		
		t.Tick("Secret share data");
		
// 		std::cout << "Shuffling shares" << std::endl;
// 		mpc.MPC_Sort(sz1, sz2, Fx2, sx, nBytes, com);
		if(party == Charlie)
		{
			mpc.MPC_Sort(sz1, sz2, Fx2, sx, nBytes, com);
			Fx.resize(Fx2.size());
			for(int idx = 0; idx < Fx.size(); idx++)
			{
				conv(Fx[idx], Fx2[idx]);
			}
		}
		else
		{
			std::vector<__uint128_t> temp;
			mpc.MPC_Sort(sz1, sz2, temp, sx, nBytes, com);
		}
		t.Tick("3pc sort");
		
		RepRingShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = (mpc.modulo + sz1.x[2*idx] - sz1.x[2*idx + 1]) % mpc.modulo;
			psi.y[idx] = (mpc.modulo + sz1.y[2*idx] - sz1.y[2*idx + 1]) % mpc.modulo;
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("PSI proof");
		
		if(sz1.x.size()/2 == nItems) return;
		// Convert shares in Z2 to binary shares
		int size = sz2.x.size();
		
		std::cout << "****************Sending: " << com->sendingSize/1024.0/1024.0 << std::endl;
// 		com->sendingSize = 10*com->sendingSize/16;
		
		nBytes = 10;
		
// 		if(sz1.x.size()/2 == nItems) return;
		
		// Prove union bound by polynomial approach
		int psiSize = psi.x.size();
		int unionSize = 2*nItems - psiSize;
		
		unsigned char *bytes = new unsigned char[nBytes*2*nItems];
		
		// Alice & Bob Compute f(k, x) and send to Charlie
		if(party == Alice || party == Bob)
		{
// 			std::cout << "Union/psi size: " << unionSize << " " << psiSize << " " << (unionSize + psiSize) <<  std::endl;
			
			std::vector<unsigned __int128> coeffs = rng->GetUInt128Array(unionSize);
			
			ZZ_pX f;
			f.SetLength(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				f[idx] = coeffs[idx];
			}
			
			coeffs.resize(0);
			coeffs.shrink_to_fit();
			
			printPolynomial(f, "Random polynomial F(x)", 1);
			
// 			std::cout << "Evaluate the polynomial over " << nItems << " points" << std::endl;
			// Sending data for the union
			ZZ_p *evals = new ZZ_p[nItems];
			evaluatePolynomial(evals, f, Fx.data(), nItems);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, evals, nBytes, nItems);
			com->SendCharlie(bytes, nItems*nBytes);
			
			// Get commitment of f[0] from Charlie
			std::vector<unsigned char> commitment(COMMITMENT_LENGTH);
			com->AwaitCharlie(commitment.data(), commitment.size());
			
			// Send f(x) to Charlie
			ArrayEncoder::ZZpX2ByteArray(bytes, f, nBytes, unionSize);
			com->SendCharlie(bytes, nBytes*unionSize);
			
			// Send rng to Charlie
			com->SendCharlie(rand.data(), rand.size());
			
			std::vector<unsigned char> seed(SEED_LENGTH), msg(nBytes);
			com->AwaitCharlie(seed.data(), seed.size());
			com->AwaitCharlie(msg.data(), msg.size());
			HashCommitment hc;
			assert(hc.Verification(commitment, msg, seed));
			
			delete rng;
			delete [] evals;
		}
		else if(party == Charlie)
		{
			// Charlie count the size of the union and send it to Alice & Bob 
			unsigned char *ptr = bytes;
			
			// Sort array according to f(k,x) and remember the permutation
// 			std::cout << "Sort data based on f(k,x)" << std::endl;
			std::vector<ZZ_p> Fx_copy = Fx;
			std::vector<ZZ> dataFx(2*nItems);
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				conv(dataFx[idx], Fx[idx]);
			}
			
			vector<int> index(2*nItems);
			for (int idx = 0 ; idx != index.size() ; idx++) {
				index[idx] = idx;
			}
			
			// Get the permutation
			std::sort(index.begin(), index.end(), [&](const int& a, const int& b) {return (dataFx[a] < dataFx[b]);});
			
			// Sort
			std::sort(dataFx.begin(), dataFx.end());
			
			// Remove duplicated items
			dataFx.erase(std::unique(dataFx.begin(), dataFx.end()), dataFx.end());
			Fx.resize(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				conv(Fx[idx], dataFx[idx]);
			}
			
			///-------------------------------------Union Interpolation-----------------------------------------------///
			Timer t2;
			
			///----------------Build Tree From F(k1,x)-----------------///
			
			// Build tree and evaluate derivative polynomial
			PolynomialTree *tree = new PolynomialTree();
			buildTree(tree, Fx.data(), unionSize);
			t2.Tick("Build tree");
			ZZ_pX derivativePolynomial;
			derivativePolynomial.SetLength(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				derivativePolynomial[idx] = (tree->val[idx + 1]*(idx + 1));
			}
			
			ZZ_p *derivativeEvals = new ZZ_p[unionSize];
			
			downtreeDivision(derivativeEvals, tree, derivativePolynomial, 0, unionSize);
			
			ZZ_p *evals = new ZZ_p[2*nItems];
			
			///----------------Interpolate from p1(F(k1,x))-----------------///
			
			com->AwaitAlice(bytes, nItems*nBytes);
			com->AwaitBob((bytes + nItems*nBytes), nItems*nBytes);
			
			ptr = bytes;
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				conv(evals[idx], ZZFromBytes(ptr, nBytes));
				ptr += nBytes;
			}
			
// 			std::cout << "Verify (x, p(x))" << std::endl;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(Fx_copy[index[idx - 1]] == Fx_copy[index[idx]])
				{
					assert(evals[index[idx - 1]] == evals[index[idx]]);
				}
			}
			
			ZZ_p *unionEvals = new ZZ_p[unionSize];
			
			int count = 0; 
			int count2 = 0;
			
			unionEvals[count] = evals[index[0]]; count++;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(Fx_copy[index[idx - 1]] != Fx_copy[index[idx]])
				{
					unionEvals[count] = evals[index[idx]];
					count++;
				}
			}
			
			assert(count == unionSize);
			
			reconstruct(unionEvals, derivativeEvals, tree, 0, unionSize);
			
			ZZ_pX f = tree->interp;
			
			printPolynomial(f, "Interpolated polynomial", 1);
			
			t2.Tick("Reconstruction");
			
			// Send Alice & Bob f[0]
			std::vector<unsigned char> seed;
			
			std::vector<unsigned char> msg(nBytes);
			BytesFromZZ(msg.data(), rep(f[0]), nBytes);
			
			HashCommitment hc;
			
			std::vector<unsigned char> commitment = hc.Commit(msg, seed);
			
			com->SendAlice(commitment.data(), commitment.size());
			com->SendBob(commitment.data(), commitment.size());
			
			// Receiver p(x) from Alice & Bob 
			unsigned char *pA = new unsigned char[nBytes*unionSize];
			unsigned char *pB = new unsigned char[nBytes*unionSize];
			unsigned char *pC = new unsigned char[nBytes*psiSize];
			unsigned char *pD = new unsigned char[nBytes*psiSize];
			com->AwaitAlice(pA, nBytes*unionSize);
			com->AwaitBob(pB, nBytes*unionSize);
			
			assert(0 == std::memcmp(pA, pB, nBytes*unionSize));
			
			// Verify that pA = f
			
			ptr = pA;
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				ZZ temp;
				ZZFromBytes(temp, ptr, nBytes);
// 				assert(temp == rep(f[idx]));
				ptr += nBytes;
			}
			
			// Receive k2 from PA and PB
			std::vector<unsigned char> randAlice(SEED_LENGTH/2), randBob(SEED_LENGTH/2);
			com->AwaitAlice(randAlice.data(), SEED_LENGTH/2);
			com->AwaitBob(randBob.data(), SEED_LENGTH/2);
			
			assert(randAlice == randBob);
			
			delete [] pA;
			delete [] pB;
			delete [] pC;
			delete [] pD;
			delete [] derivativeEvals;
			delete [] unionEvals;
			
			com->SendAlice(seed.data(), seed.size());
			com->SendBob(seed.data(), seed.size());
			
			com->SendAlice(msg.data(), msg.size());
			com->SendBob(msg.data(), msg.size());
		}
		
		t.Tick("Verify Union Bound");
	}
	
	void PSI_CA_Hybrid2(int party, int nItems, Communicator *com, Machine *machine, float rate)
	{
		Timer t;
		int fieldSize = 80; int nBytes = 10;
		
		MPCOperations mpc;
		
		RepArithShares sx, sy, sz1, sz2;
		
		std::vector<ZZ_p> Fx;
		AESRNG *rng;
		std::vector<unsigned char> rand;
		
		if(party != Charlie)
		{
			std::vector<uint64_t> data;
			GenerateInput(party, data, nItems, rate);
			rand  = com->getSharedRandomSeed(party, partner);
			rng = new AESRNG(rand.data());
			rng->ComputeFx(Fx, data);
		}
		
		mpc.MPC_ShareData(sx, Alice, Fx, nBytes, nItems, com);
		mpc.MPC_ShareData(sy, Bob, Fx, nBytes, nItems, com);
		
		sx.x.insert(sx.x.end(), std::make_move_iterator(sy.x.begin()), std::make_move_iterator(sy.x.end()));
		sx.y.insert(sx.y.end(), std::make_move_iterator(sy.y.begin()), std::make_move_iterator(sy.y.end()));
		
		t.Tick("Secret share data");
		
// 		std::cout << "Shuffling shares" << std::endl;
		if(party == Charlie)
		{
			mpc.MPC_Sort(sz1, sz2, Fx, sx, nBytes, com);
		}
		else
		{
			std::vector<ZZ_p> temp;
			mpc.MPC_Sort(sz1, sz2, temp, sx, nBytes, com);
		}
		
		
		t.Tick("3pc sort");
		
		RepArithShares psi; psi.x.resize(sz1.x.size()/2); psi.y.resize(sz1.x.size()/2);
		
		for(int idx = 0; idx < sz1.x.size()/2; idx++)
		{
			psi.x[idx] = sz1.x[2*idx] - sz1.x[2*idx + 1];
			psi.y[idx] = sz1.y[2*idx] - sz1.y[2*idx + 1];
		}
		
		mpc.MPC_CheckZero(psi, nBytes, com);
		
		t.Tick("Verify PSI");
		
// 		if(sz1.x.size()/2 == nItems) return;
		
		// Prove union bound by polynomial approach
		int psiSize = psi.x.size();
		int unionSize = 2*nItems - psiSize;
		
		unsigned char *bytes = new unsigned char[nBytes*2*nItems];
		
		// Alice & Bob Compute f(k, x) and send to Charlie
		if(party == Alice || party == Bob)
		{
// 			std::cout << "Union/psi size: " << unionSize << " " << psiSize << " " << (unionSize + psiSize) <<  std::endl;
			
			std::vector<unsigned __int128> coeffs = rng->GetUInt128Array(unionSize);
			
			ZZ_pX f;
			f.SetLength(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				f[idx] = coeffs[idx];
			}
			
			coeffs.resize(0);
			coeffs.shrink_to_fit();
			
			printPolynomial(f, "Random polynomial F(x)", 1);
			
// 			std::cout << "Evaluate the polynomial over " << nItems << " points" << std::endl;
			// Sending data for the union
			ZZ_p *evals = new ZZ_p[nItems];
			evaluatePolynomial(evals, f, Fx.data(), nItems);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, evals, nBytes, nItems);
			com->SendCharlie(bytes, nItems*nBytes);
			
			// Get commitment of f[0] from Charlie
			std::vector<unsigned char> commitment(COMMITMENT_LENGTH);
			com->AwaitCharlie(commitment.data(), commitment.size());
			
			// Send f(x) to Charlie
			ArrayEncoder::ZZpX2ByteArray(bytes, f, nBytes, unionSize);
			com->SendCharlie(bytes, nBytes*unionSize);
			
			// Send rng to Charlie
			com->SendCharlie(rand.data(), rand.size());
			
			std::vector<unsigned char> seed(SEED_LENGTH), msg(nBytes);
			com->AwaitCharlie(seed.data(), seed.size());
			com->AwaitCharlie(msg.data(), msg.size());
			HashCommitment hc;
			assert(hc.Verification(commitment, msg, seed));
			
			delete rng;
			delete [] evals;
		}
		else if(party == Charlie)
		{
			// Charlie count the size of the union and send it to Alice & Bob 
			unsigned char *ptr = bytes;
			
			// Sort array according to f(k,x) and remember the permutation
// 			std::cout << "Sort data based on f(k,x)" << std::endl;
			std::vector<ZZ_p> Fx_copy = Fx;
			std::vector<ZZ> dataFx(2*nItems);
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				conv(dataFx[idx], Fx[idx]);
			}
			
			vector<int> index(2*nItems);
			for (int idx = 0 ; idx != index.size() ; idx++) {
				index[idx] = idx;
			}
			
			// Get the permutation
			std::sort(index.begin(), index.end(), [&](const int& a, const int& b) {return (dataFx[a] < dataFx[b]);});
			
			// Sort
			std::sort(dataFx.begin(), dataFx.end());
			
			// Remove duplicated items
			dataFx.erase(std::unique(dataFx.begin(), dataFx.end()), dataFx.end());
			Fx.resize(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				conv(Fx[idx], dataFx[idx]);
			}
			
			///-------------------------------------Union Interpolation-----------------------------------------------///
			Timer t2;
			
			///----------------Build Tree From F(k1,x)-----------------///
			
			// Build tree and evaluate derivative polynomial
			PolynomialTree *tree = new PolynomialTree();
			buildTree(tree, Fx.data(), unionSize);
			t2.Tick("Build tree");
			ZZ_pX derivativePolynomial;
			derivativePolynomial.SetLength(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				derivativePolynomial[idx] = (tree->val[idx + 1]*(idx + 1));
			}
			
			ZZ_p *derivativeEvals = new ZZ_p[unionSize];
			
			downtreeDivision(derivativeEvals, tree, derivativePolynomial, 0, unionSize);
			
			ZZ_p *evals = new ZZ_p[2*nItems];
			
			///----------------Interpolate from p1(F(k1,x))-----------------///
			
			com->AwaitAlice(bytes, nItems*nBytes);
			com->AwaitBob((bytes + nItems*nBytes), nItems*nBytes);
			
			ptr = bytes;
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				conv(evals[idx], ZZFromBytes(ptr, nBytes));
				ptr += nBytes;
			}
			
// 			std::cout << "Verify (x, p(x))" << std::endl;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(Fx_copy[index[idx - 1]] == Fx_copy[index[idx]])
				{
					assert(evals[index[idx - 1]] == evals[index[idx]]);
				}
			}
			
			ZZ_p *unionEvals = new ZZ_p[unionSize];
			
			int count = 0; 
			int count2 = 0;
			
			unionEvals[count] = evals[index[0]]; count++;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(Fx_copy[index[idx - 1]] != Fx_copy[index[idx]])
				{
					unionEvals[count] = evals[index[idx]];
					count++;
				}
			}
			
			assert(count == unionSize);
			
			reconstruct(unionEvals, derivativeEvals, tree, 0, unionSize);
			
			ZZ_pX f = tree->interp;
			
			printPolynomial(f, "Interpolated polynomial", 1);
			
			t2.Tick("Reconstruction");
			
			// Send Alice & Bob f[0]
			std::vector<unsigned char> seed;
			
			std::vector<unsigned char> msg(nBytes);
			BytesFromZZ(msg.data(), rep(f[0]), nBytes);
			
			HashCommitment hc;
			
			std::vector<unsigned char> commitment = hc.Commit(msg, seed);
			
			com->SendAlice(commitment.data(), commitment.size());
			com->SendBob(commitment.data(), commitment.size());
			
			// Receiver p(x) from Alice & Bob 
			unsigned char *pA = new unsigned char[nBytes*unionSize];
			unsigned char *pB = new unsigned char[nBytes*unionSize];
			unsigned char *pC = new unsigned char[nBytes*psiSize];
			unsigned char *pD = new unsigned char[nBytes*psiSize];
			com->AwaitAlice(pA, nBytes*unionSize);
			com->AwaitBob(pB, nBytes*unionSize);
			
			assert(0 == std::memcmp(pA, pB, nBytes*unionSize));
			
			// Verify that pA = f
			
			ptr = pA;
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				ZZ temp;
				ZZFromBytes(temp, ptr, nBytes);
				assert(temp == rep(f[idx]));
				ptr += nBytes;
			}
			
			// Receive k2 from PA and PB
			std::vector<unsigned char> randAlice(SEED_LENGTH/2), randBob(SEED_LENGTH/2);
			com->AwaitAlice(randAlice.data(), SEED_LENGTH/2);
			com->AwaitBob(randBob.data(), SEED_LENGTH/2);
			
			assert(randAlice == randBob);
			
			delete [] pA;
			delete [] pB;
			delete [] pC;
			delete [] pD;
			delete [] derivativeEvals;
			delete [] unionEvals;
			
			com->SendAlice(seed.data(), seed.size());
			com->SendBob(seed.data(), seed.size());
			
			com->SendAlice(msg.data(), msg.size());
			com->SendBob(msg.data(), msg.size());
		}
		
		t.Tick("Verify Union Bound");
	}
	
	void polynomialTest(int party, int nItems)
	{
// 		SetNumThreads(4);
		Timer t;
		ZZ_pX f;
		
		ZZ n = (ZZ(1) << 20) ;
		
// 		f.SetLength(nItems);
// 		
// 		for(int idx = 0; idx < nItems - 1; idx++)
// 		{
// 			f[idx] = 0;
// 		}
// 		
// 		f[nItems - 1] = 1;
		
		ZZ_p *evals = new ZZ_p[nItems];
		std::vector<ZZ_p> dataFx(nItems);
		
		for(int idx = 0; idx < nItems; idx++)
		{
			random(dataFx[idx]); 
		}
		
		Timer NewRSA;
		// Build tree to get g(x) = (x - x_1)...(x - x_n)
		PolynomialTree *tree = new PolynomialTree();
		buildTree(tree, dataFx.data(), nItems);
		t.Tick("Build tree");
		
		// Compute r(x) = f(x) mod g(x)
		ZZ_pXModulus g(tree->val);
		polynomialMod(f, n, g);
		t.Tick("r(x) = f(x) mod g(x)");
		
		downtreeDivision(evals, tree, f, 0, nItems);
		t.Tick("Evaluate");

		NewRSA.Tick("New RSA");
		std::cout << "evals[0]: " << evals[0] << std::endl;
		
		ZZ_p temp;
		
		for(int idx = 0; idx < nItems; idx++)
		{
			power(temp, dataFx[idx], n);
			if(idx == 0) std::cout << "\nevals'[0]: " << temp << std::endl;
		}
		
		t.Tick("Regular RSA");
		
// 		evaluatePolynomial(evals, f, dataFx.data(), nItems);
// 		pt.Tick("Batched RSA");
	}
	
	void TestRSAMultiEval(int nItems, int fieldSize)
	{
		Timer t;
		ZZ_p temp; 
		random(temp);
		ZZ pp = (ZZ(1) << fieldSize) - 1;;
		for(int idx = 0; idx < 1000; idx++)
		{
			random(temp);
			power(temp, temp, pp);
		}
		t.Tick("RSA");
		
		polynomialTest(party, 1000);
		
		t.Tick("Poly RSA");
	}
	
	int unionCountLowerBound(int party, int nItems, Communicator *com, Machine *machine, bool doPSI)
	{
		SetNumThreads(8);
		
		// Bytes required for each f(x)
		int nBytes = ceil(40 + 2*20)/8;
		
		std::cout << "nItems: " << nItems << "\tnBytes: " << nBytes << "\t" << ZZ_p::ModulusSize() << std::endl;
		
		std::vector<uint64_t> data;
		std::vector<ZZ_p> dataFx;
		
		std::cout << "Generate input" << std::endl;
		GenerateInput(party, data, nItems);
		
		Timer t;
		
		int unionSize;
		
		// Alice & Bob Compute f(k, x) and send to Charlie
		if(party == Alice || party == Bob)
		{
			std::vector<unsigned char> rand  = machine->getSharedRandomSeed(party, partner, com);
			AESRNG *rng = new AESRNG(rand.data());
			rng->ComputeFx(dataFx, data);
			
			t.Tick("AES");
			
			unsigned char *bytes = new unsigned char[nBytes*nItems];
			unsigned char *ptr = bytes;
			
			ArrayEncoder::ZZArray2ByteArray(bytes, dataFx.data(), nBytes, nItems);
			
			std::cout << "Sending Charlie: " << (nBytes*nItems) << " bytes" << std::endl;
			com->SendCharlie(bytes, nBytes*nItems);
			
// 			int unionSize;
			
			com->AwaitCharlie((unsigned char *)(&unionSize), sizeof(int));
			
			std::cout << "Union size from Charlie: " << unionSize << std::endl;
			
			// Sample a random polynomial of degree unionSize
			
			std::cout << "Sample a random polynomial: " << (nBytes*unionSize) << " bytes" << std::endl;
			ZZ_pX f, g;
			
			random(f, unionSize);
			
			unsigned char *polyBytes = new unsigned char[nBytes*unionSize];
			unsigned char *theirBytes = new unsigned char[nBytes*unionSize];
			
			ArrayEncoder::ZZpX2ByteArray(polyBytes, f, nBytes, unionSize);
			
			
			if(party == Alice)
			{
				com->SendEvaluationPartner(polyBytes, nBytes*unionSize);
				com->AwaitEvaluationPartner(theirBytes, nBytes*unionSize);
			}
			else if (party == Bob)
			{
				com->AwaitEvaluationPartner(theirBytes, nBytes*unionSize);
				com->SendEvaluationPartner(polyBytes, nBytes*unionSize);
			}
			
			ArrayEncoder::ByteArray2ZZpX(g, theirBytes, nBytes, unionSize);
			
			delete [] theirBytes;
			
			add(f, f, g);
			
			printPolynomial(f, "Random polynomial", 10);
			
			std::cout << "Evaluate the polynomial over " << nItems << " points" << std::endl;
			// Evaluate f over data points
			ZZ_p *evals = new ZZ_p[nItems];
			evaluatePolynomial(evals, f, dataFx.data(), nItems);
			
			ArrayEncoder::ZZArray2ByteArray(bytes, evals, nBytes, nItems);
			
			com->SendCharlie((unsigned char *)bytes, nItems*nBytes);
			
			// Get commitment of f[0] from Charlie
			std::vector<unsigned char> commitment(COMMITMENT_LENGTH);
			com->AwaitCharlie(commitment.data(), commitment.size());
			
			ArrayEncoder::ZZpX2ByteArray(polyBytes, f, nBytes, unionSize);
			
			// Send p(x) to Charlie
			com->SendCharlie(polyBytes, nBytes*unionSize);
			
			std::vector<unsigned char> seed(SEED_LENGTH), msg(nBytes);
			
			com->AwaitCharlie(seed.data(), seed.size());
			com->AwaitCharlie(msg.data(), msg.size());
			
			HashCommitment hc;
			
			assert(hc.Verification(commitment, msg, seed));
			
			if(doPSI)
			{
				unsigned char *psi = new unsigned char[nBytes*(2*nItems - unionSize)];
				
				com->AwaitCharlie(psi, nBytes*(2*nItems - unionSize));
				
				std::vector<unsigned char> myHash = CryptoUtility::ComputeHash(psi, nBytes*(2*nItems - unionSize));
				
				std::vector<unsigned char> theirHash(myHash.size());
				
				if(party == Alice)
				{
					com->SendEvaluationPartner(myHash.data(), myHash.size());
					com->AwaitEvaluationPartner(theirHash.data(), theirHash.size());
				}
				else if(party == Bob)
				{
					com->AwaitEvaluationPartner(theirHash.data(), theirHash.size());
					com->SendEvaluationPartner(myHash.data(), myHash.size());
				}
				
				assert(myHash == theirHash);
				
				ZZ *psiFx = new ZZ[2*nItems - unionSize];
				
				ArrayEncoder::ByteArray2ZZArray(psiFx, psi, nBytes, (2*nItems - unionSize));
				
				delete [] psi;
			}
		}
		else if(party == Charlie)
		{
			// Charlie count the size of the union and send it to Alice & Bob 
			std::cout << "Receiving data from Alice and Bob: " << (nBytes*2*nItems) << " bytes" << std::endl;
			unsigned char *bytes = new unsigned char[nBytes*2*nItems];
			unsigned char *ptr = bytes;
			
			com->AwaitAlice(bytes, nBytes*nItems);
			com->AwaitBob((bytes + nBytes*nItems), nBytes*nItems);
			
			dataFx.resize(2*nItems);
			std::vector<ZZ> dataFx2(2*nItems);
			
			ArrayEncoder::ByteArray2ZZArray(dataFx2.data(), bytes, nBytes, 2*nItems);
			
			// Sort array according to f(k,x) and remember the permutation
			std::cout << "Sort data based on f(k,x)" << std::endl;
			
			std::vector<ZZ> dataFxBak = dataFx2;
			
			vector<int> index(2*nItems);
			for (int idx = 0 ; idx != index.size() ; idx++) {
				index[idx] = idx;
			}
			
			// Get the permutation
			std::sort(index.begin(), index.end(), [&](const int& a, const int& b) {return (dataFx2[a] < dataFx2[b]);});
			
			// Sort
			std::sort(dataFx2.begin(), dataFx2.end());
			
			// Remove duplicated items
			dataFx2.erase(std::unique(dataFx2.begin(), dataFx2.end()), dataFx2.end());
			
			std::cout << "union size: " << dataFx2.size() << std::endl;
			
			unionSize = dataFx2.size();
			
			dataFx.resize(unionSize);
			for(int idx = 0; idx < unionSize; idx++)
			{
				 conv(dataFx[idx], dataFx2[idx]);
			}
			
			com->SendAlice((unsigned char *)(&unionSize), sizeof(int));
			com->SendBob((unsigned char *)(&unionSize), sizeof(int));
			
			Timer t2;
			std::cout << "Build tree: " << std::endl;
			// Build tree and evaluate derivative polynomial
			PolynomialTree *tree = new PolynomialTree();
			buildTree(tree, dataFx.data(), unionSize);
			t2.Tick("Build tree");
			ZZ_pX derivativePolynomial;
			derivativePolynomial.SetLength(unionSize);
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				derivativePolynomial[idx] = (tree->val[idx + 1]*(idx + 1));
			}
			
			ZZ_p *derivativeEvals = new ZZ_p[unionSize];
			
			downtreeDivision(derivativeEvals, tree, derivativePolynomial, 0, unionSize);
			
			t2.Tick("Downtree division");
			unsigned char *arr = new unsigned char[2*nItems*nBytes];
			
			ZZ_p *evals = new ZZ_p[2*nItems];
			
			com->AwaitAlice(arr, nItems*nBytes);
			com->AwaitBob((arr + nItems*nBytes), nItems*nBytes);
			
			ptr = arr;
			
			for(int idx = 0; idx < 2*nItems; idx++)
			{
				conv(evals[idx], ZZFromBytes(ptr, nBytes));
				ptr += nBytes;
			}
			
			delete [] arr;
			
			std::cout << "Verify (x, p(x))" << std::endl;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(dataFxBak[index[idx - 1]] == dataFxBak[index[idx]])
				{
					assert(evals[index[idx - 1]] == evals[index[idx]]);
				}
			}
			
			ZZ_p *unionEvals = new ZZ_p[unionSize];
			ZZ *intersection = new ZZ[2*nItems - unionSize];
			
			int count = 0; 
			int count2 = 0;
			
			unionEvals[count] = evals[index[0]]; count++;
			
			for(int idx = 1; idx < 2*nItems; idx++)
			{
				if(dataFxBak[index[idx - 1]] != dataFxBak[index[idx]])
				{
					unionEvals[count] = evals[index[idx]];
					count++;
				}
				else
				{
					intersection[count2] = dataFxBak[index[idx]];
					count2++;
				}
			}
			
			assert(count2 == (2*nItems - unionSize));
			
			std::cout << "count: " << count << "\tunion size: " << unionSize << std::endl;
			
			reconstruct(unionEvals, derivativeEvals, tree, 0, unionSize);
			
			ZZ_pX f = tree->interp;
			
			printPolynomial(f, "Interpolated polynomial", 10);
			
			t2.Tick("Reconstruction");
			
			// Send Alice & Bob f[0]
			std::vector<unsigned char> seed;
			
			std::vector<unsigned char> msg(nBytes);
			BytesFromZZ(msg.data(), rep(f[0]), nBytes);
			
			HashCommitment hc;
			
			std::vector<unsigned char> commitment = hc.Commit(msg, seed);
			
			com->SendAlice(commitment.data(), commitment.size());
			com->SendBob(commitment.data(), commitment.size());
			
			// Receiver p(x) from Alice & Bob 
			unsigned char *pA = new unsigned char[nBytes*unionSize];
			unsigned char *pB = new unsigned char[nBytes*unionSize];
			
			com->AwaitAlice(pA, nBytes*unionSize);
			com->AwaitBob(pB, nBytes*unionSize);
			
			assert(0 == std::memcmp(pA, pB, nBytes*unionSize));
			
			// Verify that pA = f
			
			ptr = pA;
			
			for(int idx = 0; idx < unionSize; idx++)
			{
				ZZ temp;
				ZZFromBytes(temp, ptr, nBytes);
				assert(temp == rep(f[idx]));
				ptr += nBytes;
			}
			
			delete [] pA;
			delete [] pB;
			
			com->SendAlice(seed.data(), seed.size());
			com->SendBob(seed.data(), seed.size());
			
			com->SendAlice(msg.data(), msg.size());
			com->SendBob(msg.data(), msg.size());
			
			if(doPSI)
			{
				// Charlie sends f(k,x) in the intersection to Alice and Bob
				int intersectionSize = 2*nItems - unionSize;
				
				unsigned char *psi = new unsigned char[nBytes*intersectionSize];
				ptr = psi;
				
				for(int idx = 0; idx < intersectionSize; idx++)
				{
					BytesFromZZ(ptr, intersection[idx], nBytes);
					ptr += nBytes;
				}
				
				com->SendAlice(psi, nBytes*intersectionSize);
				com->SendBob(psi, nBytes*intersectionSize);
				
				delete [] psi;
			}
		}
			
		t.Tick("");
		
		std::cout << "Communication (MBytes): " << com->sendingSize/1024.0/1024.0 << "\t" << com->receivingSize/1024.0/1024.0 << "\t" << (com->sendingSize + com->receivingSize)/1024.0/1024.0 << std::endl;
		
		return unionSize;
	}

	void GenerateInput(int party, std::vector<uint64_t>& data, int nItems, float rate = 0.5)
	{
		if(party == Alice)
		{
			data.resize(nItems);
			
			for(int idx = 0; idx < nItems; idx++)
			{
				data[idx] = idx + 1;
			}
		}
		else if(party == Bob)
		{
// 			srand (time(NULL));
			int intersectionSize;
			
			if(rate <= 1) intersectionSize = nItems*rate; //rand() % (nItems);
			else intersectionSize = rate;
			
// 			if(rate == 1) intersectionSize -= 64;
			
// 			if(intersectionSize == 0) intersectionSize = 64;
			std::cout << "intersection size: " << intersectionSize << std::endl;
			
			data.resize(nItems);
			
			for(int idx = 0; idx < nItems; idx++)
			{
				data[idx] = (idx + 1 + nItems - intersectionSize);
			}
		}
	}
	
	void sse_trans(uint8_t *inp, uint8_t *out, int nrows, int ncols)
	{
	    #   define INP(x,y) inp[(x)*ncols/8 + (y)/8]
	    #   define OUT(x,y) out[(y)*nrows/8 + (x)/8]
		int rr, cc, i, h;
		union { __m128i x; uint8_t b[16]; } tmp;
		assert(nrows % 8 == 0 && ncols % 8 == 0);

		// Do the main body in 16x8 blocks:
		for (rr = 0; rr <= nrows - 16; rr += 16) {
			for (cc = 0; cc < ncols; cc += 8) {
				for (i = 0; i < 16; ++i)
					tmp.b[i] = INP(rr + i, cc);
				for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
					*(uint16_t*)&OUT(rr,cc+i)= _mm_movemask_epi8(tmp.x);
			}
		}
		if (rr == nrows) return;

		// The remainder is a block of 8x(16n+8) bits (n may be 0).
		//  Do a PAIR of 8x8 blocks in each step:
		for (cc = 0; cc <= ncols - 16; cc += 16) {
			for (i = 0; i < 8; ++i) {
				tmp.b[i] = h = *(uint16_t const*)&INP(rr + i, cc);
				tmp.b[i + 8] = h >> 8;
			}
			for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1)) {
				OUT(rr, cc + i) = h = _mm_movemask_epi8(tmp.x);
				OUT(rr, cc + i + 8) = h >> 8;
			}
		}
		if (cc == ncols) return;

		//  Do the remaining 8x8 block:
		for (i = 0; i < 8; ++i)
			tmp.b[i] = INP(rr + i, cc);
		for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
			OUT(rr, cc + i) = _mm_movemask_epi8(tmp.x);
	}
	
	// Incorrect implementation
	int unionPSICount(int party, int nItems, Communicator *com, Machine *machine, Triplets& tps, bool doPSI)
	{
		return 0;
	}
};

#endif



	





