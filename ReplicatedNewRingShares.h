#ifndef REPLICATED_NEW_RING_SHARES_H__
#define REPLICATED_NEW_RING_SHARES_H__

#include <cassert>
#include <emmintrin.h>
#include "Utility/CryptoUtility.h"
#include "Utility/ISecureRNG.h"
#include "Utility/Communicator.h"

#include <emp-tool/utils/block.h>
#include <emp-tool/utils/prp.h>
#include <emp-tool/utils/prg.h>

#if !defined (__AES__)
    #error "AES-NI instructions not enabled"
#endif

using namespace Utility;

// Replicated shares: A: x1, x3 B: x2, x1 C: x3, x2
class RepNewRingShares 
{
public:
	/// Constructors
	RepNewRingShares(){
		x.resize(0);
		y.resize(0);
	}
	
	RepNewRingShares(const RepNewRingShares& other):x(other.x), y(other.y)/*, rx(other.rx), ry(other.ry)*/{}	
	RepNewRingShares(RepNewRingShares&& other):x(std::move(other.x)), y(std::move(other.y))/*, rx(std::move(other.rx)), ry(std::move(other.ry))*/{}
	
	// Copy/Move Constructors
	RepNewRingShares(const std::vector<__uint128_t>& x, const std::vector<__uint128_t>& y) : x(x), y(y){}
	RepNewRingShares(const std::vector<__uint128_t>&& x, const std::vector<__uint128_t>&& y) : x(x), y(y){}
	
	// Copy/Move Assignments
	RepNewRingShares& operator=(const RepNewRingShares& other)
	{
		if(this != &other)
		{
			this->x = other.x;
			this->y = other.y;
		}
		return *this;
	}
	
	RepNewRingShares& operator=(RepNewRingShares&& other)
	{
		this->x = std::move(other.x);
		this->y = std::move(other.y);
		return *this;
	}
	
	std::vector<__uint128_t> x;
	std::vector<__uint128_t> y;
};

class MPCNewRingOperations
{
public:
	MPCNewRingOperations()
	{
		modulo = 1; modulo = (modulo << 80); modulo -= 65;
		mm = 1; mm = (mm << 80); mm -= 1;
		mask = 1; mask = (mask << 40) -1;
	}
	
	MPCNewRingOperations(__uint128_t modulo){
		this->modulo = modulo;
		mm = 1; mm = (mm << 80); mm -= 1;
		mask = 1; mask = (mask << 40) -1;
	}
	
	__uint128_t ringMult(__uint128_t x, __uint128_t y)
	{
		return x*y;
	}
	
	void GenerateRandomShares(RepNewRingShares& s, int nBytes, int length, Communicator *com)
	{
		s.x.resize(length); s.y.resize(length);
		
#if !defined (__AES__)
// 		for(int idx = 0; idx < length; idx++)
// 		{
// 			random(s.x[idx]);
// 		}
#else
		emp::PRG prg(emp::fix_key);
		emp::block *data = new emp::block[length];
		prg.random_block(data, length);
		for(int idx = 0; idx < length; idx++)
		{
			s.x[idx] = *((__uint128_t *)(data + idx));
		}
		
		delete [] data;
#endif
		circularForwarding(s.y, s.x, nBytes, length, com);
	}
	
	void MPC_Sort(RepNewRingShares& sz1, RepNewRingShares& sz2, std::vector<__uint128_t>& Fx, RepNewRingShares& shares, int nBytes, Communicator *com)
	{
		Timer t;
		int length = shares.x.size();
		
		std::vector<__uint128_t> dataShares, macShares;
		std::vector<unsigned char> seed;
		std::vector<__uint128_t> random;
		
		unsigned char *bytes = new unsigned char[nBytes*length];
		
		RepNewRingShares MacKey; GenerateRandomShares(MacKey, nBytes, 1, com);
		RepNewRingShares MACs;
		
		RepNewRingShares vShares, vMACs;
		MACs = repMul(shares, MacKey, nBytes, com);
		t.Tick("repMul");
		
		// Store shares here for verification at the end
		vShares = shares;
		vMACs   = MACs;
		
		// Open f(k,.) to P3 (shares)
		std::vector<__uint128_t> data = openOneParty(shares, nBytes, com, Charlie);
		t.Tick("OpenOneParty");
		if(com->party == Charlie) Fx = data;
		
// 		std::cout << "Alice - Charlie Shuffling" << std::endl;
		// Alice & Charlie get 2-out-of-2 shares and shuffle and re-randomize
		if(com->party == Alice)
		{
			dataShares = std::move(shares.x);
			macShares  = std::move(MACs.x);
			
			seed = com->getSharedRandomSeed(com->party, Charlie);
			AESRNG *rng = new AESRNG(seed.data());
			
			rng->GetRandomRingFx(random, nBytes, 2*length);
			t.Tick("rng->GetRandomRingFx");
			
			for(int idx = 0; idx < length; idx++)
			{
				dataShares[idx] = (random[idx] + dataShares[idx]);
				macShares[idx]  = (random[idx + length] + macShares[idx]);
			}
			
			t.Tick("re-randomize");
			
			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
			
			t.Tick("Shuffle shares");
			
			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
			
			t.Tick("Shuffle macs");
			
// 			ArrayEncoder::ZZArray2ByteArray(bytes, dataShares.data(), nBytes, length);
			com->SendBob((unsigned char *)(dataShares.data()), nBytes*length);
			
// 			ArrayEncoder::ZZArray2ByteArray(bytes, macShares.data(), nBytes, length);
			com->SendBob((unsigned char *)(macShares.data()), nBytes*length);
			t.Tick("Send data to Bob");
		}
		else if(com->party == Bob)
		{
			dataShares.resize(length); macShares.resize(length);
			
			com->AwaitAlice((unsigned char *)(dataShares.data()), nBytes*length);
// 			ArrayEncoder::ByteArray2ZZArray(dataShares.data(), bytes, nBytes, length);
			
			com->AwaitAlice((unsigned char *)(macShares.data()), nBytes*length);
// 			ArrayEncoder::ByteArray2ZZArray(macShares.data(), bytes, nBytes, length);
			t.Tick("Receive data from Alice");
		}
		else if(com->party == Charlie)
		{
			dataShares = std::move(shares.x);
			macShares  = std::move(MACs.x);
			
			for(int idx = 0; idx < length; idx++)
			{
				dataShares[idx] = (dataShares[idx] + shares.y[idx]);
				macShares[idx]  = (macShares[idx] + MACs.y[idx]);
			}
			
			t.Tick("re-randomize");
			
			seed = com->getSharedRandomSeed(com->party, Alice);
			AESRNG *rng = new AESRNG(seed.data());
			
			rng->GetRandomRingFx(random, nBytes, 2*length);
			t.Tick("rng->GetRandomRingFx");
			
			for(int idx = 0; idx < length; idx++)
			{
				dataShares[idx] = (dataShares[idx] - random[idx]);
				macShares[idx]  = (macShares[idx] - random[idx + length]);
			}
			t.Tick("re-randomize");
			
			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
			
			t.Tick("Shuffle shares");
			
			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
			t.Tick("Shuffle macs");
			
			std::shuffle(data.begin(), data.end(), std::mt19937(perm[0]));
			t.Tick("Shuffle data");
		}
		
// 		t.Tick("Alie - Charlie Shuffle");
		
// 		std::cout << "Charlie computes the permutation for Bob - Charlie shuffling" << std::endl;
		// Bob & Charlie get 2-out-of-2 shares and shuffle and re-randomize
		std::vector<int> index;
		std::vector<int> Z1Index, Z2Index;
		std::vector<__uint128_t> z1, z2, macz1, macz2;
		std::vector<__uint128_t> sortedDataShares, sortedMacShares;
		
		RepNewRingShares smacz1, smacz2;
		
		int count = 0;
		
		if(com->party == Alice)
		{
			com->AwaitCharlie((unsigned char *)(&count), sizeof(int));
		}
		else if(com->party == Bob)
		{
			index.resize(length);
			com->AwaitCharlie((unsigned char *)(index.data()), length*sizeof(int));
			
			com->AwaitCharlie((unsigned char *)(&count), sizeof(int));
			
			Z1Index.resize(count); Z2Index.resize(length - count);
			com->AwaitCharlie((unsigned char *)(Z1Index.data()), count*sizeof(int));
			if(count < length) com->AwaitCharlie((unsigned char *)(Z2Index.data()), (length - count)*sizeof(int));
			
			
		}
		else if(com->party == Charlie)
		{
			index.resize(length);
			std::iota(index.begin(), index.end(), 0);
			std::sort(index.begin(), index.end(), [&](const int& a, const int& b) {return ((data[a]) < (data[b]));});
			
			com->SendBob((unsigned char *)(index.data()), length*sizeof(int));
			
			for(int idx = 0; idx < (length - 1); idx++)
			{
				if(data[index[idx]] == data[index[idx + 1]])
				{
					count+= 2;
					Z1Index.push_back(idx);
					Z1Index.push_back(idx + 1);
					idx += 1;
				}
				else
				{
					Z2Index.push_back(idx);
				}
			}
			
			if(data[index[length - 2]] != data[index[length - 1]]) Z2Index.push_back(length - 1);
			
			assert(Z1Index.size() + Z2Index.size() == length);
			
			com->SendBob((unsigned char *)(&count), sizeof(int));
			com->SendBob((unsigned char *)(Z1Index.data()), count*sizeof(int));
			com->SendBob((unsigned char *)(Z2Index.data()), (length - count)*sizeof(int));
			
			com->SendAlice((unsigned char *)(&count), sizeof(int));
		}
		
		t.Tick("Bob - Charlie Shuffle");
		
		RepNewRingShares sz, smacz;
		
// 		std::cout << "Bob - Charlie Shuffling" << std::endl;
		sz1.x.resize(length); sz1.y.resize(length); smacz1.x.resize(length); smacz1.y.resize(length);
		sz2.x.resize(length - count); sz2.y.resize(length - count); smacz2.x.resize(length - count); smacz2.y.resize(length - count);
		
		t.Tick("Resize data");
		
		if(com->party == Alice)
		{
			std::vector<unsigned char> seed1 = com->getSharedRandomSeed(com->party, Bob);
			std::vector<unsigned char> seed2 = com->getSharedRandomSeed(com->party, Charlie);
			
			AESRNG *r1 = new AESRNG(seed1.data());
			AESRNG *r2 = new AESRNG(seed2.data());
			
			std::vector<__uint128_t> r1Vec, r2Vec;
			
			r1->GetRandomRingFx(r1Vec, nBytes, 2*length);
			r2->GetRandomRingFx(r2Vec, nBytes, 2*length);
			
			t.Tick("r2->GetRandomRingFx");
			
			for(int idx = 0; idx < length; idx++)
			{
				sz1.x[idx] = r1Vec[idx];
				sz1.y[idx] = r2Vec[idx];
				smacz1.x[idx] = r1Vec[idx + length];
				smacz1.y[idx] = r2Vec[idx + length];
			}
			
			t.Tick("re-rand sz/smacz");
		}
		else
		{
			z1.resize(count); macz1.resize(count);
			z2.resize(length - count); macz2.resize(length - count);
			
			for(int idx = 0; idx < Z1Index.size(); idx++)
			{
				z1[idx] = dataShares[index[Z1Index[idx]]];
				macz1[idx] = macShares[index[Z1Index[idx]]];
			}
			
			for(int idx = 0; idx < Z2Index.size(); idx++)
			{
				z2[idx] = dataShares[index[Z2Index[idx]]];
				macz2[idx] = macShares[index[Z2Index[idx]]];
			}
			
			z1.insert(z1.end(), z2.begin(), z2.end());
			macz1.insert(macz1.end(), macz2.begin(), macz2.end());
			
			t.Tick("resize data");
			
			std::vector<__uint128_t> rVec;
			
			std::vector<unsigned char> seed = com->getSharedRandomSeed(com->party, Alice);
			AESRNG *r = new AESRNG(seed.data());
			
			r->GetRandomRingFx(rVec, nBytes, 2*length);
			t.Tick("r->GetRandomRingFx");
			// Alice - Bob: r1		Alice - Charlie: r2
			std::vector<__uint128_t> temp(length);
			if(com->party == Bob)
			{
				for(int idx = 0; idx < length; idx++)
				{
					sz1.x[idx] = (z1[idx] - rVec[idx]);
					sz1.y[idx] = rVec[idx];
					smacz1.x[idx] = (macz1[idx] - rVec[idx + length]);
					smacz1.y[idx] = rVec[idx + length];
				}
				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, sz1.x.data(), nBytes, length);
				com->SendCharlie((unsigned char *)(sz1.x.data()), nBytes*length);
				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, smacz1.x.data(), nBytes, length);
				com->SendCharlie((unsigned char *)(smacz1.x.data()), nBytes*length);
				
				// Receive x - r1 - r2
				com->AwaitCharlie((unsigned char *)(temp.data()), nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
				
				for(int idx = 0; idx < length; idx++)
				{
					sz1.x[idx] = temp[idx];
				}
				
				com->AwaitCharlie((unsigned char *)(temp.data()), nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
				
				for(int idx = 0; idx < length; idx++)
				{
					smacz1.x[idx] = temp[idx];
				}
			}
			else
			{
				for(int idx = 0; idx < length; idx++)
				{
					sz1.y[idx] = (z1[idx] - rVec[idx]);
					sz1.x[idx] = rVec[idx];
					smacz1.y[idx] = (macz1[idx] - rVec[idx + length]);
					smacz1.x[idx] = rVec[idx + length];
				}
				
				// Receive x1 - r1
				com->AwaitBob((unsigned char *)(temp.data()), nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
				
				for(int idx = 0; idx < length; idx++)
				{
					sz1.y[idx] = (sz1.y[idx] + temp[idx]);
				}
				
				com->AwaitBob((unsigned char *)(temp.data()), nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
				
				for(int idx = 0; idx < length; idx++)
				{
					smacz1.y[idx] = (smacz1.y[idx] + temp[idx]);
				}
				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, sz1.y.data(), nBytes, length);
				com->SendBob((unsigned char *)(sz1.y.data()), nBytes*length);
				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, smacz1.y.data(), nBytes, length);
				com->SendBob((unsigned char *)(smacz1.y.data()), nBytes*length);
			}
		}
		
		for(int idx = count; idx < length; idx++)
		{
			sz2.x[idx - count] = sz1.x[idx];
			sz2.y[idx - count] = sz1.y[idx];
			smacz2.x[idx - count] = smacz1.x[idx];
			smacz2.y[idx - count] = smacz1.y[idx];
		}
		
		sz1.x.resize(count); sz1.y.resize(count); smacz1.x.resize(count); smacz1.y.resize(count);
		t.Tick("Before verifying MACs");
		
		std::vector<__uint128_t> alpha = openAllParties(MacKey, nBytes, com);
		
		std::cout << "Verify MACs" << std::endl;
		
		RepNewRingShares zeros; zeros.x.resize(2*vShares.x.size()); zeros.y.resize(2*vShares.x.size());
		
		int ref = 0;
		for(int idx = 0; idx < vShares.x.size(); idx++)
		{
			zeros.x[idx] = (ringMult(alpha[0], vShares.x[idx]) - vMACs.x[idx]);
			zeros.y[idx] = (ringMult(alpha[0], vShares.y[idx]) - vMACs.y[idx]);
		}
		
		ref += vShares.x.size();
		
		for(int idx = 0; idx < sz1.x.size(); idx++)
		{
			zeros.x[idx + ref] = (ringMult(alpha[0], sz1.x[idx]) - smacz1.x[idx]);
			zeros.y[idx + ref] = (ringMult(alpha[0], sz1.y[idx]) - smacz1.y[idx]);
		}
		
		ref += sz1.x.size();
		
		for(int idx = 0; idx < sz2.x.size(); idx++)
		{
			zeros.x[idx + ref] = (ringMult(alpha[0], sz2.x[idx]) - smacz2.x[idx]);
			zeros.y[idx + ref] = (ringMult(alpha[0], sz2.y[idx]) - smacz2.y[idx]);
		}
// 		t.Tick("Prepare MAC check");
		MPC_CheckZero(zeros, nBytes, com);
		
		delete [] bytes;
		
		t.Tick("Check Zeros");
		std::cout << "End shuffling" << std::endl;
	}
	
// 	void MPC_Shuffle(RepNewRingShares& shares, int nBytes, Communicator *com)
// 	{
// 		int length = shares.x.size();
// 		
// 		std::vector<__uint128_t> dataShares, macShares;
// 		std::vector<unsigned char> seed;
// 		std::vector<__uint128_t> random;
// 		
// 		unsigned char *bytes = new unsigned char[nBytes*length];
// 		
// 		RepNewRingShares MacKey; GenerateRandomShares(MacKey, nBytes, 1, com);
// 		RepNewRingShares MACs;
// 		
// 		RepNewRingShares vShares, vMACs;
// 		
// // 		MPC_MAC(MACs, rMACs, shares, MacKey, nBytes, com);
// // 		std::cout << "Compute MAC..." << std::endl;
// 		MACs = repMul(shares, MacKey, nBytes, com);
// 		
// 		// Store shares here for verification at the end
// 		vShares = shares;
// 		vMACs   = MACs;
// 		
// // 		std::cout << "Alice - Charlie Shuffling" << std::endl;
// 		// Alice & Charlie get 2-out-of-2 shares and shuffle and re-randomize
// 		if(com->party == Alice)
// 		{
// 			dataShares = (shares.x);
// 			macShares  = (MACs.x);
// 			
// 			seed = com->getSharedRandomSeed(com->party, Charlie);
// 			AESRNG *rng = new AESRNG(seed.data());
// 			
// 			rng->GetRandomRingFx(random, nBytes, 2*length);
// 			
// 			for(int idx = 0; idx < length; idx++)
// 			{
// 				dataShares[idx] += random[idx];
// 				macShares[idx]  += random[idx + length];
// 			}
// 			
// 			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
// 			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
// 			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
// 			
// 			ArrayEncoder::ZZArray2ByteArray(bytes, dataShares.data(), nBytes, length);
// 			com->SendBob(bytes, nBytes*length);
// 			
// 			ArrayEncoder::ZZArray2ByteArray(bytes, macShares.data(), nBytes, length);
// 			com->SendBob(bytes, nBytes*length);
// 		}
// 		else if(com->party == Bob)
// 		{
// 			dataShares.resize(length); macShares.resize(length);
// 			
// 			com->AwaitAlice(bytes, nBytes*length);
// 			ArrayEncoder::ByteArray2ZZArray(dataShares.data(), bytes, nBytes, length);
// 			
// 			com->AwaitAlice(bytes, nBytes*length);
// 			ArrayEncoder::ByteArray2ZZArray(macShares.data(), bytes, nBytes, length);
// 		}
// 		else if(com->party == Charlie)
// 		{
// 			dataShares = (shares.x);
// 			macShares  = (MACs.x);
// 			
// 			for(int idx = 0; idx < length; idx++)
// 			{
// 				dataShares[idx] += shares.y[idx];
// 				macShares[idx]  += MACs.y[idx];
// 			}
// 			
// 			seed = com->getSharedRandomSeed(com->party, Alice);
// 			AESRNG *rng = new AESRNG(seed.data());
// 			
// 			rng->GetRandomRingFx(random, nBytes, 2*length);
// 			
// 			for(int idx = 0; idx < length; idx++)
// 			{
// 				dataShares[idx] -= random[idx];
// 				macShares[idx]  -= random[idx + length];
// 			}
// 			
// 			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
// 			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
// 			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
// 		}
// 		
// 		
// 		std::cout << "Bob - Charlie Shuffling" << std::endl;
// 		// Alice & Charlie get 2-out-of-2 shares and shuffle and re-randomize
// 		if(com->party == Alice)
// 		{
// 			
// 		}
// 		else if(com->party == Bob)
// 		{
// 			seed = com->getSharedRandomSeed(com->party, Charlie);
// 			AESRNG *rng = new AESRNG(seed.data());
// 			
// 			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
// 			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
// 			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
// 		}
// 		else if(com->party == Charlie)
// 		{
// 			seed = com->getSharedRandomSeed(com->party, Bob);
// 			AESRNG *rng = new AESRNG(seed.data());
// 			
// 			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
// 			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
// 			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
// 		}
// 		
// 		std::cout << "Convert to replicated shares" << std::endl;
// 		if(com->party == Alice)
// 		{
// 			std::vector<unsigned char> seed1 = com->getSharedRandomSeed(com->party, Bob);
// 			std::vector<unsigned char> seed2 = com->getSharedRandomSeed(com->party, Charlie);
// 			
// 			AESRNG *r1 = new AESRNG(seed1.data());
// 			AESRNG *r2 = new AESRNG(seed2.data());
// 			
// 			std::vector<__uint128_t> r1Vec, r2Vec;
// 			
// 			r1->GetRandomRingFx(r1Vec, nBytes, 2*length);
// 			r2->GetRandomRingFx(r2Vec, nBytes, 2*length);
// 			
// 			for(int idx = 0; idx < length; idx++)
// 			{
// 				shares.x[idx] = r1Vec[idx];
// 				shares.y[idx] = r2Vec[idx];
// 				MACs.x[idx] = r1Vec[idx + length];
// 				MACs.y[idx] = r2Vec[idx + length];
// 			}
// 		}
// 		else
// 		{
// 			std::vector<__uint128_t> rVec;
// 			
// 			std::vector<unsigned char> seed = com->getSharedRandomSeed(com->party, Alice);
// 			AESRNG *r = new AESRNG(seed.data());
// 			
// 			r->GetRandomRingFx(rVec, nBytes, 2*length);
// 			
// 			// Alice - Bob: r1		Alice - Charlie: r2
// 			std::vector<__uint128_t> temp(length);
// 			if(com->party == Bob)
// 			{
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					shares.x[idx] = dataShares[idx] - rVec[idx];
// 					shares.y[idx] = rVec[idx];
// 					MACs.x[idx] = macShares[idx] - rVec[idx + length];
// 					MACs.y[idx] = rVec[idx + length];
// 				}
// 				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, shares.x.data(), nBytes, length);
// 				com->SendCharlie(bytes, nBytes*length);
// 				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, MACs.x.data(), nBytes, length);
// 				com->SendCharlie(bytes, nBytes*length);
// 				
// 				// Receive x - r1 - r2
// 				com->AwaitCharlie(bytes, nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
// 				
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					shares.x[idx] = temp[idx];
// 				}
// 				
// 				com->AwaitCharlie(bytes, nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
// 				
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					MACs.x[idx] = temp[idx];
// 				}
// 			}
// 			else
// 			{
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					shares.y[idx] = dataShares[idx] - rVec[idx];
// 					shares.x[idx] = rVec[idx];
// 					MACs.y[idx] = macShares[idx] - rVec[idx + length];
// 					MACs.x[idx] = rVec[idx + length];
// 				}
// 				
// 				// Receive x1 - r1
// 				com->AwaitBob(bytes, nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
// 				
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					shares.y[idx] += temp[idx];
// 				}
// 				
// 				com->AwaitBob(bytes, nBytes*length);
// 				ArrayEncoder::ByteArray2ZZArray(temp.data(), bytes, nBytes, length);
// 				
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					MACs.y[idx] += temp[idx];
// 				}
// 				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, shares.y.data(), nBytes, length);
// 				com->SendBob(bytes, nBytes*length);
// 				
// 				ArrayEncoder::ZZArray2ByteArray(bytes, MACs.y.data(), nBytes, length);
// 				com->SendBob(bytes, nBytes*length);
// 			}
// 		}
// 		
// 		std::vector<__uint128_t> alpha = openAllParties(MacKey, nBytes, com);
// 		
// 		std::cout << "Verify MACs" << std::endl;
// 		
// 		RepNewRingShares zeros; zeros.x.resize(2*vShares.x.size()); zeros.y.resize(2*vShares.x.size());
// 		
// 		int ref = 0;
// 		for(int idx = 0; idx < vShares.x.size(); idx++)
// 		{
// 			zeros.x[idx] = alpha[0]*vShares.x[idx] - vMACs.x[idx];
// 			zeros.y[idx] = alpha[0]*vShares.y[idx] - vMACs.y[idx];
// 		}
// 		
// 		ref += vShares.x.size();
// 		
// 		for(int idx = 0; idx < shares.x.size(); idx++)
// 		{
// 			zeros.x[idx + ref] = alpha[0]*shares.x[idx] - MACs.x[idx];
// 			zeros.y[idx + ref] = alpha[0]*shares.y[idx] - MACs.y[idx];
// 		}
// 		
// 		MPC_CheckZero(zeros, nBytes, com);
// 		
// 		delete [] bytes;
// 		
// 		std::cout << "End shuffling" << std::endl;
// 	}
// 	
// 	void MPC_MAC(RepNewRingShares& MACs, RepNewRingShares& rMACs, RepNewRingShares& shares, const RepNewRingShares& key, int nBytes, Communicator *com)
// 	{
// 		// Get random r
// 		RepNewRingShares r;
// 		GenerateRandomShares(r, nBytes, 1, com);
// 		
// 		int length = shares.x.size();
// 		
// 		// Compute rx, ry
// 		RepNewRingShares rShares = repMul(shares, r, nBytes, com);
// 		// Compute Macs.x and Macs.rx
// 		MACs = repMul(shares, key, nBytes, com);
// 		rMACs = repMul(rShares, key, nBytes, com);
// 		
// 		// Verify multiplications
// 		std::vector<unsigned char> seed = com->getCommonSeed();
// 		AESRNG *rng = new AESRNG(seed.data());
// 		
// 		std::vector<__uint128_t> randCoeffs(2*length);
// 		rng->GetRandomRingFx(randCoeffs, nBytes, 2*length);
// 		
// 		RepNewRingShares sum, rsum;
// 		sum.x.resize(1); sum.y.resize(1); rsum.x.resize(1); rsum.y.resize(1);
// 		for(int idx = 0; idx < length; idx++)
// 		{
// 			sum.x[0] += (randCoeffs[idx]*shares.x[idx] + randCoeffs[length + idx]*MACs.x[idx]);
// 			sum.y[0] += (randCoeffs[idx]*shares.y[idx] + randCoeffs[length + idx]*MACs.y[idx]);
// 			
// 			rsum.x[0] += (randCoeffs[idx]*rShares.x[idx] + randCoeffs[length + idx]*rMACs.x[idx]);
// 			rsum.y[0] += (randCoeffs[idx]*rShares.y[idx] + randCoeffs[length + idx]*rMACs.y[idx]);
// 		}
// 		
// 		std::vector<__uint128_t> rval = openAllParties(r, nBytes, com);
// 		
// 		sum.x[0] = rval[0]*sum.x[0] - rsum.x[0];
// 		sum.y[0] = rval[0]*sum.y[0] - rsum.y[0];
// 		
// 		MPC_CheckZero(sum, nBytes, com);
// 	}
// 	
	void MPC_ShareData(RepNewRingShares& shares, int partyId, const std::vector<__uint128_t>& data, int nBytes, int length, Communicator *com)
	{
		GenerateRandomShares(shares, nBytes, length, com);
		
// 		std::cout << "shares[0]: " << shares.x[0] << " " << shares.y[0] << std::endl;
		
		std::vector<__uint128_t> r = openOneParty(shares, nBytes, com, partyId);
		
// 		std::cout << "AAA" << std::endl;
		
		unsigned char *bytes = new unsigned char[nBytes*length];
		
		if(com->party == partyId)
		{
			for(int idx = 0; idx < length; idx++)
			{
				if(data[idx] >= r[idx]) 
				{
					r[idx] = data[idx] - r[idx];
				}
				else
				{
					r[idx] = data[idx] - r[idx];
				}
			}
// 			std::cout << "AAA" << std::endl;
			assert(length == r.size());
			
// 			ArrayEncoder::ZZArray2ByteArray(bytes, r.data(), nBytes, length);
			com->SendNext((unsigned char *)(r.data()), nBytes*length);
			com->SendPrevious((unsigned char *)(r.data()), nBytes*length);
// 			std::cout << "AAA" << std::endl;
		}
		else
		{
			assert(r.size() == 0);
			r.resize(length);
			
			std::vector<unsigned char> myView;
			std::vector<unsigned char> theirView;
			
			if(com->party == ((partyId + 1) % 3))
			{
				com->AwaitPrevious(bytes, nBytes*length);
				
				myView = CryptoUtility::ComputeHash(bytes, nBytes*length);
				theirView.resize(myView.size());
				
				com->SendNext(myView.data(), myView.size());
				com->AwaitNext(theirView.data(), theirView.size());
			}
			else if(com->party == ((partyId + 2) % 3))
			{
				com->AwaitNext(bytes, nBytes*length);
				
				myView = CryptoUtility::ComputeHash(bytes, nBytes*length);
				theirView.resize(myView.size());
				
				com->AwaitPrevious(theirView.data(), theirView.size());
				com->SendPrevious(myView.data(), myView.size());
			}
// 			std::cout << "AAA" << std::endl;
			// Comparing view
			assert(myView == theirView);
			
// 			ArrayEncoder::ByteArray2ZZArray(r.data(), bytes, nBytes, length);
			for(int idx = 0; idx < length; idx++)
			{
				r[idx] = *((__uint128_t *)(bytes + nBytes*idx));
			}
		}
		
		delete [] bytes;
// 		std::cout << "AAA" << std::endl;
		for(int idx = 0; idx < length; idx++)
		{
			if(com->party == Alice) shares.x[idx] = (shares.x[idx] + r[idx]);
			if(com->party == Bob)   shares.y[idx] = (shares.y[idx] + r[idx]);
		}
// 		std::cout << "AAA" << std::endl;
	}
	
	void MPC_CheckZero(const RepNewRingShares& shares, int nBytes, Communicator *com)
	{
		int length = shares.x.size();
		
		std::vector<unsigned char> seed = com->getCommonSeed();
		AESRNG *rng = new AESRNG(seed.data());
		
		std::vector<__uint128_t> randCoeffs(length);
		Timer t;
		rng->GetRandomRingFx(randCoeffs, nBytes, length);
// 		std::cout << "length: " << length << std::endl;
// 		t.Tick("Checkzero random generation");
		RepNewRingShares sum; sum.x.resize(1); sum.y.resize(1);
		sum.x[0] = 0; sum.y[0] = 0;
		
		for(int idx = 0; idx < length; idx++)
		{
			sum.x[0] += ringMult(randCoeffs[idx], shares.x[idx]);
			sum.y[0] += ringMult(randCoeffs[idx], shares.y[idx]);
		}
		
		sum.x[0] = (sum.x[0]);
		sum.y[0] = (sum.y[0]);
		
// 		t.Tick("Linear combination");
		std::vector<__uint128_t> zeroVals = openAllParties(sum, nBytes, com);
		
		std::cout << "sum: " << zeroVals[0] << std::endl;
		
		if((zeroVals[0]) == 0) std::cout << "MPC Zero: ********************************" << std::endl;
	}
	
	RepNewRingShares repMul(const RepNewRingShares& s1, const RepNewRingShares& s2, int nBytes, Communicator *com)
	{
		int length = s1.x.size();
		assert(length >= s2.x.size());
		
		RepNewRingShares res; 
		res.x.resize(length);
		res.y.resize(length);
		
		Timer t;
		
// 		__uint128_t temp;
		if(length == s2.x.size())
		{
			for(int idx = 0; idx < length; idx++)
			{
// 				random(temp);
				res.x[idx] = (ringMult(s2.x[idx], s1.x[idx]) + ringMult(s2.x[idx], s1.y[idx]) + ringMult(s2.y[idx],s1.x[idx]));
			}
			
			t.Tick("Mult operations");
			circularForwarding(res.y, res.x, nBytes, length, com);
			t.Tick("Sending data");
		}
		else if(1 == s2.x.size())
		{
			for(int idx = 0; idx < length; idx++)
			{
// 				random(temp);
				res.x[idx] = (ringMult(s2.x[0], s1.x[idx]) + ringMult(s2.x[0], s1.y[idx]) + ringMult(s2.y[0], s1.x[idx]));
			}
			t.Tick("Mult operations");
			circularForwarding(res.y, res.x, nBytes, length, com);
			t.Tick("Sending data");
		}
		
		return res;
	}
	
	RepNewRingShares repConstMul(const RepNewRingShares& s1, const __uint128_t& c, int nBytes, Communicator *com)
	{
		int length = s1.x.size();
		
		RepNewRingShares res; res.x.resize(length); res.y.resize(length);
		
		for(int idx = 0; idx < length; idx++)
		{
			res.x[idx] = ringMult(c, s1.x[idx]);
			res.y[idx] = ringMult(c, s1.y[idx]);
		}
		
		return res;
	}
	
	std::vector<__uint128_t> openOneParty(RepNewRingShares& shares, int nBytes, Communicator *com, int partyId)
	{
		int length = shares.x.size();
		
		std::vector<__uint128_t> res;
		
		if(com->party == partyId)
		{
			unsigned char *rPrev = new unsigned char[nBytes*length];
			
			std::vector<unsigned char> pView, nView;
			
			// Charlie - Alice - Bob
			com->AwaitPrevious(rPrev, nBytes*length);
			pView = CryptoUtility::ComputeHash(rPrev, nBytes*length);
			
			nView.resize(pView.size());
			
			com->AwaitNext(nView.data(), nView.size());
			
			assert(nView == pView);
			
			res.resize(length);
// 			ArrayEncoder::ByteArray2ZZArray(res.data(), rPrev, nBytes, length);
			
			for(int idx = 0; idx < length; idx++)
			{
				res[idx] = *((__uint128_t *)(rPrev + nBytes*idx));
				res[idx] = (res[idx] + shares.x[idx] + shares.y[idx]);
			}
			
			delete [] rPrev;
		}
		else if(com->party == ((partyId + 1) % 3))
		{
			// Alice - Bob - Charlie
			unsigned char *sPrev = (unsigned char *)(shares.x.data());;
// 			ArrayEncoder::ZZArray2ByteArray(sPrev, shares.x.data(), nBytes, length);
			std::vector<unsigned char> pView = CryptoUtility::ComputeHash(sPrev, nBytes*length);
			com->SendPrevious(pView.data(), pView.size());
// 			delete [] sPrev;
		}
		else if(com->party == ((partyId + 2) % 3))
		{
			// Bob - Charlie - Alice
			unsigned char *sNext = (unsigned char *)(shares.y.data());
// 			ArrayEncoder::ZZArray2ByteArray(sNext, shares.y.data(), nBytes, length);
			com->SendNext(sNext, nBytes*length);
// 			delete [] sNext;
		}
		
		return res;
	}
	
	std::vector<__uint128_t> openAllParties(RepNewRingShares& shares, int nBytes, Communicator *com)
	{
		int length = shares.x.size();
		
		unsigned char *sPrev = (unsigned char *)(shares.x.data());
		unsigned char *sNext = (unsigned char *)(shares.y.data());
		
		unsigned char *rPrev = new unsigned char[nBytes*length];
		unsigned char *rNext = new unsigned char[nBytes*length];
		
// 		ArrayEncoder::ZZArray2ByteArray(sPrev, shares.x.data(), nBytes, length);
// 		ArrayEncoder::ZZArray2ByteArray(sNext, shares.y.data(), nBytes, length);
		
		std::vector<unsigned char> sView, rView;
		sView = CryptoUtility::ComputeHash(sPrev, nBytes*length);
		rView.resize(sView.size());
		
		if(com->party == Alice)
		{
			// Charlie - Alice - Bob
			com->SendNext(sNext, nBytes*length);
			com->SendPrevious(sView.data(), sView.size());
			com->AwaitNext(rView.data(), rView.size());
			com->AwaitPrevious(rPrev, nBytes*length);
		}
		else if(com->party == Bob)
		{
			// Alice - Bob - Charlie
			com->AwaitPrevious(rPrev, nBytes*length);
			com->SendNext(sNext, nBytes*length);
			com->SendPrevious(sView.data(), sView.size());
			com->AwaitNext(rView.data(), rView.size());
		}
		else if(com->party == Charlie)
		{
			// Bob - Charlie - Alice
			com->AwaitNext(rView.data(), rView.size());
			com->AwaitPrevious(rPrev, nBytes*length);
			com->SendNext(sNext, nBytes*length);
			com->SendPrevious(sView.data(), sView.size());
		}
		
		sView = CryptoUtility::ComputeHash(rPrev, nBytes*length);
		assert(sView == rView);
		
		std::vector<__uint128_t> res(length);
		
// 		ArrayEncoder::ByteArray2ZZArray(res.data(), rPrev, nBytes, length);
		
		for(int idx = 0; idx < length; idx++)
		{
			res[idx] = *(__uint128_t *)(rPrev + nBytes*idx);
			res[idx] = (res[idx] + shares.x[idx] + shares.y[idx]);
		}
		
// 		delete [] sPrev;
// 		delete [] sNext;
		delete [] rPrev;
		delete [] rNext;
		
		return res;
	}
	
	void circularForwarding(std::vector<__uint128_t>& recv, const std::vector<__uint128_t>& send, int nBytes, int length, Communicator *com)
	{
// 		unsigned char *sendData = new unsigned char[nBytes*length];
// 		unsigned char *recvData = new unsigned char[nBytes*length];
		
// 		ArrayEncoder::ZZArray2ByteArray(sendData, send.data(), nBytes, length);
		
		if(com->party == Alice)
		{
			com->SendNext((uint8_t *)(send.data()), nBytes*length);
			com->AwaitPrevious((uint8_t *)(recv.data()), nBytes*length);
		}
		else if(com->party == Bob)
		{
			com->AwaitPrevious((uint8_t *)(recv.data()), nBytes*length);
			com->SendNext((uint8_t *)(send.data()), nBytes*length);
		}
		else if(com->party == Charlie)
		{
			com->AwaitPrevious((uint8_t *)(recv.data()), nBytes*length);
			com->SendNext((uint8_t *)(send.data()), nBytes*length);
		}
		
// 		ArrayEncoder::ByteArray2ZZArray(recv.data(), recvData, nBytes, length);
		
// 		delete [] sendData;
// 		delete [] recvData;
	}
	
	__uint128_t modulo;
	__uint128_t mask;
	__uint128_t mm;
	__uint128_t one = 1;
};

#endif
