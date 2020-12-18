#ifndef REPLICATED_RING_SHARES_H__
#define REPLICATED_RING_SHARES_H__

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
class RepRingShares 
{
public:
	/// Constructors
	RepRingShares(){
		x.resize(0);
		y.resize(0);
	}
	
	RepRingShares(const RepRingShares& other):x(other.x), y(other.y)/*, rx(other.rx), ry(other.ry)*/{}	
	RepRingShares(RepRingShares&& other):x(std::move(other.x)), y(std::move(other.y))/*, rx(std::move(other.rx)), ry(std::move(other.ry))*/{}
	
	// Copy/Move Constructors
	RepRingShares(const std::vector<__uint128_t>& x, const std::vector<__uint128_t>& y) : x(x), y(y){}
	RepRingShares(const std::vector<__uint128_t>&& x, const std::vector<__uint128_t>&& y) : x(x), y(y){}
	
	// Copy/Move Assignments
	RepRingShares& operator=(const RepRingShares& other)
	{
		if(this != &other)
		{
			this->x = other.x;
			this->y = other.y;
		}
		return *this;
	}
	
	RepRingShares& operator=(RepRingShares&& other)
	{
		this->x = std::move(other.x);
		this->y = std::move(other.y);
		return *this;
	}
	
	std::vector<__uint128_t> x;
	std::vector<__uint128_t> y;
};

class MPCRingOperations
{
public:
	MPCRingOperations()
	{
		modulo = 1; modulo = (modulo << 80); modulo -= 65;
		mm = 1; mm = (mm << 80); mm -= 1;
		mask = 1; mask = (mask << 40) -1;
		zeroMask = -1; zeroMask -= mm;
	}
	
	MPCRingOperations(__uint128_t modulo){
		this->modulo = modulo;
		mm = 1; mm = (mm << 80); mm -= 1;
		mask = 1; mask = (mask << 40) -1;
		zeroMask = -1; zeroMask -= mm;
	}
	
	__uint128_t ringMult(__uint128_t x, __uint128_t y)
	{
		__uint128_t x1, x2, y1, y2, res;
		x1 = x & mask; x2 = (x >> 40) & mask;
		y1 = y & mask; y2 = (y >> 40) & mask;
		
		res = (65*x2*y2 + x1*y1) + ((x1*y2+x2*y1) << 40);
		
		res = (res >> 80) + ((res >> 80) << 6) + (res & mm);
// 		res = (res >> 80) + ((res >> 80) << 6) + (res & mm);
		while((res >= modulo)) res -= modulo;
		return res;
	}
	
	void modp(__uint128_t& x)
	{
		x = (x >> 80) + ((x >> 80) << 6) + (x & mm);
		while((x >= modulo)) x -= modulo;
	}
	
	__uint128_t ringMult2(__uint128_t x, __uint128_t y)
	{
		__uint128_t x1, x2, y1, y2, res;
		x1 = x & mask; 
		x = x >> 40; 
		x2 = x & mask;
		y1 = y & mask; 
		y = y >> 40; 
		y2 = y & mask;
		res = (65*x2*y2 + x1*y1 + ((x1*y2+x2*y1)*(one << 40) % modulo)) % modulo;
		return res;
	}
	
	void GenerateRandomShares(RepRingShares& s, int nBytes, int length, Communicator *com)
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
			s.x[idx] = *((__uint128_t *)(data + idx)) >> 48;
			s.y[idx] = *((__uint128_t *)(data + idx)) >> 48;
		}
		
// 		std::cout << "GenerateRandomShares: " << s.x[0] << std::endl;
		delete [] data;
#endif
// 		circularForwarding(s.y, s.x, nBytes, length, com);
	}
	
	void MPC_Sort(RepRingShares& sz1, RepRingShares& sz2, std::vector<__uint128_t>& Fx, RepRingShares& shares, int nBytes, Communicator *com)
	{
		Timer t;
		nBytes = 10;
		int length = shares.x.size();
		
		std::vector<__uint128_t> dataShares, macShares;
		std::vector<unsigned char> seed;
		std::vector<__uint128_t> random;
		
		unsigned char *bytes = new unsigned char[nBytes*length];
		
		RepRingShares MacKey; GenerateRandomShares(MacKey, nBytes, 1, com);
		RepRingShares MACs;
		
		RepRingShares vShares, vMACs;
		MACs = repMul(shares, MacKey, nBytes, com);
		
		// Store shares here for verification at the end
		vShares = shares;
		vMACs   = MACs;
		
		// Open f(k,.) to P3 (shares)
		std::vector<__uint128_t> data = openOneParty(shares, nBytes, com, Charlie);
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
			
			for(int idx = 0; idx < length; idx++)
			{
				dataShares[idx] = (random[idx] + dataShares[idx]); modp(dataShares[idx]);
				macShares[idx]  = (random[idx + length] + macShares[idx]); modp(macShares[idx]);
			}
			
			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
			std::shuffle(dataShares.begin(), dataShares.end(), std::mt19937(perm[0]));
			std::shuffle(macShares.begin(), macShares.end(), std::mt19937(perm[0]));
			
			compactData(bytes, dataShares, nBytes);
			com->SendBob(bytes, nBytes*length);
			
			compactData(bytes, macShares, nBytes);
			com->SendBob(bytes, nBytes*length);
		}
		else if(com->party == Bob)
		{
			dataShares.resize(length); macShares.resize(length);
			com->AwaitAlice(bytes, nBytes*length);
			uncompactData(dataShares, bytes, nBytes);
			
			com->AwaitAlice(bytes, nBytes*length);
			uncompactData(macShares, bytes, nBytes);
		}
		else if(com->party == Charlie)
		{
			dataShares = std::move(shares.x);
			macShares  = std::move(MACs.x);
			
			for(int idx = 0; idx < length; idx++)
			{
				dataShares[idx] = (dataShares[idx] + shares.y[idx]); modp(dataShares[idx]);
				macShares[idx]  = (macShares[idx] + MACs.y[idx]); modp(macShares[idx]);
			}
			
			seed = com->getSharedRandomSeed(com->party, Alice);
			AESRNG *rng = new AESRNG(seed.data());
			
			rng->GetRandomRingFx(random, nBytes, 2*length);
			
			for(int idx = 0; idx < length; idx++)
			{
				dataShares[idx] = (modulo + dataShares[idx] - random[idx]); modp(dataShares[idx]);
				macShares[idx]  = (modulo + macShares[idx] - random[idx + length]); modp(macShares[idx]);
			}
// 			t.Tick("re-randomize");
			
			std::vector<uint64_t> perm = rng->GetUInt64Array(1);
			std::vector<int> index(data.size());
			std::iota(index.begin(), index.end(), 0);
			
			std::shuffle(index.begin(), index.end(), std::mt19937(perm[0]));
			
// 			t.Tick("init index");
			std::vector<__uint128_t> temp(data.size());
			for(int idx = 0; idx < data.size(); idx++)
			{
				temp[idx] = data[index[idx]];
			}
			
			data = (temp);
			
			for(int idx = 0; idx < data.size(); idx++)
			{
				temp[idx] = dataShares[index[idx]];
			}
			
			dataShares = (temp);
			
			for(int idx = 0; idx < data.size(); idx++)
			{
				temp[idx] = macShares[index[idx]];
			}
			
			macShares = temp;
			
			
		}
		
// 		t.Tick("Alie - Charlie Shuffle");
		
// 		std::cout << "Charlie computes the permutation for Bob - Charlie shuffling" << std::endl;
		// Bob & Charlie get 2-out-of-2 shares and shuffle and re-randomize
		std::vector<int> index;
		std::vector<int> Z1Index, Z2Index;
		std::vector<__uint128_t> z1, z2, macz1, macz2;
		std::vector<__uint128_t> sortedDataShares, sortedMacShares;
		
		RepRingShares smacz1, smacz2;
		
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
		
// 		t.Tick("Bob - Charlie Shuffle");
		
		RepRingShares sz, smacz;
		
// 		std::cout << "Bob - Charlie Shuffling" << std::endl;
		sz1.x.resize(length); sz1.y.resize(length); smacz1.x.resize(length); smacz1.y.resize(length);
		sz2.x.resize(length - count); sz2.y.resize(length - count); smacz2.x.resize(length - count); smacz2.y.resize(length - count);
		
// 		t.Tick("Resize data");
		
		if(com->party == Alice)
		{
			std::vector<unsigned char> seed1 = com->getSharedRandomSeed(com->party, Bob);
			std::vector<unsigned char> seed2 = com->getSharedRandomSeed(com->party, Charlie);
			
			AESRNG *r1 = new AESRNG(seed1.data());
			AESRNG *r2 = new AESRNG(seed2.data());
			
			std::vector<__uint128_t> r1Vec, r2Vec;
			
			r1->GetRandomRingFx(r1Vec, nBytes, 2*length);
			r2->GetRandomRingFx(r2Vec, nBytes, 2*length);
			
// 			t.Tick("r2->GetRandomRingFx");
			
			for(int idx = 0; idx < length; idx++)
			{
				sz1.x[idx] = r1Vec[idx];
				sz1.y[idx] = r2Vec[idx];
				smacz1.x[idx] = r1Vec[idx + length];
				smacz1.y[idx] = r2Vec[idx + length];
			}
			
// 			t.Tick("re-rand sz/smacz");
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
			
// 			t.Tick("resize data");
			
			std::vector<__uint128_t> rVec;
			
			std::vector<unsigned char> seed = com->getSharedRandomSeed(com->party, Alice);
			AESRNG *r = new AESRNG(seed.data());
			
			r->GetRandomRingFx(rVec, nBytes, 2*length);
// 			t.Tick("r->GetRandomRingFx");
			// Alice - Bob: r1		Alice - Charlie: r2
			std::vector<__uint128_t> temp(length);
			if(com->party == Bob)
			{
				for(int idx = 0; idx < length; idx++)
				{
					sz1.x[idx] = (modulo + z1[idx] - rVec[idx]); modp(sz1.x[idx]);
					sz1.y[idx] = rVec[idx];
					smacz1.x[idx] = (modulo + macz1[idx] - rVec[idx + length]); modp(smacz1.x[idx]);
					smacz1.y[idx] = rVec[idx + length];
				}
				
				compactData(bytes, sz1.x, nBytes);
				com->SendCharlie(bytes, nBytes*length);
				
				compactData(bytes, smacz1.x, nBytes);
				com->SendCharlie(bytes, nBytes*length);
				
// 				com->SendCharlie((unsigned char *)(sz1.x.data()), nBytes*length);
// 				com->SendCharlie((unsigned char *)(smacz1.x.data()), nBytes*length);
				
				// Receive x - r1 - r2
				com->AwaitCharlie(bytes, nBytes*length);
				uncompactData(sz1.x, bytes, nBytes);
				
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					sz1.x[idx] = temp[idx];
// 				}
				
				com->AwaitCharlie(bytes, nBytes*length);
				uncompactData(smacz1.x, bytes, nBytes);
// 				for(int idx = 0; idx < length; idx++)
// 				{
// 					smacz1.x[idx] = temp[idx];
// 				}
			}
			else
			{
				for(int idx = 0; idx < length; idx++)
				{
					sz1.y[idx] = (modulo + z1[idx] - rVec[idx]); modp(sz1.y[idx]);
					sz1.x[idx] = rVec[idx];
					smacz1.y[idx] = (modulo + macz1[idx] - rVec[idx + length]); modp(smacz1.y[idx]);
					smacz1.x[idx] = rVec[idx + length];
				}
				
				// Receive x1 - r1
// 				com->AwaitBob((unsigned char *)(temp.data()), nBytes*length);
				com->AwaitBob(bytes, nBytes*length);
				uncompactData(temp, bytes, nBytes);
				
				for(int idx = 0; idx < length; idx++)
				{
					sz1.y[idx] = (sz1.y[idx] + temp[idx]); modp(sz1.y[idx]);
				}
				
// 				com->AwaitBob((unsigned char *)(temp.data()), nBytes*length);
				com->AwaitBob(bytes, nBytes*length);
				uncompactData(temp, bytes, nBytes);
				
				for(int idx = 0; idx < length; idx++)
				{
					smacz1.y[idx] = (smacz1.y[idx] + temp[idx]); modp(smacz1.y[idx]);
				}
				
				compactData(bytes, sz1.y, nBytes);
				com->SendBob(bytes, nBytes*length);
				
				compactData(bytes, smacz1.y, nBytes);
				com->SendBob(bytes, nBytes*length);
				
// 				com->SendBob((unsigned char *)(sz1.y.data()), nBytes*length);
// 				com->SendBob((unsigned char *)(smacz1.y.data()), nBytes*length);
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
// 		t.Tick("Before verifying MACs");
		
		std::vector<__uint128_t> alpha = openAllParties(MacKey, nBytes, com);
		
// 		std::cout << "Verify MACs" << std::endl;
		
		RepRingShares zeros; zeros.x.resize(2*vShares.x.size()); zeros.y.resize(2*vShares.x.size());
		
		int ref = 0;
		for(int idx = 0; idx < vShares.x.size(); idx++)
		{
			zeros.x[idx] = (ringMult(alpha[0], vShares.x[idx]) + modulo - vMACs.x[idx]); modp(zeros.x[idx]);
			zeros.y[idx] = (ringMult(alpha[0], vShares.y[idx]) + modulo - vMACs.y[idx]); modp(zeros.y[idx]);
		}
		
		ref += vShares.x.size();
		
		for(int idx = 0; idx < sz1.x.size(); idx++)
		{
			zeros.x[idx + ref] = (ringMult(alpha[0], sz1.x[idx]) + modulo - smacz1.x[idx]); modp(zeros.x[idx + ref]);
			zeros.y[idx + ref] = (ringMult(alpha[0], sz1.y[idx]) + modulo - smacz1.y[idx]); modp(zeros.y[idx + ref]);
		}
		
		ref += sz1.x.size();
		
		for(int idx = 0; idx < sz2.x.size(); idx++)
		{
			zeros.x[idx + ref] = (ringMult(alpha[0], sz2.x[idx]) + modulo - smacz2.x[idx]); modp(zeros.x[idx + ref]);
			zeros.y[idx + ref] = (ringMult(alpha[0], sz2.y[idx]) + modulo - smacz2.y[idx]); modp(zeros.y[idx + ref]);
		}
		
// 		t.Tick("Prepare MAC check");
		MPC_CheckZero(zeros, nBytes, com);
		
		delete [] bytes;
		
// 		t.Tick("Check Zeros");
// 		std::cout << "End shuffling" << std::endl;
	}
	
	void MPC_ShareData(RepRingShares& shares, int partyId, const std::vector<__uint128_t>& data, int nBytes, int length, Communicator *com)
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
				r[idx] = data[idx] + modulo - r[idx]; modp(r[idx]);
			}
// 			std::cout << "AAA" << std::endl;
			assert(length == r.size());
			compactData(bytes, r, nBytes);
			com->SendNext(bytes, nBytes*length);
			com->SendPrevious(bytes, nBytes*length);
			
// 			com->SendNext((unsigned char *)(r.data()), nBytes*length);
// 			com->SendPrevious((unsigned char *)(r.data()), nBytes*length);
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
			uncompactData(r, bytes, nBytes);
			
// 			for(int idx = 0; idx < length; idx++)
// 			{
// 				r[idx] = *((__uint128_t *)(bytes + nBytes*idx));
// 			}
		}
		
		delete [] bytes;
// 		std::cout << "AAA" << std::endl;
		for(int idx = 0; idx < length; idx++)
		{
			if(com->party == Alice)
			{
				shares.x[idx] = (shares.x[idx] + r[idx]); modp(shares.x[idx]);
			}
			if(com->party == Bob) 
			{
				shares.y[idx] = (shares.y[idx] + r[idx]); modp(shares.y[idx]);
			}
		}
// 		std::cout << "AAA" << std::endl;
	}
	
	void MPC_CheckZero(const RepRingShares& shares, int nBytes, Communicator *com)
	{
		int length = shares.x.size();
		
		std::vector<unsigned char> seed = com->getCommonSeed();
		AESRNG *rng = new AESRNG(seed.data());
		
		std::vector<__uint128_t> randCoeffs(length);
		Timer t;
		rng->GetRandomRingFx(randCoeffs, nBytes, length);
// 		std::cout << "length: " << length << std::endl;
// 		t.Tick("Checkzero random generation");
		RepRingShares sum; sum.x.resize(1); sum.y.resize(1);
		sum.x[0] = 0; sum.y[0] = 0;
		
		for(int idx = 0; idx < length; idx++)
		{
			sum.x[0] += ringMult(randCoeffs[idx], shares.x[idx]);
			sum.y[0] += ringMult(randCoeffs[idx], shares.y[idx]);
		}
		
		sum.x[0] = (sum.x[0] % modulo);
		sum.y[0] = (sum.y[0] % modulo);
		
// 		t.Tick("Linear combination");
		std::vector<__uint128_t> zeroVals = openAllParties(sum, nBytes, com);
		
		std::cout << "sum: " << zeroVals[0] << std::endl;
		
		if((zeroVals[0] % modulo) == 0) std::cout << "MPC Zero: ********************************" << std::endl;
	}
	
	RepRingShares repMul(const RepRingShares& s1, const RepRingShares& s2, int nBytes, Communicator *com)
	{
		nBytes = 10;
		int length = s1.x.size();
		assert(length >= s2.x.size());
		
		RepRingShares res; 
		res.x.resize(length);
		res.y.resize(length);
		
		Timer t;
		
// 		__uint128_t temp;
		if(length == s2.x.size())
		{
			for(int idx = 0; idx < length; idx++)
			{
// 				random(temp);
				res.x[idx] = ringMult(s2.x[idx], s1.x[idx] + s1.y[idx]) + ringMult(s2.y[idx],s1.x[idx]); modp(res.x[idx]);
// 				res.x[idx] = (res.x[idx] >> 80)*65 + (res.x[idx] & mm);
			}
			
// 			t.Tick("Mult operations");
			circularForwarding(res.y, res.x, nBytes, length, com);
// 			t.Tick("Sending data");
		}
		else if(1 == s2.x.size())
		{
			for(int idx = 0; idx < length; idx++)
			{
// 				random(temp);
				res.x[idx] = ringMult(s2.x[0], s1.x[idx] + s1.y[idx]) + ringMult(s2.y[0], s1.x[idx]); modp(res.x[idx]);
// 				res.x[idx] = (res.x[idx] >> 80)*65 + (res.x[idx] & mm);
			}
// 			t.Tick("Mult operations");
			circularForwarding(res.y, res.x, nBytes, length, com);
// 			t.Tick("Sending data");
		}
		
		return res;
	}
	
	RepRingShares repConstMul(const RepRingShares& s1, const __uint128_t& c, int nBytes, Communicator *com)
	{
		int length = s1.x.size();
		
		RepRingShares res; res.x.resize(length); res.y.resize(length);
		
		for(int idx = 0; idx < length; idx++)
		{
			res.x[idx] = ringMult(c, s1.x[idx]);
			res.y[idx] = ringMult(c, s1.y[idx]);
		}
		
		return res;
	}
	
	std::vector<__uint128_t> openOneParty(RepRingShares& shares, int nBytes, Communicator *com, int partyId)
	{
		int length = shares.x.size();
		
		std::vector<__uint128_t> res;
		nBytes = 10;
		
		unsigned char *bytes = new unsigned char[nBytes*length];
		
		if(com->party == partyId)
		{
			unsigned char *rPrev = bytes;
			
			std::vector<unsigned char> pView, nView;
			
			// Charlie - Alice - Bob
			com->AwaitPrevious(rPrev, nBytes*length);
			pView = CryptoUtility::ComputeHash(rPrev, nBytes*length);
			
			nView.resize(pView.size());
			
			com->AwaitNext(nView.data(), nView.size());
			
			assert(nView == pView);
			
			res.resize(length);
			
			uncompactData(res, rPrev, nBytes);
			
			for(int idx = 0; idx < length; idx++)
			{
				res[idx] = (res[idx] + shares.x[idx] + shares.y[idx]); modp(res[idx]);
			}
			
// 			delete [] rPrev;
		}
		else if(com->party == ((partyId + 1) % 3))
		{
			// Alice - Bob - Charlie
			unsigned char *sPrev = bytes;
			compactData(sPrev, shares.x, nBytes);
			std::vector<unsigned char> pView = CryptoUtility::ComputeHash(sPrev, nBytes*length);
			com->SendPrevious(pView.data(), pView.size());
		}
		else if(com->party == ((partyId + 2) % 3))
		{
			// Bob - Charlie - Alice
			unsigned char *sNext = bytes;
			compactData(sNext, shares.y, nBytes);
			com->SendNext(sNext, nBytes*length);
		}
		
		delete [] bytes;
		
		return res;
	}
	
	std::vector<__uint128_t> openAllParties(RepRingShares& shares, int nBytes, Communicator *com)
	{
		nBytes = 10;
		int length = shares.x.size();
		
		unsigned char *rPrev = new unsigned char[nBytes*length];
// 		unsigned char *rNext = new unsigned char[nBytes*length];
		
		unsigned char *bytes = new unsigned char[nBytes*length];
		
		unsigned char *sPrev = bytes;
		unsigned char *sNext = bytes;
		
		std::vector<unsigned char> sView, rView;
		compactData(sPrev, shares.x, nBytes);
		sView = CryptoUtility::ComputeHash(sPrev, nBytes*length);
		rView.resize(sView.size());
		
		compactData(sNext, shares.y, nBytes);
		
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
		uncompactData(res, rPrev, nBytes);
		
		for(int idx = 0; idx < length; idx++)
		{
// 			res[idx] = *(__uint128_t *)(rPrev + nBytes*idx);
			res[idx] = (res[idx] + shares.x[idx] + shares.y[idx]); modp(res[idx]);
		}
		
// 		delete [] sPrev;
// 		delete [] sNext;
		delete [] rPrev;
// 		delete [] rNext;
		delete [] bytes;
		return res;
	}
	
	void circularForwarding(std::vector<__uint128_t>& recv, const std::vector<__uint128_t>& send, int nBytes, int length, Communicator *com)
	{
		unsigned char *bytes = new unsigned char[nBytes*length];
		
		if(com->party == Alice)
		{
			compactData(bytes, send, nBytes);
			com->SendNext(bytes, nBytes*length);
			
			com->AwaitPrevious(bytes, nBytes*length);
			uncompactData(recv, bytes, nBytes);
		}
		else if(com->party == Bob)
		{
			com->AwaitPrevious(bytes, nBytes*length);
			uncompactData(recv, bytes, nBytes);
			
			compactData(bytes, send, nBytes);
			com->SendNext(bytes, nBytes*length);
		}
		else if(com->party == Charlie)
		{
			com->AwaitPrevious(bytes, nBytes*length);
			uncompactData(recv, bytes, nBytes);
			
			compactData(bytes, send, nBytes);
			com->SendNext(bytes, nBytes*length);
		}
		
		delete [] bytes;
	}
	
	__uint128_t zeroMask;
	__uint128_t modulo;
	__uint128_t mask;
	__uint128_t mm;
	__uint128_t one = 1;
};

#endif
