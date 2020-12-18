#ifndef SHARES_H__
#define SHARES_H__

#include <cassert>
#include <emmintrin.h>
#include "Utility/CryptoUtility.h"
#include "Utility/ISecureRNG.h"
#include "Utility/Communicator.h"

#include <NTL/ZZ.h>

using namespace Utility;

class Shares 
{
public:
	// P1 (s_1, s_3^s_1) P2(s_2, s_1^s_2) P3(s_3, s_2^s_3)
	Shares(){
		s.resize(0);
		t.resize(0);
	}
	
	Shares(const Shares& other):s(other.s), t(other.t){}
	
	Shares(Shares&& other):s(std::move(other.s)), t(std::move(other.t)){}
	
	Shares(const std::vector<uint64_t>& s, const std::vector<uint64_t>& t) : s(s), t(t){}
	
	Shares(std::vector<uint64_t>&& s, std::vector<uint64_t>&& t) : s(std::move(s)), t(std::move(t)){}
	
	Shares& operator=(const Shares& other)
	{
		if(this != &other)
		{
			this->s = other.s;
			this->t = other.t;
		}
		return *this;
	}
	
	Shares& operator=(Shares&& other)
	{
		this->s = std::move(other.s);
		this->t = std::move(other.t);
		
		return *this;
	}
	
	void GenerateRandomShares(const uint64_t numTriples, Communicator *com)
	{
		// Generate random key for the PRNG
		std::vector<unsigned char> mySeed = com->getCommonSeed();
		
		uint64_t mySeedVal = ((uint64_t *)mySeed.data())[0];
		uint64_t theirSeedVal;
		
		// Forward key to party (i + 1)
		if(com->party == Alice)
		{
// 			com->SendNext(mySeed.data(), mySeed.size());
// 			com->AwaitPrevious(theirSeed.data(), theirSeed.size());
			theirSeedVal = mySeedVal + 2;
		}
		else 
		{
			if(com->party == Bob)
			{
				theirSeedVal = mySeedVal;
				mySeedVal += 1;
			}
			else
			{
				theirSeedVal = mySeedVal + 1;
				mySeedVal += 2;
			}
// 			com->AwaitPrevious(theirSeed.data(), theirSeed.size());
// 			com->SendNext(mySeed.data(), mySeed.size());
		}
		
		int length = (numTriples + 63)/64;
		
		s.resize(length); t.resize(length);
		
		
		std::mt19937_64 myGen(mySeedVal);
		std::mt19937_64 theirGen(theirSeedVal);
		std::uniform_int_distribution<unsigned long long> dis;
		
		for(int idx = 0; idx < length; idx++)
		{
			s[idx] = dis(myGen);
			t[idx] = dis(theirGen)^s[idx];
		}
	}
	
	void SecretShare(const std::vector<uint64_t>& v, int numTriples, Communicator *com, const int id)
	{
		this->GenerateRandomShares(numTriples, com);
		
		std::vector<uint64_t> ra = this->OpenOneParty(com, id);		
		if(com->party == id)
		{
			assert(ra.size() == v.size());
			
			// Compute a \xor v
			for(int idx = 0; idx < ra.size(); idx++)
			{
				ra[idx] ^= v[idx];
			}
			
			com->SendNext((unsigned char *)(ra.data()), ra.size()*sizeof(uint64_t));
			com->SendPrevious((unsigned char *)(ra.data()), ra.size()*sizeof(uint64_t));
		}
		else
		{
			assert(ra.size() == 0);
			ra.resize(this->s.size());
			
			std::vector<unsigned char> myView;
			std::vector<unsigned char> theirView;
			
			if(com->party == ((id + 1) % 3))
			{
				com->AwaitPrevious((unsigned char *)(ra.data()), ra.size()*sizeof(uint64_t));
				myView = ArrayEncoder::Hash(ra);
				theirView.resize(myView.size());
				
				com->SendNext(myView.data(), myView.size());
				com->AwaitNext(theirView.data(), theirView.size());
			}
			else if(com->party == ((id + 2) % 3))
			{
				com->AwaitNext((unsigned char *)(ra.data()), ra.size()*sizeof(uint64_t));
				myView = ArrayEncoder::Hash(ra);
				theirView.resize(myView.size());
				
				com->AwaitPrevious(theirView.data(), theirView.size());
				com->SendPrevious(myView.data(), myView.size());
			}
			
			// Comparing view
			assert(myView == theirView);
		}
		
		for(int idx = 0; idx < this->s.size(); idx++)
		{
			this->s[idx] ^= ra[idx];
		}
	}
	
	
	std::vector<uint64_t> OpenOneParty(Communicator *com, int id)
	{
		std::vector<uint64_t> t_prev, t_next, res;
		std::vector<unsigned char> mv, tv;
		
		if(com->party == id)
		{
			t_prev.resize(t.size());
			t_next.resize(t.size());
			
			// Charlie - Alice - Bob
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
			
			std::vector<uint64_t> temp(t.size());
			
			for(int idx = 0; idx < t.size(); idx++)
			{
				temp[idx] = t[idx]^t_prev[idx];
			}
			
			mv = CryptoUtility::ComputeHash((unsigned char *)(temp.data()), temp.size()*sizeof(uint64_t));
			tv.resize(mv.size());
			
			com->AwaitNext((unsigned char *)(tv.data()), tv.size());
			assert(mv == tv);
		
			res.resize(t.size());
		
			for(int idx = 0; idx < t.size(); idx++)
			{
				res[idx] = s[idx]^t_prev[idx];
			}
		}
		else if(com->party == ((id + 1) % 3))
		{
			// Alice - Bob - Charlie
			std::vector<unsigned char> mv = CryptoUtility::ComputeHash((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			
			com->SendPrevious(mv.data(), mv.size());
// 			com->SendPrevious((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
		}
		else if(com->party == ((id + 2) % 3))
		{
			// Bob - Charlie - Alice
		  
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
		}
		
		return res;
	}
	
	std::vector<uint64_t> OpenAllParty(Communicator *com) const
	{
		std::vector<uint64_t> t_prev, t_next;
		t_prev.resize(t.size());
		t_next.resize(t.size());
		
		std::vector<unsigned char> mv, tv;
		mv = CryptoUtility::ComputeHash((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
		tv.resize(mv.size());
		
		if(com->party == Alice)
		{
			// Charlie - Alice - Bob
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->SendPrevious(mv.data(), mv.size());
			com->AwaitNext(tv.data(), tv.size());
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
		}
		else if(com->party == Bob)
		{
			// Alice - Bob - Charlie
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->SendPrevious(mv.data(), mv.size());
			com->AwaitNext(tv.data(), tv.size());
		}
		else if(com->party == Charlie)
		{
			// Bob - Charlie - Alice
			com->AwaitNext(tv.data(), tv.size());
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->SendPrevious(mv.data(), mv.size());
		}
		
		std::vector<uint64_t> res(t.size());
		std::vector<uint64_t> temp(t.size());
		for(int idx = 0; idx < t.size(); idx++) temp[idx] = t[idx]^t_prev[idx];
		mv = CryptoUtility::ComputeHash((unsigned char *)(temp.data()), temp.size()*sizeof(uint64_t));
		
		assert(mv == tv);
		
		for(int idx = 0; idx < t.size(); idx++)
		{
// 			assert(t[idx] == t_prev[idx]^t_next[idx]);
			res[idx] = s[idx]^t_prev[idx];
		}
		
		return res;
	}
	
	std::vector<uint64_t> OpenAllParty1(Communicator *com) const
	{
		std::vector<uint64_t> t_prev, t_next;
		t_prev.resize(t.size());
		t_next.resize(t.size());
		
		if(com->party == Alice)
		{
			// Charlie - Alice - Bob
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->SendPrevious((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->AwaitNext((unsigned char *)(t_next.data()), t.size()*sizeof(uint64_t));
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
		}
		else if(com->party == Bob)
		{
			// Alice - Bob - Charlie
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->SendPrevious((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->AwaitNext((unsigned char *)(t_next.data()), t.size()*sizeof(uint64_t));
		}
		else if(com->party == Charlie)
		{
			// Bob - Charlie - Alice
			com->AwaitNext((unsigned char *)(t_next.data()), t.size()*sizeof(uint64_t));
			com->AwaitPrevious((unsigned char *)(t_prev.data()), t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
			com->SendPrevious((unsigned char *)(t.data()), t.size()*sizeof(uint64_t));
		}
		
		std::vector<uint64_t> res(t.size());
		
		for(int idx = 0; idx < t.size(); idx++)
		{
			assert(t[idx] == t_prev[idx]^t_next[idx]);
			res[idx] = s[idx]^t_prev[idx];
		}
		
		return res;
	}
	
	void CompareView(const std::vector<uint64_t>& view, Communicator *com)
	{
		std::vector<unsigned char> myView = ArrayEncoder::Hash(view);
		std::vector<unsigned char> theirView(myView.size());
		
		if(com->party == Alice)
		{
			com->SendNext(myView.data(), myView.size());
			com->AwaitPrevious(theirView.data(), theirView.size());
		}
		else
		{
			com->AwaitPrevious(theirView.data(), theirView.size());
			com->SendNext(myView.data(), myView.size());
		}
		
		assert(myView == theirView);
	}
	
	void CompareView2(const std::vector<uint64_t>& view, Communicator *com)
	{
		std::vector<unsigned char> myView = ArrayEncoder::Hash(this->s);
		std::vector<unsigned char> myView2 = ArrayEncoder::Hash(this->t);
		std::vector<unsigned char> mv = ArrayEncoder::Hash(view);
		
		std::vector<unsigned char> theirView(myView.size() + mv.size());
		
		myView2.insert(myView2.end(), mv.begin(), mv.end());
		myView.insert(myView.end(), mv.begin(), mv.end());
		
		if(com->party == Alice)
		{
			com->SendNext(myView2.data(), myView2.size());
			com->AwaitPrevious(theirView.data(), theirView.size());
		}
		else
		{
			com->AwaitPrevious(theirView.data(), theirView.size());
			com->SendNext(myView2.data(), myView2.size());
		}
		
		assert(myView == theirView);
	}
	
	void TripleVerificationWithOpening(const Shares& a, const Shares& b, const Shares& c, Communicator *com)
	{
		int length = a.s.size();
		
		Shares z = a;
		z.s.insert(z.s.end(), b.s.begin(), b.s.end());
		z.s.insert(z.s.end(), c.s.begin(), c.s.end());
		z.t.insert(z.t.end(), b.t.begin(), b.t.end());
		z.t.insert(z.t.end(), c.t.begin(), c.t.end());
		
		std::vector<uint64_t> rz = z.OpenAllParty(com);
		
		for(int idx = 0; idx < length; idx++)
		{
			assert(rz[idx + 2*length] == (rz[idx + length] & rz[idx]));
		}
		
// 		std::vector<uint64_t> ra = a.OpenAllParty(com);
// 		std::vector<uint64_t> rb = b.OpenAllParty(com);
// 		std::vector<uint64_t> rc = c.OpenAllParty(com);
// 		
// 		for(int idx = 0; idx < ra.size(); idx++)
// 		{
// 			assert(rc[idx] == (rb[idx] & ra[idx]));
// 		}
	}
	
	void TripleVerificationWithoutOpening(const Shares& a, const Shares& b, const Shares& c, const Shares& x, const Shares& y, const Shares& z, Communicator *com)
	{
		Shares rho = ComputeXOR(x, a);
		Shares sigma = ComputeXOR(y, b);
		int length = rho.s.size();
		
		rho.s.insert(rho.s.end(), (sigma.s.begin()), (sigma.s.end()));
		rho.t.insert(rho.t.end(), (sigma.t.begin()), (sigma.t.end()));
		
// 		std::cout << "Reveal rho" << std::endl;
		std::vector<uint64_t> reveal_rho = rho.OpenAllParty(com);
		
		std::vector<uint64_t> reveal_sigma, view;
		view = reveal_rho;
		reveal_sigma.insert(reveal_sigma.end(), std::make_move_iterator(reveal_rho.begin() + length), std::make_move_iterator(reveal_rho.end()));
		reveal_rho.resize(length);
		
		Shares res = ComputeXOR(ComputeXOR(z, c), ComputeXOR(ComputeConstAND(a, reveal_sigma), ComputeConstAND(b, reveal_rho)));
		
		for(int idx = 0; idx < res.s.size(); idx++)
		{
			reveal_rho[idx] &= reveal_sigma[idx];
		}
		
		res.ComputeConstXOR(reveal_rho);
		
		res.CompareView2(view, com);
	}
	
	static Shares ComputeXOR(const Shares& x, const Shares& y)
	{
// 		Timer ttt;
		Shares ret;
		
		ret.s.resize(x.s.size());
		ret.t.resize(x.t.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			ret.s[idx] = x.s[idx] ^ y.s[idx];
			ret.t[idx] = x.t[idx] ^ y.t[idx];
		}
		
		return ret;
	}
	
	static Shares ComputeXOR2(const Shares& x, const Shares& y)
	{
		Shares ret;
		
		ret.s.resize(x.s.size());
		ret.t.resize(x.t.size());
		
		for(int idx = 0; idx < x.s.size()/2; idx += 2)
		{
			ret.s[2+idx] = x.s[2*idx] ^ y.s[2*idx];
			ret.s[2+idx+1] = x.s[2*idx+1] ^ y.s[2*idx+1];
			ret.t[2*idx] = x.t[2*idx] ^ y.t[2*idx];
			ret.t[2*idx+1] = x.t[2*idx+1] ^ y.t[2*idx+1];
		}
		
		return ret;
	}
	
	static Shares ComputeXOR(Shares&& x, Shares&& y)
	{
		Shares ret;
		
		ret.s.resize(x.s.size());
		ret.t.resize(x.t.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			ret.s[idx] = x.s[idx] ^ y.s[idx];
			ret.t[idx] = x.t[idx] ^ y.t[idx];
		}
		
		return ret;
	}
	
	void ComputeXOR(const Shares& x)
	{
		assert(this->s.size() == x.s.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			this->s[idx] ^= x.s[idx];
			this->t[idx] ^= x.t[idx];
		}
	}
	
	void ComputeXOR(Shares&& x)
	{
		assert(this->s.size() == x.s.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			this->s[idx] ^= x.s[idx];
			this->t[idx] ^= x.t[idx];
		}
	}
	
	static Shares ComputeConstXOR(const Shares& x, std::vector<uint64_t> sigma)
	{
		assert(x.s.size() == sigma.size());
		
		Shares ret;
		
		ret.s.resize(x.s.size());
		
		for(int idx = 0; idx < sigma.size(); idx++)
		{
			ret.s[idx] = x.s[idx] ^ sigma[idx];
		}
		
		ret.t = x.t;
		
		return ret;
	}
	
	static Shares ComputeConstXOR(const Shares& x, uint64_t sigma)
	{
// 		assert(x.s.size() == sigma.size());
		
		Shares ret;
		
		ret.s.resize(x.s.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			ret.s[idx] = x.s[idx] ^ sigma;
		}
		
		ret.t = x.t;
		
		return ret;
	}
	
	void ComputeConstXOR(const std::vector<uint64_t>& sigma)
	{
		assert(this->s.size() == sigma.size());
		
		for(int idx = 0; idx < sigma.size(); idx++)
		{
			this->s[idx] ^= sigma[idx];
		}
	}
	
	void ComputeConstXOR(std::vector<uint64_t>&& sigma)
	{
		assert(this->s.size() == sigma.size());
		
		for(int idx = 0; idx < sigma.size(); idx++)
		{
			this->s[idx] ^= sigma[idx];
		}
	}
	
	static Shares ComputeConstAND(const Shares& x, std::vector<uint64_t> sigma)
	{
		Shares ret;
		assert(x.s.size() == sigma.size());
		
		ret.s.resize(x.s.size());
		ret.t.resize(x.t.size());
		
		for(int idx = 0; idx < sigma.size(); idx++)
		{
			ret.s[idx] = x.s[idx] & sigma[idx];
			ret.t[idx] = x.t[idx] & sigma[idx];
		}
		
		return ret;
	}
	
	static Shares ComputeConstAND(const Shares& x, uint64_t sigma)
	{
		Shares ret;
// 		assert(x.s.size() == sigma.size());
		
		ret.s.resize(x.s.size());
		ret.t.resize(x.t.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			ret.s[idx] = x.s[idx] & sigma;
			ret.t[idx] = x.t[idx] & sigma;
		}
		
		return ret;
	}
	
	void ComputeConstAND(std::vector<uint64_t> sigma)
	{
		assert(this->s.size() == sigma.size());
		
		for(int idx = 0; idx < sigma.size(); idx++)
		{
			this->s[idx] &= sigma[idx];
			this->t[idx] &= sigma[idx];
		}
	}
	
	static Shares ComputeComplement(const Shares& x)
	{
		Shares ret;
		ret.s.resize(x.s.size());
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			ret.s[idx] = x.s[idx] ^ 0xFFFFFFFFFFFFFFULL;
		}
		
		ret.t = x.t;
		
		return ret;
	}
	
	void ComputeComplement()
	{
		for(int idx = 0; idx < this->s.size(); idx++)
		{
			this->s[idx] ^= 0xFFFFFFFFFFFFFFULL;
		}
	}
	
	void ComputeAND(const Shares& x, Communicator *com)
	{
		this->s.resize(x.s.size());
		this->t.resize(x.s.size());
		
		Shares a;
		a.GenerateRandomShares(x.s.size()*64, com);
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			this->s[idx] = (x.t[idx] & this->t[idx])^(x.s[idx] & this->s[idx])^a.t[idx];
		}
		
		if(com->party == Alice)
		{
			com->SendNext((unsigned char *)(this->s.data()), this->s.size()*sizeof(uint64_t));
			com->AwaitPrevious((unsigned char *)(this->t.data()), this->t.size()*sizeof(uint64_t));
		}
		else
		{
			com->AwaitPrevious((unsigned char *)(this->t.data()), this->t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(this->s.data()), this->s.size()*sizeof(uint64_t));
		}
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			this->t[idx] ^= this->s[idx];
		}
	}
	
	void ComputeAND(const Shares& x, const Shares& y, Communicator *com)
	{
		this->s.resize(x.s.size());
		this->t.resize(x.s.size());
		
		Shares a;
		a.GenerateRandomShares(x.s.size()*64, com);
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			this->s[idx] = (x.t[idx] & y.t[idx])^(x.s[idx] & y.s[idx])^a.t[idx];
		}
		
		if(com->party == Alice)
		{
			com->SendNext((unsigned char *)(this->s.data()), this->s.size()*sizeof(uint64_t));
			com->AwaitPrevious((unsigned char *)(this->t.data()), this->t.size()*sizeof(uint64_t));
		}
		else
		{
			com->AwaitPrevious((unsigned char *)(this->t.data()), this->t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(this->s.data()), this->s.size()*sizeof(uint64_t));
		}
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			this->t[idx] ^= this->s[idx];
		}
	}
	
	static Shares ComputeAND2(const Shares& x, const Shares& y, Communicator *com)
	{
		Shares z;
		z.s.resize(x.s.size());
		z.t.resize(x.s.size());
		
		Shares a;
		a.GenerateRandomShares(x.s.size()*64, com);
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			z.s[idx] = (x.t[idx] & y.t[idx])^(x.s[idx] & y.s[idx])^a.t[idx];
		}
		
		if(com->party == Alice)
		{
			com->SendNext((unsigned char *)(z.s.data()), z.s.size()*sizeof(uint64_t));
			com->AwaitPrevious((unsigned char *)(z.t.data()), z.t.size()*sizeof(uint64_t));
		}
		else
		{
			com->AwaitPrevious((unsigned char *)(z.t.data()), z.t.size()*sizeof(uint64_t));
			com->SendNext((unsigned char *)(z.s.data()), z.s.size()*sizeof(uint64_t));
		}
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			z.t[idx] ^= z.s[idx];
		}
		
		return z;
	}
	
	static Shares ComputeNOT(const Shares& x)
	{
		Shares res;
		
		uint64_t mask = 0xFFFFFFFFFFFFFFFFULL;
		
		res.t = x.t;
		res.s.resize(x.s.size());
		
		for(int idx = 0; idx < x.s.size(); idx++)
		{
			res.s[idx] = x.s[idx] ^ mask;
		}
		
		return res;
	}
	
	
	// Left Shift with 1 padding
	Shares LeftShiftOne()
	{
		Shares x;
		
		int size = this->s.size();
		x.s.resize(size);
		x.t.resize(size);
		
		for(int idx = 0; idx < size - 1; idx++)
		{
			x.s[idx] = (this->s[idx] >> 1);
			x.t[idx] = (this->t[idx] >> 1);
			x.s[idx] += ((this->s[idx + 1] & 1) << 63);
			x.t[idx] += ((this->t[idx + 1] & 1) << 63);
		}
		
		x.s[size - 1] = (this->s[size - 1] >> 1);
		x.t[size - 1] = (this->t[size - 1] >> 1);
		x.s[size - 1] += (1ULL << 63);
		
		return x;
	}
	
	// Right Shift with 1 padding
	Shares RightShiftOne()
	{
		Shares x;
		uint64_t mask = 1;
		mask <<= 63;
		
		int size = this->s.size();
		x.s.resize(size);
		x.t.resize(size);
		
		for(int idx = size-1; idx > 0; idx--)
		{
			x.s[idx] = (this->s[idx] << 1);
			x.t[idx] = (this->t[idx] << 1);
			x.s[idx] += ((this->s[idx - 1] & mask) >> 63);
			x.t[idx] += ((this->t[idx - 1] & mask) >> 63);
		}
		
		x.s[0] = (this->s[0] << 1) + 1ULL;
		x.t[0] = (this->t[0] << 1);
		
		return x;
	}
	
	static std::vector<Shares> RippleCarryAdder(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		assert(x.size() >= y.size());
		
		std::vector<Shares> output(x.size() + 1);
		std::vector<Shares> sum(x.size()), prod(x.size());
		
		Shares carry;
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
			prod[idx] = Shares::ComputeAND2(x[idx], y[idx], com);
		}
		
		std::cout << std::endl;
		
		carry = prod[0];
		output[0] = sum[0];
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			output[idx] = Shares::ComputeXOR(sum[idx], carry);
			carry = Shares::ComputeXOR(prod[idx], Shares::ComputeAND2(sum[idx], carry, com));
		}
		
		output[x.size()] = carry;
		
		return output;
	}
	
// 	static bool IsNonDecreasing(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
// 	{
// 		std::vector<Shares> sum(x.size()), prod1(x.size()), prod2(x.size());
// 		
// 		Shares gt, eq;
// 		gt.s.resize(x[0].s.size());
// 		gt.t.resize(x[0].t.size());
// 		
// 		// Init: t_0 = 0
// 		for(int idx = 0; idx < gt.s.size(); idx++)
// 		{
// 			gt.s[idx] = 0;
// 			gt.t[idx] = 0;
// 		}
// 		
// 		for(int idx = 0; idx < x.size(); idx++)
// 		{
// 			sum[idx] = Shares::ComputeNOT(Shares::ComputeXOR(x[idx], y[idx]));
// 			prod1[idx] = Shares::ComputeAND2(x[idx], Shares::ComputeNOT(y[idx]), com);
// 		}
// 		
// 		eq = sum[0];
// 		
// 		for(int idx = 0; idx < x.size(); idx++)
// 		{
// 			gt = Shares::ComputeXOR(prod1[idx], Shares::ComputeAND2(sum[idx], gt, com));
// 			
// 			if(idx != 0)
// 			{
// 				eq.ComputeAND(sum[idx], com);
// 			}
// 		}
// 		
// 		std::vector<uint64_t> gtVec = gt.OpenAllParty(com);
// 		std::vector<uint64_t> eqVec = eq.OpenAllParty(com);
// 
// 		// Init: t_0 = 0
// 		std::vector<uint64_t> ge(gtVec.size());
// 		for(int idx = 0; idx < ge.size(); idx++)
// 		{
// 			ge[idx] = (gtVec[idx] | ((gtVec[idx] ^ 0xFFFFFFFFFFFFFFFFULL) & eqVec[idx])) ^ 0xFFFFFFFFFFFFFFFFULL;
// 			assert(ge[idx] == 0);
// 		}
// 		
// 		return true;
// 	}
// 	
	std::vector<uint64_t> s;
	std::vector<uint64_t> t;
};

#endif
