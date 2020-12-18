#ifndef TRIPLET__H__
#define TRIPLET__H__

#include <iostream>
#include "Shares.h"
#include "Utility/Timer.h"
#include "Utility/ISecureRNG.h"

using namespace std;
using namespace Utility;

#define BLOCKSIZE 6 // 2^6 = 64

class Triplets
{
public:
	Triplets(){}
	
	Triplets(const Shares& a, const Shares& b, const Shares& c) : a(a), b(b), c(c){}
	
	Triplets(Shares&& a, Shares&& b, Shares&& c) : a(std::move(a)), b(std::move(b)), c(std::move(c)){}
	
	Triplets(const Triplets& other) : a(other.a), b(other.b), c(other.c){}
	
	Triplets(Triplets&& other): a(std::move(other.a)), b(std::move(other.b)), c(std::move(other.c)){}
	
	Triplets& operator=(const Triplets& other)
	{
		if(&other != this)
		{
			this->a = other.a;
			this->b = other.b;
			this->c = other.c;
		}
		
		return *this;
	}
	
	Triplets& operator=(Triplets&& other)
	{
		this->a = std::move(other.a);
		this->b = std::move(other.b);
		this->c = std::move(other.c);
		
		return *this;
	}
	
	~Triplets(){}
	
	void reset(){a.s.resize(0); b.s.resize(0); c.s.resize(0); a.t.resize(0); b.t.resize(0); c.t.resize(0);}
	
	void append(const Shares& x, const Shares y, const Shares& z)
	{
		this->a.s.insert(this->a.s.end(), x.s.begin(), x.s.end());
		this->a.t.insert(this->a.t.end(), x.t.begin(), x.t.end());
		
		this->b.s.insert(this->b.s.end(), y.s.begin(), y.s.end());
		this->b.t.insert(this->b.t.end(), y.t.begin(), y.t.end());
		
		this->c.s.insert(this->c.s.end(), z.s.begin(), z.s.end());
		this->c.t.insert(this->c.t.end(), z.t.begin(), z.t.end());
		
		assert(a.s.size() == a.t.size());
		assert(b.s.size() == b.t.size());
		assert(c.s.size() == c.t.size());
		assert(a.s.size() == b.t.size());
		assert(b.s.size() == c.t.size());
	}
	
	Shares ComputeAND(const Shares& x, const Shares& y, Communicator *com)
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
		
		this->append(x, y, z);
		
		return z;
	}
	
	std::vector<Shares> ComputeAND(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		assert(x.size() == y.size());
		
		int nShares = x.size();
		int length = x[0].s.size();
		
		std::vector<Shares> z(nShares);
		
		// Init size for z
		for(int idx = 0; idx < nShares; idx++)
		{
			z[idx].s.resize(length);
			z[idx].t.resize(length);
		}
		
		Shares a;
		a.GenerateRandomShares(nShares*length*64, com);
		
		for(int sdx = 0; sdx < nShares; sdx++)
		{
			for(int idx = 0; idx < length; idx++)
			{
				z[sdx].s[idx] = (x[sdx].t[idx] & y[sdx].t[idx])^(x[sdx].s[idx] & y[sdx].s[idx])^a.t[sdx*length + idx];
			}
		}
		
		std::vector<uint64_t> sData(0);
		std::vector<uint64_t> tData(nShares*length);
		
		// Copy data to sData
		for(int sdx = 0; sdx < nShares; sdx++) sData.insert(sData.end(), z[sdx].s.begin(), z[sdx].s.end());
		
		if(com->party == Alice)
		{
			com->SendNext((unsigned char *)(sData.data()), nShares*length*sizeof(uint64_t));
			com->AwaitPrevious((unsigned char *)(tData.data()), nShares*length*sizeof(uint64_t));
		}
		else
		{
			com->AwaitPrevious((unsigned char *)(tData.data()), nShares*length*sizeof(uint64_t));
			com->SendNext((unsigned char *)(sData.data()), nShares*length*sizeof(uint64_t));
		}
		
		// Convert back to z[sdx].t
		for(int sdx = 0; sdx < nShares; sdx++)
		{
			for(int idx = 0; idx < length; idx++)
			{
				z[sdx].t[idx] = tData[sdx*length + idx];
			}
		}
		
		for(int sdx = 0; sdx < nShares; sdx++)
		{
			for(int idx = 0; idx < length; idx++)
			{
				z[sdx].t[idx] ^= z[sdx].s[idx];
			}
		}
		
		for(int sdx = 0; sdx < nShares; sdx++) this->append(x[sdx], y[sdx], z[sdx]);
		
		return z;
	}
	
	std::vector<Shares> RippleCarryAdder2(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		assert(x.size() >= y.size());
		
		std::vector<Shares> output(x.size() + 1);
		std::vector<Shares> sum(x.size()), prod(x.size());
		
		Shares carry;
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
		}
		
		prod = this->ComputeAND(x, y, com);
		
		std::cout << std::endl;
		
		carry = prod[0];
		output[0] = sum[0];
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			output[idx] = Shares::ComputeXOR(sum[idx], carry);
			carry = Shares::ComputeXOR(prod[idx], this->ComputeAND(sum[idx], carry, com));
		}
		
		output[x.size()] = carry;
		
		return output;
	}
	
	std::vector<Shares> RippleCarryAdder(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		assert(x.size() >= y.size());
		
		std::vector<Shares> output(x.size() + 1);
		std::vector<Shares> sum(x.size())/*, prod(x.size())*/;
		
		Shares carry;
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
		}
		
// 		prod = this->ComputeAND(x, y, com);
		
		carry = this->ComputeAND(x[0], y[0], com);
		
// 		std::cout << std::endl;
		
// 		carry = prod[0];
		output[0] = sum[0];
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			output[idx] = Shares::ComputeXOR(sum[idx], carry);
			carry = Shares::ComputeXOR(carry, this->ComputeAND(Shares::ComputeXOR(x[idx], carry), Shares::ComputeXOR(y[idx], carry), com));
		}
		
		output[x.size()] = carry;
		
		return output;
	}
	
	std::vector<Shares> RippleCarryAdder(const std::vector<Shares>& x, const std::vector<Shares>& y, const std::vector<Shares>& z, Communicator *com)
	{
		assert(x.size() == y.size());
		assert(x.size() == z.size());
		
		int length = x[0].s.size();
		
		std::vector<Shares> c(x.size() + 1), s(x.size()+1);
		
		// c[0]: shares of zeros
		c[0].s.resize(length);
		c[0].t.resize(length);
		s[x.size()].s.resize(length);
		s[x.size()].t.resize(length);
		
		std::vector<Shares> sum(x.size());
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
			s[idx] = Shares::ComputeXOR(sum[idx], z[idx]);
		}
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			c[idx + 1] = Shares::ComputeXOR(z[idx], this->ComputeAND(Shares::ComputeXOR(z[idx], x[idx]), Shares::ComputeXOR(z[idx], y[idx]), com));
		}
		
		std::cout << "Phase 2 adder" << std::endl;
		std::vector<Shares> output = this->RippleCarryAdder(c, s, com);
		
		return output;
	}
	
	std::vector<Shares> RippleCarryAdder2(const std::vector<Shares>& x, const std::vector<Shares>& y, const std::vector<Shares>& z, Communicator *com)
	{
		assert(x.size() == y.size());
		assert(x.size() == z.size());
		
		int length = x[0].s.size();
		
		std::vector<Shares> c(x.size() + 1), s(x.size()+1);
		
		// c[0]: shares of zeros
		c[0].s.resize(length);
		c[0].t.resize(length);
		s[x.size()].s.resize(length);
		s[x.size()].t.resize(length);
		
		std::vector<Shares> sum(x.size());
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
			s[idx] = Shares::ComputeXOR(sum[idx], z[idx]);
		}
		
		sum.insert(sum.end(), x.begin(), x.end());
		std::vector<Shares> zz = z;
		zz.insert(zz.end(), y.begin(), y.end());
		
		std::vector<Shares> temp = this->ComputeAND(sum, zz, com);
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			c[idx + 1] = Shares::ComputeXOR(temp[idx], temp[idx + x.size()]);
		}
		
		std::cout << "Phase 2 adder" << std::endl;
		std::vector<Shares> output = this->RippleCarryAdder(c, s, com);
		
		return output;
	}
	
	void modp(std::vector<Shares>& shares, const ZZ& p, Communicator *com)
	{
		int length = shares[0].s.size();
		
		ZZ p2 = p;
		
		Shares gt = IsGreaterThanConst(shares, p, com);
		Shares zero; zero.s.resize(length, 0); zero.t.resize(length, 0);
		
		std::vector<Shares> sp(shares.size());
		
		for(int idx = 0; idx < shares.size(); idx++)
		{
			if(IsZero(p2 & 1))
			{
				sp[idx] = zero;
			}
			else
			{
				sp[idx] = gt;
			}
			
			p2 >>= 1;
		}
		
		shares = RippleBorrowSubtractor(shares, sp, com);
		shares.resize(shares.size() - 1);
		
		gt = IsGreaterThanConst(shares, p, com);
		
		p2 = p;
		for(int idx = 0; idx < shares.size(); idx++)
		{
			if(IsZero(p2 & 1))
			{
				sp[idx] = zero;
			}
			else
			{
				sp[idx] = gt;
			}
			
			p2 >>= 1;
		}
		
		shares = RippleBorrowSubtractor(shares, sp, com);
		
		shares.resize(shares.size() - 1);
	}
	
	void modp_fast(std::vector<Shares>& shares, const ZZ& p, Communicator *com)
	{
		int length = shares[0].s.size();
		
		ZZ p2 = p;
		
		Shares gt = shares[shares.size() - 2];
		Shares zero; zero.s.resize(length, 0); zero.t.resize(length, 0);
		
		std::vector<Shares> sp(shares.size()), s2p(shares.size());
		
		for(int idx = 0; idx < shares.size(); idx++)
		{
			if(IsOdd(p >> idx))
			{
				sp[idx] = gt;
			}
			else
			{
				sp[idx] = zero;
			}
		}
		
		s2p[0] = zero;
		for(int idx = 0; idx < shares.size()-1; idx++)
		{
			if(IsOdd(p >> idx))
			{
				s2p[idx+1] = gt;
			}
			else
			{
				s2p[idx+1] = zero;
			}
		}
		
		for(int idx = 0; idx < shares.size(); idx++) sp[idx] = Shares::ComputeXOR(sp[idx], s2p[idx]);
		
		shares = RippleBorrowSubtractor(shares, sp, com);
		shares.resize(shares.size() - 2);
	}
	
	std::vector<Shares> RippleBorrowSubtractor(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		assert(x.size() >= y.size());
		
		std::vector<Shares> output(x.size());
		std::vector<Shares> sum(x.size()), prod(x.size());
		
		Shares borrow; 
		
		output[0] = Shares::ComputeXOR(x[0], y[0]);
		borrow = this->ComputeAND(y[0], Shares::ComputeNOT(x[0]), com);
		
		for(int idx = 1; idx < x.size()-2; idx++)
		{
			output[idx] = Shares::ComputeXOR(Shares::ComputeXOR(x[idx], borrow), y[idx]);
			borrow = Shares::ComputeXOR(borrow, this->ComputeAND(Shares::ComputeNOT(Shares::ComputeXOR(x[idx], borrow)), Shares::ComputeXOR(y[idx], borrow), com));
		}
		
		return output;
	}
	
	std::vector<Shares> RippleBorrowSubtractor2(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		assert(x.size() >= y.size());
		
		std::vector<Shares> output(x.size());
		std::vector<Shares> sum(x.size()), prod(x.size());
		
		Shares borrow;
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
		}
		
		std::vector<Shares> ix(x.size()); 
		for(int idx = 0; idx < x.size(); idx++) ix[idx] = Shares::ComputeNOT(x[idx]);
		
		prod = this->ComputeAND(ix, y, com);
		
		borrow = prod[0];
		output[0] = sum[0];
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			output[idx] = Shares::ComputeXOR(sum[idx], borrow);
			borrow = Shares::ComputeXOR(prod[idx], this->ComputeAND(Shares::ComputeNOT(sum[idx]), borrow, com));
		}
		
		return output;
	}
	
	Shares IsGreaterThanConst(const std::vector<Shares>& x, const ZZ& p, Communicator *com)
	{
		int length = x.size();
		ZZ p2 = p;
		
		std::vector<uint64_t> pVec(length);
		
		for(int idx = 0; idx < length; idx++)
		{
			if(IsZero(p2 & 1))
			{
				pVec[idx] = 0;
			}
			else
			{
				pVec[idx] = 0xFFFFFFFFFFFFFFFFULL;
			}
			
			p2 >>= 1;
		}
		
		Shares gt;
		std::vector<Shares> sum(length), prod(length);
		gt.s.resize(x[0].s.size());
		gt.t.resize(x[0].t.size());
		
		int count = 0;
		for(int idx = 0; idx < length; idx++)
		{
			sum[idx] = Shares::ComputeNOT(Shares::ComputeConstXOR(x[idx], pVec[count]));
			prod[idx] = Shares::ComputeConstAND(x[idx], (pVec[count] ^ 0xFFFFFFFFFFFFFFFFULL));
			count++;
		}
		
		for(int idx = 0; idx < length; idx++)
		{
			gt = Shares::ComputeXOR(prod[idx], this->ComputeAND(sum[idx], gt, com));
		}
		
		return gt;
	}
	
	Shares IsGreaterThanConst2(const std::vector<Shares>& x, const ZZ& p, Communicator *com)
	{
		int length = x.size();
		ZZ p2 = p;
		
		std::vector<uint64_t> pVec(length);
		
		for(int idx = 0; idx < length; idx++)
		{
			if(IsZero(p2 & 1))
			{
				pVec[idx] = 0;
			}
			else
			{
				pVec[idx] = 0xFFFFFFFFFFFFFFFFULL;
			}
			
			p2 >>= 1;
		}
		
		Shares gt;
		std::vector<Shares> sum(length), prod(length);
		gt.s.resize(x[0].s.size());
		gt.t.resize(x[0].t.size());
		
		int count = 0;
		for(int idx = 0; idx < length; idx++)
		{
			sum[idx] = Shares::ComputeNOT(Shares::ComputeConstXOR(x[idx], pVec[count]));
			prod[idx] = Shares::ComputeConstAND(x[idx], (pVec[count] ^ 0xFFFFFFFFFFFFFFFFULL));
			count++;
		}
		
		for(int idx = 0; idx < length; idx++)
		{
			gt = Shares::ComputeXOR(prod[idx], this->ComputeAND(sum[idx], gt, com));
		}
		
		return gt;
	}
	
// 	bool IsNonDecreasing(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
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
// 			prod1[idx] = this->ComputeAND(x[idx], Shares::ComputeNOT(y[idx]), com);
// 		}
// 		
// 		eq = sum[0];
// 		
// 		for(int idx = 0; idx < x.size(); idx++)
// 		{
// 			gt = Shares::ComputeXOR(prod1[idx], this->ComputeAND(sum[idx], gt, com));
// 			
// 			if(idx != 0)
// 			{
// 				eq = this->ComputeAND(eq, sum[idx], com);
// 			}
// 		}
// 		
// 		std::vector<uint64_t> gtVec = gt.OpenAllParty(com);
// 		std::vector<uint64_t> eqVec = eq.OpenAllParty(com);
// 
// 		// Init: t_0 = 0
// 		bool res = true;
// 		std::vector<uint64_t> ge(gtVec.size());
// 		for(int idx = 0; idx < ge.size(); idx++)
// 		{
// 			ge[idx] = (gtVec[idx] | ((gtVec[idx] ^ 0xFFFFFFFFFFFFFFFFULL) & eqVec[idx])) ^ 0xFFFFFFFFFFFFFFFFULL;
// 			if(ge[idx] != 0) res = false;
// 		}
// 		
// 		return res;
// 	}
// 	
	// c_(i+1) = c_i \XOR [(c_i \XOR x_i) AND NOT(c_i \XOR y_i)]
	bool IsIncreasing(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		std::vector<Shares> sum(x.size());
		
		Shares gt;
		gt.s.resize(x[0].s.size());
		gt.t.resize(x[0].t.size());
		
		// Init: t_0 = 0
		for(int idx = 0; idx < gt.s.size(); idx++)
		{
			gt.s[idx] = 0;
			gt.t[idx] = 0;
		}
		
		gt = this->ComputeAND(x[0], Shares::ComputeNOT(y[0]), com);
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			gt = Shares::ComputeXOR(gt, this->ComputeAND(Shares::ComputeXOR(gt, x[idx]), Shares::ComputeNOT(Shares::ComputeXOR(gt, y[idx])), com));
		}
		
		std::vector<uint64_t> gtVec = gt.OpenAllParty(com);

		// Init: t_0 = 0
		bool res = true;
		std::vector<uint64_t> ge(gtVec.size());
		for(int idx = 0; idx < ge.size(); idx++)
		{
			gtVec[idx] ^= 0xFFFFFFFFFFFFFFFFULL;
			if(gtVec[idx] != 0) res = false;
		}
		
		return res;
	}
	
	bool IsIncreasing2(const std::vector<Shares>& x, const std::vector<Shares>& y, Communicator *com)
	{
		std::vector<Shares> sum(x.size()), prod1(x.size());
		
		Shares gt;
		gt.s.resize(x[0].s.size());
		gt.t.resize(x[0].t.size());
		
		// Init: t_0 = 0
		for(int idx = 0; idx < gt.s.size(); idx++)
		{
			gt.s[idx] = 0;
			gt.t[idx] = 0;
		}
		
		std::vector<Shares> iy(y.size());
		for(int idx = 0; idx < y.size(); idx++) iy[idx] = Shares::ComputeNOT(y[idx]);
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeNOT(Shares::ComputeXOR(x[idx], y[idx]));
// 			prod1[idx] = this->ComputeAND(x[idx], Shares::ComputeNOT(y[idx]), com);
		}
		
		prod1 = this->ComputeAND(x, iy, com);
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			gt = Shares::ComputeXOR(prod1[idx], this->ComputeAND(sum[idx], gt, com));
		}
		
		std::vector<uint64_t> gtVec = gt.OpenAllParty(com);

		// Init: t_0 = 0
		bool res = true;
		std::vector<uint64_t> ge(gtVec.size());
		for(int idx = 0; idx < ge.size(); idx++)
		{
			gtVec[idx] ^= 0xFFFFFFFFFFFFFFFFULL;
			if(gtVec[idx] != 0) res = false;
		}
		
		return res;
	}
	
	void CompareAndSwap(std::vector<Shares>& x, std::vector<Shares>& y, Communicator *com)
	{
		// Compare
		Shares gt;
		gt.s.resize(x[0].s.size());
		gt.t.resize(x[0].t.size());
		
		// Init: t_0 = 0
		for(int idx = 0; idx < gt.s.size(); idx++)
		{
			gt.s[idx] = 0;
			gt.t[idx] = 0;
		}
		
		gt = this->ComputeAND(x[0], Shares::ComputeNOT(y[0]), com);
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			gt = Shares::ComputeXOR(gt, this->ComputeAND(Shares::ComputeXOR(gt, x[idx]), Shares::ComputeNOT(Shares::ComputeXOR(gt, y[idx])), com));
		}
		
		// Swap
		// x <-- (x+y)s + x	y <-- (x+y)s + y
		std::vector<Shares> sum(x.size());
		std::vector<Shares> sign(x.size());
		std::vector<Shares> prod;
// 		CompareAndSwap
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
			sign[idx] = gt;
		}
		
		prod = this->ComputeAND(sum, sign, com);
		for(int idx = 0; idx < x.size(); idx++)
		{
			x[idx] = Shares::ComputeXOR(x[idx], prod[idx]);
			y[idx] = Shares::ComputeXOR(y[idx], prod[idx]);
		}
	}
	
	void PairwiseComparison(std::vector<Shares>& x, std::vector<Shares>& y, Communicator *com)
	{
		std::vector<Shares> sum(x.size());
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
		}
		
		Shares eq = sum[0];
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			eq = this->ComputeAND(eq, sum[idx], com);
		}
		
		std::vector<Shares> eqs(x.size());
		for(int idx = 0; idx < x.size(); idx++)
		{
			eqs[idx] = eq;
		}
		
		x = this->ComputeAND(x, eqs, com);
	}
	
	void PairwiseComparison2(Shares& eq, std::vector<Shares>& x, std::vector<Shares>& y, Communicator *com)
	{
		std::vector<Shares> sum(x.size());
		
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
		}
		
		eq = sum[0];
		
		for(int idx = 1; idx < x.size(); idx++)
		{
			eq = this->ComputeAND(eq, sum[idx], com);
		}
	}
	
	void Swapping(std::vector<Shares>& x, std::vector<Shares>& y, Communicator *com)
	{
		// Compare
		Shares gt;
// 		gt.s.resize(x[0].s.size());
// 		gt.t.resize(x[0].t.size());
// 		
		gt.GenerateRandomShares(64ULL*x[0].s.size(), com);
// 		// Init: t_0 = 0
// 		for(int idx = 0; idx < gt.s.size(); idx++)
// 		{
// 			gt.s[idx] = 0;
// 			gt.t[idx] = 0;
// 		}
		
// 		gt = this->ComputeAND(x[0], Shares::ComputeNOT(y[0]), com);
// 		
// 		for(int idx = 1; idx < x.size(); idx++)
// 		{
// 			gt = Shares::ComputeXOR(gt, this->ComputeAND(Shares::ComputeXOR(gt, x[idx]), Shares::ComputeNOT(Shares::ComputeXOR(gt, y[idx])), com));
// 		}
		
		// Swap
		// x <-- (x+y)s + x	y <-- (x+y)s + y
		std::vector<Shares> sum(x.size());
		std::vector<Shares> sign(x.size());
		std::vector<Shares> prod;
// 		CompareAndSwap
		for(int idx = 0; idx < x.size(); idx++)
		{
			sum[idx] = Shares::ComputeXOR(x[idx], y[idx]);
			sign[idx] = gt;
		}
		
		prod = this->ComputeAND(sum, sign, com);
		for(int idx = 0; idx < x.size(); idx++)
		{
			x[idx] = Shares::ComputeXOR(x[idx], prod[idx]);
			y[idx] = Shares::ComputeXOR(y[idx], prod[idx]);
		}
	}
	
	void TripleVerificationWithOpening(Communicator *com)
	{
		int length = a.s.size();
		
		Shares z = this->a;
		z.s.insert(z.s.end(), this->b.s.begin(), this->b.s.end());
		z.s.insert(z.s.end(), this->c.s.begin(), this->c.s.end());
		z.t.insert(z.t.end(), this->b.t.begin(), this->b.t.end());
		z.t.insert(z.t.end(), this->c.t.begin(), this->c.t.end());
		
		std::vector<uint64_t> rz = z.OpenAllParty(com);
		
		for(int idx = 0; idx < length; idx++)
		{
			assert(rz[idx + 2*length] == (rz[idx + length] & rz[idx]));
		}
	}
	
	void TripleVerificationWithoutOpening(const Triplets& xyz, Communicator *com)
	{
		Timer t;
		Shares rho = Shares::ComputeXOR(xyz.a, this->a);
		Shares sigma = Shares::ComputeXOR(xyz.b, this->b);
		
// 		t.Tick("Compute XOR");
		
		int length = rho.s.size();
		
		rho.s.insert(rho.s.end(), (sigma.s.begin()), (sigma.s.end()));
		rho.t.insert(rho.t.end(), (sigma.t.begin()), (sigma.t.end()));
		
// 		std::cout << "Reveal rho" << std::endl;
		std::vector<uint64_t> reveal_rho = rho.OpenAllParty(com);
		
// 		t.Tick("Reveal rho & sigma");
		
		Shares res;
		
		std::vector<uint64_t> reveal_sigma, view;
		view = reveal_rho;
		reveal_sigma.insert(reveal_sigma.end(), std::make_move_iterator(reveal_rho.begin() + length), std::make_move_iterator(reveal_rho.end()));
		reveal_rho.resize(length);
		
// 		t.Tick("CompareView");
		
		res = res.ComputeXOR(res.ComputeXOR(xyz.c, this->c), res.ComputeXOR(res.ComputeConstAND(this->a, reveal_sigma), res.ComputeConstAND(this->b, reveal_rho)));
		
		
		for(int idx = 0; idx < res.s.size(); idx++)
		{
			reveal_rho[idx] &= reveal_sigma[idx];
		}
		
		res.ComputeConstXOR(reveal_rho);
		
// 		t.Tick("Compute beaver");
		
		res.CompareView2(view, com);
		
// 		t.Tick("Triplet Verification");
	}
	
	void VerifyTriples(std::vector<unsigned char>& rngseed, Communicator *com)
	{
		Timer t;
		
		// Generate common rng
		AESRNG *rng = new AESRNG(rngseed.data());
		
		uint64_t size = this->a.s.size();
		uint64_t numTriples = (2*size + 1)*64;
		
		std::cout << "Number of tripples: " << this->a.s.size()*64 << std::endl;
		
		Triplets triples;
		triples.a.GenerateRandomShares(numTriples, com);
		triples.b.GenerateRandomShares(numTriples, com);
		
// 		t.Tick("Generate random triples");
		
		triples.c.ComputeAND(triples.a, triples.b, com);
		
// 		t.Tick("AND random triples");
		
		std::vector<uint64_t> seed = rng->GetUInt64Array(2);
		
		int randomLoc = seed[1] % (2*size + 1);
		
		Triplets testTriple;
		
		testTriple.a.s.resize(1); testTriple.a.t.resize(1);
		testTriple.b.s.resize(1); testTriple.b.t.resize(1);
		testTriple.c.s.resize(1); testTriple.c.t.resize(1);
		
		testTriple.a.s[0] = triples.a.s[randomLoc];
		testTriple.a.t[0] = triples.a.t[randomLoc];
		
		testTriple.b.s[0] = triples.b.s[randomLoc];
		testTriple.b.t[0] = triples.b.t[randomLoc];
		
		testTriple.c.s[0] = triples.c.s[randomLoc];
		testTriple.c.t[0] = triples.c.t[randomLoc];
		
		testTriple.TripleVerificationWithOpening(com);
		
		triples.a.s.erase(triples.a.s.begin() + randomLoc);
		triples.b.s.erase(triples.b.s.begin() + randomLoc);
		triples.c.s.erase(triples.c.s.begin() + randomLoc);
		triples.a.t.erase(triples.a.t.begin() + randomLoc);
		triples.b.t.erase(triples.b.t.begin() + randomLoc);
		triples.c.t.erase(triples.c.t.begin() + randomLoc);
		
// 		t.Tick("Check random triple");
		
// 		std::vector<int> index(1<<20);
// 		std::iota (std::begin(index), std::end(index), 0);
// 		std::shuffle(index.begin(), index.end(), std::mt19937(seed[0]));
		
// 		t.Tick("Init shuffle index");
		
		if(triples.a.s.size() <= (1 << 20))
		{
			std::shuffle(triples.a.s.begin(), triples.a.s.end(), std::mt19937(seed[0]));
			std::shuffle(triples.a.t.begin(), triples.a.t.end(), std::mt19937(seed[0]));
			
			std::shuffle(triples.b.s.begin(), triples.b.s.end(), std::mt19937(seed[0]));
			std::shuffle(triples.b.t.begin(), triples.b.t.end(), std::mt19937(seed[0]));
			
			std::shuffle(triples.c.s.begin(), triples.c.s.end(), std::mt19937(seed[0]));
			std::shuffle(triples.c.t.begin(), triples.c.t.end(), std::mt19937(seed[0]));
		}
		else
		{
			int groupSize = triples.a.s.size()/(1 << 20);
			std::vector<std::vector<uint64_t> > temp(1 << 20);
			
			for(int idx = 0; idx < temp.size(); idx++) temp[idx].resize(groupSize);
			
// 			t.Tick("Init temp");
			
			// Process a.s
			int count;
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					temp[idx][kdx] = triples.a.s[count]; count++;
				}
			}
			
			std::shuffle(temp.begin(), temp.end(), std::mt19937(seed[0]));
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					triples.a.s[count] = temp[idx][kdx]; count++;
				}
			}
			
// 			t.Tick("Round 1: a.s");
			
			// a.t
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					temp[idx][kdx] = triples.a.t[count]; count++;
				}
			}
			
			std::shuffle(temp.begin(), temp.end(), std::mt19937(seed[0]));
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					triples.a.t[count] = temp[idx][kdx]; count++;
				}
			}
			
// 			t.Tick("Round 1: a.t");
			
			// b.s & b.t
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					temp[idx][kdx] = triples.b.s[count]; count++;
				}
			}
			
			std::shuffle(temp.begin(), temp.end(), std::mt19937(seed[0]));
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					triples.b.s[count] = temp[idx][kdx]; count++;
				}
			}
			
// 			t.Tick("Round 2: b.s");
			
			// b.t
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					temp[idx][kdx] = triples.b.t[count]; count++;
				}
			}
			
			std::shuffle(temp.begin(), temp.end(), std::mt19937(seed[0]));
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					triples.b.t[count] = temp[idx][kdx]; count++;
				}
			}
			
// 			t.Tick("Round 2: b.t");
			
			// c.s & c.t
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					temp[idx][kdx] = triples.c.s[count]; count++;
				}
			}
			
			std::shuffle(temp.begin(), temp.end(), std::mt19937(seed[0]));
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					triples.c.s[count] = temp[idx][kdx]; count++;
				}
			}
			
// 			t.Tick("Round 2: c.s");
			
			// c.t
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					temp[idx][kdx] = triples.c.t[count]; count++;
				}
			}
			
			std::shuffle(temp.begin(), temp.end(), std::mt19937(seed[0]));
			
			count = 0;
			for(int idx = 0; idx < temp.size(); idx++)
			{
				for(int kdx = 0; kdx < groupSize; kdx++)
				{
					triples.c.t[count] = temp[idx][kdx]; count++;
				}
			}
// 			t.Tick("Round 2: c.t");
		}
		
// 		t.Tick("Shuffle");
		
		
		
		Triplets verify1, verify2;
		
// 		t.Tick("Open triples");
		
		verify2.a.s.insert(verify2.a.s.end(), std::make_move_iterator(triples.a.s.begin() + size), std::make_move_iterator(triples.a.s.begin() + 2*size));
		verify2.b.s.insert(verify2.b.s.end(), std::make_move_iterator(triples.b.s.begin() + size), std::make_move_iterator(triples.b.s.begin() + 2*size));
		verify2.c.s.insert(verify2.c.s.end(), std::make_move_iterator(triples.c.s.begin() + size), std::make_move_iterator(triples.c.s.begin() + 2*size));
		verify2.a.t.insert(verify2.a.t.end(), std::make_move_iterator(triples.a.t.begin() + size), std::make_move_iterator(triples.a.t.begin() + 2*size));
		verify2.b.t.insert(verify2.b.t.end(), std::make_move_iterator(triples.b.t.begin() + size), std::make_move_iterator(triples.b.t.begin() + 2*size));
		verify2.c.t.insert(verify2.c.t.end(), std::make_move_iterator(triples.c.t.begin() + size), std::make_move_iterator(triples.c.t.begin() + 2*size));
		
		triples.a.s.resize(size); triples.a.t.resize(size);
		triples.b.s.resize(size); triples.b.t.resize(size);
		triples.c.s.resize(size); triples.c.t.resize(size);
		
		verify1 = std::move(triples);
		
// 		t.Tick("Copy time");
		
		this->TripleVerificationWithoutOpening(verify1, com);
		this->TripleVerificationWithoutOpening(verify2, com);
		
		t.Tick("Verifying Triplets");
		
		delete rng;
	}
	
	void swap(int i, int j)
	{
		int r1 = i >> BLOCKSIZE;
		int c1 = i & ((1 << BLOCKSIZE) - 1);
		
		int r2 = j >> BLOCKSIZE;
		int c2 = j & ((1 << BLOCKSIZE) - 1);
		
		// swap a
		if(((a.s[r1] >> c1) & 1) != ((a.s[r2] >> c2) & 1))
		{
			a.s[r1] ^= (1ULL << c1);
			a.s[r2] ^= (1ULL << c2);
		}
		
		if(((a.t[r1] >> c1) & 1) != ((a.t[r2] >> c2) & 1))
		{
			a.t[r1] ^= (1ULL << c1);
			a.t[r2] ^= (1ULL << c2);
		}
		
		if(((b.s[r1] >> c1) & 1) != ((b.s[r2] >> c2) & 1))
		{
			b.s[r1] ^= (1ULL << c1);
			b.s[r2] ^= (1ULL << c2);
		}
		
		if(((b.t[r1] >> c1) & 1) != ((b.t[r2] >> c2) & 1))
		{
			b.t[r1] ^= (1ULL << c1);
			b.t[r2] ^= (1ULL << c2);
		}
		
		if(((c.s[r1] >> c1) & 1) != ((c.s[r2] >> c2) & 1))
		{
			c.s[r1] ^= (1ULL << c1);
			c.s[r2] ^= (1ULL << c2);
		}
		
		if(((c.t[r1] >> c1) & 1) != ((c.t[r2] >> c2) & 1))
		{
			c.t[r1] ^= (1ULL << c1);
			c.t[r2] ^= (1ULL << c2);
		}
	}
	
	Shares a;
	Shares b;
	Shares c;
};

#endif
