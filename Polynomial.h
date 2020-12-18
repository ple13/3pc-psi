#ifndef POLYNOMIAL_H__
#define POLYNOMIAL_H__

#include <cstddef>
#include <cassert>
#include <vector>
#include <chrono> 
#include <algorithm>
#include <cstring>
#include <math.h>
#include <sys/time.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_limbs.h>
#include <NTL/BasicThreadPool.h>
#include <future>
#include <thread>
#include "Utility/Timer.h"

#define VERBOSE 0
#define PARALLEL 1

using namespace std;
using namespace NTL;

void polynomialMod(ZZ_pX& res, ZZ&& n, const ZZ_pXModulus& g);
void polynomialMod(ZZ_pX& res, ZZ& n, const ZZ_pXModulus& g);
void generateInput(Vec<ZZ_p>& vals, ZZ_pX& f);
void printPolynomial(const Vec<ZZ_p>& f, const std::string& s, int k);
void printPolynomial(const ZZ_pX& f, const std::string& s, int k);
void uptree_multiplication(std::vector<std::vector<ZZ_pX> >& tree, const Vec<ZZ_p>& vals, bool interpolation);
void downtree_division(std::vector<std::vector<ZZ_pX> >& downtree, std::vector<std::vector<ZZ_pX> >& uptree, const ZZ_pX& input);
void reconstructPolynomial(std::vector<ZZ_pX>& coeffs, const Vec<ZZ_p>& evals,  const std::vector<std::vector<ZZ_pX> >& interpolationUpTree, const std::vector<std::vector<ZZ_pX> >& interpolationDownTree);
void evaluatePolynomial(Vec<ZZ_p>& evals, const ZZ_pX& f, const Vec<ZZ_p>& vals);
void interpolatePolynomial(ZZ_pX& f, const Vec<ZZ_p>& vals, const Vec<ZZ_p>& evals);

class PolynomialTree
{
public:
	PolynomialTree()
	{
		left = nullptr;
		right = nullptr;
	}
	
	~PolynomialTree()
	{
		if(left != nullptr)
		{
			delete left;
		}
		
		if(right != nullptr)
		{
			delete right;
		}
	}
	
	void printTree()
	{
		printPolynomial(val, "val", 10);
	}
	
	ZZ_pX val;
	ZZ_pX interp;
	
	PolynomialTree *left;
	PolynomialTree *right;
};

void buildTree(PolynomialTree *tree, const ZZ_p *roots, int size);
void downtreeDivisionReconstruct(ZZ_p *evals, PolynomialTree *tree, const ZZ_pX& f, int pos, int size);
void downtreeDivision(ZZ_p *evals, PolynomialTree *tree, const ZZ_pX& f, int pos, int size);
void reconstruct(ZZ_p *evals, ZZ_p *derivativeEvals, PolynomialTree *tree, int pos, int size);
void evaluatePolynomial(ZZ_p *evals, const ZZ_pX& f, const ZZ_p *roots, int size);
void interpolatePolynomial(ZZ_pX& f, ZZ_p *roots, ZZ_p *evals, int size);
int test();

#endif
