///  g++ -march=native -mfpmath=sse -O3 mail.cpp -lntl -lgmp -std=c++11 -pthread
#include "Polynomial.h"

// res = x^n mod g
void polynomialMod(ZZ_pX& res, ZZ&& n, const ZZ_pXModulus& g)
{
	if((n + 2 - deg(g)) < 0)
	{
		int m = conv<int>(n);
		res.SetLength(m+1);
		res[m] = 1;
	}
	else
	{
		
		if (n % 2 == 0)
		{
			polynomialMod(res, (n >> 1), g);
			MulMod(res, res, res, g);
		}
		else
		{
			polynomialMod(res, (n >> 1), g);
			ZZ_pX temp;
			temp.SetLength(2);
			temp[0] = 0; temp[1] = 1;
			mul(temp, temp, res);
			MulMod(res, res, temp, g);
		}
	}
}

void polynomialMod(ZZ_pX& res, ZZ& n, const ZZ_pXModulus& g)
{
	if((n + 2 - deg(g)) < 0)
	{
		int m = conv<int>(n);
		res.SetLength(m+1);
		res[m] = 1;
	}
	else
	{
		
		if (n % 2 == 0)
		{
			polynomialMod(res, (n >> 1), g);
			MulMod(res, res, res, g);
		}
		else
		{
			polynomialMod(res, (n >> 1), g);
			ZZ_pX temp;
			temp.SetLength(2);
			temp[0] = 0; temp[1] = 1;
			mul(temp, temp, res);
			MulMod(res, res, temp, g);
		}
	}
}

void generateInput(Vec<ZZ_p>& vals, ZZ_pX& f)
{
	std::cout << "N points: ";
	for(int idx = 0; idx < vals.length(); idx++)
	{
		random(vals[idx]);
// 		vals[idx] = (idx + 1);
		f[idx] = vals[idx];
		if(vals.length() < 50) std::cout << vals[idx] << " ";
	}
	std::cout << std::endl;
}

void printPolynomial(const Vec<ZZ_p>& f, const std::string& s, int k)
{
	if(f.length() < k) k = f.length();
	
	std::cout << s << ": ";
	for (size_t i = 0; i < k; i++)
	{
		std::cout << f[i] << "x^" << i << " + ";
	}
	std::cout << std::endl;
}

void printPolynomial(const ZZ_pX& f, const std::string& s, int k)
{
	if(deg(f) < k) k = deg(f) + 1;
	std::cout << s << ": ";
	for (size_t i = 0; i < k; i++)
	{
		std::cout << f[i] << "x^" << i << " + ";
	}
	std::cout << std::endl;
}

void uptree_multiplication(std::vector<std::vector<ZZ_pX> >& tree, const Vec<ZZ_p>& vals, bool interpolation)
{
	int k = ceil(log2(vals.length()));
	
	// Get the full polynomial in case we do interpolation
	if(interpolation) k += 1;
	
// 	std::cout << "k: " << k << std::endl;
	tree.resize(k);
	
	// Init the bottom layer of the tree
	tree[0].resize(vals.length());
	
	for(int idx = 0; idx < vals.length(); idx++)
	{
		tree[0][idx].SetLength(2);
		tree[0][idx][0] = -vals[idx];
		tree[0][idx][1] = 1;
	}
	
	int s = 2;
	
	for(int idx = 1; idx < k; idx++)
	{
		tree[idx].resize((tree[idx-1].size() + 1) / 2);
		
		std::cout << "idx: " << idx << "\t " << tree[idx].size() << "\t";
		Timer t;
		
		if(PARALLEL)
		{
			ZZ_pContext context;
			context.save();
			
			NTL_EXEC_RANGE(tree[idx].size(), first, last) 
				context.restore();
				for(long kdx = first; kdx < last; kdx++)
				{
					mul(tree[idx][kdx], tree[idx-1][2*kdx], tree[idx-1][2*kdx + 1]);
				}
			
			NTL_EXEC_RANGE_END
		}
		else
		{
			for(int kdx = 0; kdx < tree[idx].size(); kdx++)
			{
				mul(tree[idx][kdx], tree[idx-1][2*kdx], tree[idx-1][2*kdx + 1]);
				
				if(VERBOSE) 
				{
					printPolynomial(tree[idx-1][2*kdx], "tree[idx-1][2*kdx]", 20);
					printPolynomial(tree[idx-1][2*kdx + 1], "tree[idx-1][2*kdx + 1]", 20);
					printPolynomial(tree[idx][kdx], "tree[idx][kdx]", 20);
				}
			}
		}
		
// 		t.Tick("");
	}
}

void downtree_division(std::vector<std::vector<ZZ_pX> >& downtree, std::vector<std::vector<ZZ_pX> >& uptree, const ZZ_pX& input)
{
	downtree.resize(uptree.size() + 1);
	
	downtree[downtree.size() - 1].resize(1);
	downtree[downtree.size() - 1][0].SetLength(deg(input) + 1);
	
	for(int idx = 0; idx < deg(input) + 1; idx++)
	{
		downtree[downtree.size() - 1][0][idx] = input[idx];
	}
	
	for(int idx = downtree.size() - 2; idx >= 0; idx--)
	{
		std::cout << "idx: " << idx << " " << uptree[idx].size() << "\t";
		Timer t;
		
		downtree[idx].resize(uptree[idx].size());
		
		if(PARALLEL)
		{
			ZZ_pContext context;
			context.save();
			
			NTL_EXEC_RANGE(downtree[idx].size(), first, last) 
				context.restore();
				for(long kdx = first; kdx < last; kdx++)
				{
					rem(downtree[idx][kdx], downtree[idx+1][kdx/2], uptree[idx][kdx]);
				}
			NTL_EXEC_RANGE_END
		}
		else
		{
			for(int kdx = 0; kdx < downtree[idx].size(); kdx++)
			{
				rem(downtree[idx][kdx], downtree[idx+1][kdx/2], uptree[idx][kdx]);
			}
		}
		
// 		t.Tick("");
	}
	
	if(VERBOSE) 
	{
		for(int idx = 0; idx < downtree[0].size(); idx++)
		{
			std::cout << "eval[" << idx << "]: " << downtree[0][idx][0] << std::endl;
		}
	}
}

void reconstructPolynomial(std::vector<ZZ_pX>& coeffs, const Vec<ZZ_p>& evals,  const std::vector<std::vector<ZZ_pX> >& interpolationUpTree, const std::vector<std::vector<ZZ_pX> >& interpolationDownTree)
{
	int n = evals.length();
	coeffs.resize(n);
	
	std::vector<ZZ_pX> r0, r1;
	
	for(int idx = 0; idx < n; idx++)
	{
		coeffs[idx].SetLength(1);
		coeffs[idx][0] = (evals[idx]*inv(interpolationDownTree[0][idx][0]));
	}
	
	if(PARALLEL)
	{
		for(int idx = 0; idx < interpolationUpTree.size(); idx++)
		{
			Timer t;
			r0.resize(interpolationUpTree[idx].size()/2);
			r1.resize(interpolationUpTree[idx].size()/2);
			
			ZZ_pContext context;
			context.save();
			
			NTL_EXEC_RANGE(interpolationUpTree[idx].size()/2, first, last) 
				context.restore();
				for(int kdx = first; kdx < last; kdx++)
				{
					mul(r0[kdx], interpolationUpTree[idx][2*kdx], coeffs[2*kdx + 1]);
					mul(r1[kdx], interpolationUpTree[idx][2*kdx + 1], coeffs[2*kdx]);
				}
				
			NTL_EXEC_RANGE_END
			
			
			coeffs.resize(interpolationUpTree[idx].size()/2);
			
			for(int kdx = 0; kdx < interpolationUpTree[idx].size()/2; kdx++)
			{
				add(coeffs[kdx], r0[kdx], r1[kdx]);
			}
			
// 			std::cout << "idx: " << idx << " " << coeffs.size() << " " << deg(coeffs[0]) << " ";
// 			t.Tick("");
		}
	}
	else
	{
		for(int idx = 0; idx < interpolationUpTree.size(); idx++)
		{
			Timer t;
			r0.resize(interpolationUpTree[idx].size()/2);
			r1.resize(interpolationUpTree[idx].size()/2);
			
			for(int kdx = 0; kdx < interpolationUpTree[idx].size()/2; kdx++)
			{
				if(VERBOSE) 
				{
					printPolynomial(interpolationUpTree[idx][2*kdx], "interpolationUpTree[idx][2*kdx]", 20);
					printPolynomial(coeffs[2*kdx + 1], "coeffs[2*kdx + 1]", 20);
					printPolynomial(interpolationUpTree[idx][2*kdx + 1], "interpolationUpTree[idx][2*kdx + 1]", 20);
					printPolynomial(coeffs[2*kdx], "coeffs[2*kdx]", 20);
				}
				
				mul(r0[kdx], interpolationUpTree[idx][2*kdx], coeffs[2*kdx + 1]);
				mul(r1[kdx], interpolationUpTree[idx][2*kdx + 1], coeffs[2*kdx]);
				
				if(VERBOSE) 
				{
					printPolynomial(r0[kdx], "r0", 20);
					printPolynomial(r1[kdx], "r1", 20);
				}
			}
			
			coeffs.resize(interpolationUpTree[idx].size()/2);
			for(int kdx = 0; kdx < interpolationUpTree[idx].size()/2; kdx++)
			{
				add(coeffs[kdx], r0[kdx], r1[kdx]);
			}
			
	// 		coeffs.resize(coeffs.size()/2);
			
// 			std::cout << "idx: " << idx << " " << coeffs.size() << " " << deg(coeffs[0]) << " ";
// 			t.Tick("");
		}
	}
}

void evaluatePolynomial(Vec<ZZ_p>& evals, const ZZ_pX& f, const Vec<ZZ_p>& vals)
{
	if(vals.length() == 1)
	{
		ZZ_pX g, h;
		g.SetLength(2);
		g[0] = 1;
		g[1] = -vals[0];
		
		rem(h, f, g);
		
		evals.SetLength(1);
		evals[0] = h[0];
		
		return;
	}
	else
	{
// 		Timer t;
// 		std::cout << "Uptree multiplication" << std::endl;
		
		std::vector<std::vector<ZZ_pX> > uptree, downtree;
		uptree_multiplication(uptree, vals, false);
// 		t.Tick("uptree_multiplication");
		
// 		std::cout << "Dividing down" << std::endl;
		
		downtree_division(downtree, uptree, f);
		
// 		t.Tick("downtree_division");

		for(int idx = 0; idx < downtree[0].size(); idx++)
		{
			evals[idx] = downtree[0][idx][0];
		}
	}
}

void interpolatePolynomial(ZZ_pX& f, const Vec<ZZ_p>& vals, const Vec<ZZ_p>& evals)
{
	Timer t;
	std::vector<std::vector<ZZ_pX> > interpolationUpTree, interpolationDownTree;
	int n = vals.length();
	
	// 1. Build uptree m(x) = (x - x0)...(x-x_(n-1))
	uptree_multiplication(interpolationUpTree, vals, true);
	t.Tick("uptree_multiplication");
	
	// 2. Compute the derivative. The input polynomial is at the bottom of the tree
	int k = interpolationUpTree.size();
	ZZ_pX derivativePolynomial;
	derivativePolynomial.SetLength(n);
	
	for(int idx = 0; idx < n; idx++)
	{
		derivativePolynomial[idx] = (interpolationUpTree[k-1][0][idx + 1]*(idx + 1));
	}
	
	t.Tick("derivative");
	
	if(1) printPolynomial(derivativePolynomial, "Derivative polynomial", 20);
	
	// 3. Evaluate m'(x) at x_0, ..., x_(n-1)
	
	// We don't need the full polynomial anymore
	interpolationUpTree.resize(k - 1);
	
	// Call downtree division to find the evaluation
	downtree_division(interpolationDownTree, interpolationUpTree, derivativePolynomial);
	
	t.Tick("downtree_division");
	
	// 4. Interpolate
	
	std::vector<ZZ_pX> coeffs;
	
	reconstructPolynomial(coeffs, evals, interpolationUpTree, interpolationDownTree);
	
	f = coeffs[0];
	
	t.Tick("Reconstruction");
}


void buildTree(PolynomialTree *tree, const ZZ_p *roots, int size)
{
	if(size == 1)
	{
		tree->val.SetLength(2);
		tree->val[1] = 1;
		tree->val[0] = -roots[0];
	}
	else
	{
		tree->left = new PolynomialTree();
		tree->right = new PolynomialTree();
		
		buildTree(tree->left, roots, size/2);
		buildTree(tree->right, (roots + size/2), (size - size/2));
		mul(tree->val, tree->left->val, tree->right->val);
		if(VERBOSE) tree->printTree();
	}
}

void downtreeDivisionReconstruct(ZZ_p *evals, PolynomialTree *tree, const ZZ_pX& f, int pos, int size)
{
	ZZ_pX g;
	rem(g, f, tree->val);
	
	if(tree->left != nullptr && tree->right != nullptr)
	{
		downtreeDivisionReconstruct(evals, tree->left, g, pos, size/2);
		downtreeDivisionReconstruct(evals, tree->right, g, pos + size/2, size - size/2);
		add(tree->interp, tree->left->interp*tree->right->val, tree->right->interp*tree->left->val);
	}
	else
	{
		tree->interp.SetLength(1);
		tree->interp[0] = evals[pos]*inv(g[0]);
		if(VERBOSE) std::cout << "g[" << pos << "]: " << g[0] << std::endl;
	}
}

void downtreeDivision(ZZ_p *evals, PolynomialTree *tree, const ZZ_pX& f, int pos, int size)
{
	ZZ_pX g;
	rem(g, f, tree->val);
	
	if(tree->left != nullptr)
	{
		downtreeDivision(evals, tree->left, g, pos, size/2);
	}
	else
	{
		evals[pos] = g[0];
	}
	
	
	if(tree->right != nullptr)
	{
		downtreeDivision(evals, tree->right, g, pos + size/2, size - size/2);
	}
	else
	{
		evals[pos] = g[0];
	}
}

void reconstruct(ZZ_p *evals, ZZ_p *derivativeEvals, PolynomialTree *tree, int pos, int size)
{
	if(tree->left != nullptr && tree->right != nullptr)
	{
		reconstruct(evals, derivativeEvals, tree->left, pos, size/2);
		reconstruct(evals, derivativeEvals, tree->right, pos + size/2, size - size/2);
		add(tree->interp, tree->left->interp*tree->right->val, tree->right->interp*tree->left->val);
	}
	else
	{
		tree->interp.SetLength(1);
		tree->interp[0] = inv(derivativeEvals[pos])*evals[pos];
		if(VERBOSE) std::cout << "g[" << pos << "]: " << evals[0] << std::endl;
	}
}

void evaluatePolynomial(ZZ_p *evals, const ZZ_pX& f, const ZZ_p *roots, int size)
{
	Timer t;
	PolynomialTree *tree = new PolynomialTree();
	buildTree(tree, roots, size);
	t.Tick("Build tree");
	downtreeDivision(evals, tree, f, 0, size);
	t.Tick("Evaluate");
}

void interpolatePolynomial(ZZ_pX& f, ZZ_p *roots, ZZ_p *evals, int size)
{
	PolynomialTree *tree = new PolynomialTree();
	buildTree(tree, roots, size);
	ZZ_pX derivativePolynomial;
	derivativePolynomial.SetLength(size);
	
	for(int idx = 0; idx < (size); idx++)
	{
		derivativePolynomial[idx] = (tree->val[idx + 1]*(idx + 1));
	}
	
	if(VERBOSE) printPolynomial(derivativePolynomial, "der poly", 20);
	
	downtreeDivisionReconstruct(evals, tree, derivativePolynomial, 0, size);
	
	f = tree->interp;
}

int test()
{  
// 	ZZ p = (ZZ(1) << 89) - 1;  // n = 2^24
	ZZ p = (ZZ(1) << 80) - 65; // n = 2^20
// 	ZZ p = (ZZ(1) << 72) - 93; // n = 2^16
// 	ZZ p = (ZZ(1) << 64) - 59; // n = 2^12
// 	ZZ p = (ZZ(1) << 56) - 5;  // n = 2^8
// 	ZZ p = ZZ(97);
	ZZ_p::init(p);
	
	SetNumThreads(1);
	
	int n = 100000; //(1 << 16) ;// 0x165555; //(1 << 20);
	std::cout << "n: " << n << std::endl;
// 	n = 6;
	
	Vec<ZZ_p> vals, evals;
	vals.SetLength(n);
	evals.SetLength(n);
	
	ZZ_pX f, g;
	f.SetLength(n);
	
	// Generate inputs
	generateInput(vals, f);
	
	printPolynomial(vals, "input", 10);
	
	Timer totalTime, t;
	
	ZZ_p *roots = new ZZ_p[n];
	ZZ_p *fr = new ZZ_p[n];
	
	for(int idx = 0; idx < n; idx++) {roots[idx] = vals[idx];}
	
	evaluatePolynomial(fr, f, roots, n);
	
	t.Tick("Evaluation 1");
	
	interpolatePolynomial(g, roots, fr, n);
	
	printPolynomial(g, "output", 10);
	
	t.Tick("Interpolation 1");
	
	return 0;
	
	int temp = n;
	int groupSize = 1;
	int location = 0;
	
	while(temp > 0)
	{
		if(temp & 0x1)
		{
// 			std::cout << "group size: " << groupSize << "\t " << " ";
			Vec<ZZ_p> subVals;
			Vec<ZZ_p> subEvals;
			
			subVals.SetLength(groupSize);
			subEvals.SetLength(groupSize);
			
			for(int idx = 0; idx < groupSize; idx++)
			{
				subVals[idx] = vals[idx + location];
			}
			
			evaluatePolynomial(subEvals, f, subVals);
			
			for(int idx = 0; idx < groupSize; idx++)
			{
				evals[idx + location] = subEvals[idx];
			}
			
			location += groupSize;
// 			t.Tick("");
		}
		
		temp = (temp >> 1);
		groupSize *= 2;
	}
	
	t.Tick("Evaluation 2");
	
// 	printPolynomial(evals, "evals", 20);
// 	return 0;
	
	std::cout << "---------------Polynomial Evaluation--------------" << std::endl;
	
	evaluatePolynomial(evals, f, vals);
	
	totalTime.Tick("Evaluation 3");
	
	std::cout << "---------------Polynomial interpolation--------------" << std::endl;
	
	interpolatePolynomial(g, vals, evals);
	
	totalTime.Tick("Interpolation 2");
	
	printPolynomial(g, "\n***interpolated polynomial***\n", 10);
	printPolynomial(f, "\n***original polynomial***\n", 10);
	
	for(int idx = 0; idx < n; idx++)
	{
		assert(f[idx] == g[idx]);
	}
	
	return 0;
}

