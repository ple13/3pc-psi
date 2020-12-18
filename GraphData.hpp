#ifndef GRAPH_DATA_H__
#define GRAPH_DATA_H__

#include <vector>

#define dimension 10

struct Secret_Edge 
{
	uint64_t id_u;  	// 64 bit, [1..2^64]  ~ 4 byte
	uint64_t id_v;  	// 64 bit, [1..2^64]  ~ 4 byte
	uint64_t profile_u[dimension]; 	// 10*64 bit [0..1]
	uint64_t profile_v[dimension];	// 10*64 bit, [0..1]
	uint64_t rating; 	// 8 bit, [0..5]     	~ 1 byte
	uint64_t isReal; 	// 8 bit, [true,false]  ~ 1 byte
};

typedef Secret_Edge SecretEdgeMAC;

// secret shared values
struct Half_Secret_Edge {
	uint64_t profile_u[dimension];
	uint64_t profile_v[dimension];
	uint64_t rating;
	uint64_t isReal;
	SecretEdgeMAC bMACs;
};

struct Secret_Node 
{
	uint64_t vertexID;  // 16 bit, [1..2^16-1]
	uint64_t Profile[dimension]; // 10*64 bit [0..1]
	uint64_t newProfile[dimension]; // 10*64 bit, [0..1]
	
	std::vector<Half_Secret_Edge> halfEdge;
};

// MACed secret_shared values
struct MAC_MFNode 
{
	uint32_t vertexID;  // 32 bit, [1..2^32]
	std::vector<uint64_t> Profile; // 32 bit [0..1]
	std::vector<uint64_t> newProfile; // 32 bit [0..1]
	std::vector<uint64_t> nodeData; 
};
struct MAC_MFEdge 
{
	std::bitset<40> aMAC;
	std::bitset<40> bMAC_u;
	std::bitset<40> bMAC_v;
	std::bitset<40> bMAC_data;
};

#endif
