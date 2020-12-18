class Edge
{
	int32_t id_u;
	int32_t id_v;
	uint8_t rating;
	std::vector<int64_t> profile_u;
	std::vector<int64_t> profile_v;
	bool isReal;
	uint64_t aMac_u;
	uint64_t aMac_v;
	uint64_t aMac_data;
	uint64_t bMac_data;
};

