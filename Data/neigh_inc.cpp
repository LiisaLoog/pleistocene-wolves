//-----------------------------------------------------------------------------------
// Graph defining wolf regions
//
// 0: Europe
// 1: Central_North_Eurasia
// 2: Beringia
// 3: Middle_East
// 4: South_East_Asia
// 5: Arctic_North_America
// 6: North_America
//
//-----------------------------------------------------------------------------------

const int nDemes = 7;
const int nNeigh = 16;

int neigh[nNeigh] = {
	1,          // 0
	0, 2, 3, 4, // 1
	1, 4, 5, 6, // 2
	1,          // 3
	1, 2,       // 4
	2, 6,       // 5
	2, 5        // 6
};

const int nNeighs[nDemes] = {
	1,		// Deme 0: 1
	4,		// Deme 1: 5
	4,		// Deme 2: 9
	1,		// Deme 3: 10
	2,		// Deme 4: 12
	2,		// Deme 5: 14
	2		// Deme 6: 16
};

//                                 0  1  2  3  4   5   6,  terminal
const int neighStart[nDemes+1] = { 0, 1, 5, 9, 10, 12, 14, 16 };
