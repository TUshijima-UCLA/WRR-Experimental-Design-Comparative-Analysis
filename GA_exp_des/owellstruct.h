#ifndef owellstruct_h
#define owellstruct_h
#include <vector>				//std::vector<size>


//#include <ga/ga.h>
//#include <ga/std_stream.h>
//#include <ga/GARealGenome.h>

//structure that contains a potential observation node, zone, and x and y coords
struct owell
{
	int node;
	int zone;	//real zone number the node resides in
	int zonecount;	//maps nodes to genes
	int count;
	int gene;
	double x;
	double y;
	bool highfreq;
	owell(): zone(0),x(0.0),y(0.0),highfreq(false) {};
};
//function to sort the potential owells by zone
struct by_zone
{
bool operator()(owell const &a, owell const &b)
{
	return a.zone < b.zone;
}
};
//function to sort owells by count
struct by_count_owell
{
bool operator()(owell const &a, owell const &b)
{
	return a.count > b.count;
}
};


//structure holding the robust obswells
struct robustowell
{
	int count;
	double minscore;
	std::vector<owell> optowells;
	std::vector<int> paramindex;
	std::vector<int> genome;	
	std::vector<double> maxscore;
	std::vector<double> worstcaseK;
	//GARealGenome &genome;
	//GARealGenome &worstgenome;
	bool highfreqnetwork;
	bool highfreqwells;
	bool robustdesign;
	robustowell(): highfreqnetwork(false), highfreqwells(false), robustdesign(false) {};
};
//function to sort the owells by number
struct by_count_network
{
bool operator()(robustowell const &a, robustowell const &b)
{
	return a.count > b.count;
}
};

#endif