/*
 * Data.h
 *
 *  Created on: 16 févr. 2015
 *      Author: chakro23
 */

#ifndef SHARK_DATA_H_
#define SHARK_DATA_H_
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <list>
#include <set>

using namespace std;

/**
 * A data access class to a data file
 *
 */

class Data {
public:

	class rates{
	public:
		int item;
		double rate;
	};
	// path to data file
	string dataPath;

	// store data as {user/item rate} matrix
	vector<double*> rateMatrix;
	std::map<int,int> userIds;
	std::map<int,int> itemIds;

	// scale distribution
	vector<double> scaleDist;

	void readData(string problem, std::vector<double>& rowData, std::vector<int>& rowPtr, std::vector<int>& colInd, int& numRows, int& numCols, int& numRates, double& globalMean);
	void copyCRS(vector<double*> data, vector<int*> ptr, vector<int*> idx);
	Data();
	virtual ~Data();
};

#endif /* SHARK_DATA_H_ */
