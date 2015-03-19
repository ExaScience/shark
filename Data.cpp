/*
 * Data.cpp
 *
 *  Created on: 16 févr. 2015
 *      Author: chakro23
 */

#include "Data.h"
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <set>
#include <map>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>

using namespace std;
using namespace boost;

Data::~Data()
{
	// TODO Auto-generated destructor stub
}

Data::Data()
{
}

#include <stdio.h>
#include <stdlib.h>
#include <sqlite3.h>

static int callback(void *data, int argc, char **argv, char **azColName)
{
	int i;

//   fprintf(stderr, "%s: ", (const char*)data);

   for(i=0; i<argc; i++)
   {
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }

   printf("\n");

   return 0;
}

void Data::readData(string problem, std::vector<double>& rowData, std::vector<int>& rowPtr, std::vector<int>& colInd, int& numRows, int& numCols,int& numRates, double& globalMean)
{
	// Table {row-id, col-id, rate}

	rates rating;

	std::multimap<int, rates> dataTable;
	std::multimap<int, rates>::iterator it;
	std::string line;

	double sum = 0;

	vector<string> data;

    if (problem.compare("chemo") == 0) //the data is seperated by a comma
    {
    	std::ifstream ifs;

    	ifs.open("../chembl_19_mf1/chembl-IC50-360targets.csv", std::ifstream::in);

    	while (std::getline(ifs, line))
    	{
    	 	split(data, line, is_any_of( "," ) );
        	int user = atoi(data[0].c_str());
        	int item = atoi(data[1].c_str());
        	double rate = atof(data[2].c_str());

        	numRates++;
        	sum += rate;

        	scaleDist.push_back(rate);

        	int item_size = itemIds.size();
        	int user_size = userIds.size();

        	if (userIds.count(user) <= 0) userIds.insert(std::pair<int,int>(user,user_size));
        	if (itemIds.count(item) <= 0) {itemIds.insert(std::pair<int,int>(item,item_size));}

        	rating.item = item;
        	rating.rate = rate;

        	dataTable.insert(std::pair<int,rates>(user,rating));
    	}

    	globalMean = sum / numRates;

    	ifs.close();
    }
    else if (problem.compare("MovieLens") == 0) //the data is seperated by a two double columns
    {
    	std::ifstream ifs;

    	ifs.open("../chembl_19_mf1/ratings.dat", std::ifstream::in);

    	while (std::getline(ifs, line))
    	{
    	 	split(data, line, is_any_of("::"));

    	 	int user = atoi(data[0].c_str());
        	int item = atoi(data[2].c_str());
        	double rate = atof(data[4].c_str());

        	numRates++;
        	sum += rate;

        	scaleDist.push_back(rate);

        	int item_size = itemIds.size();
        	int user_size = userIds.size();

        	if (userIds.count(user) <= 0) userIds.insert(std::pair<int,int>(user,user_size));
        	if (itemIds.count(item) <= 0) {itemIds.insert(std::pair<int,int>(item,item_size));}

        	rating.item = item;
        	rating.rate = rate;

        	dataTable.insert(std::pair<int,rates>(user,rating));
    	}

    	globalMean = sum / numRates;

    	ifs.close();
    }
    else if (problem.compare("Netflix") == 0) // the data is stored in database
	{
       #define DB "../chembl_19_mf1/netflix.sqlite"

       sqlite3 *dbfile;
	   sqlite3_stmt *statement;

	   if ( sqlite3_open(DB, &dbfile) != SQLITE_OK )
	   {
		   printf ("Error Database not opened successfully \n");
		   exit(0);
	   }

	   char *query = "select User, Movie, Rating FROM Ratings where user < 20000";

	   if ( sqlite3_prepare(dbfile, query, -1, &statement, 0 ) == SQLITE_OK )
	   {
	         int ctotal = sqlite3_column_count(statement);

	         int res = 0;

	         while ( 1 )
	         {
	             res = sqlite3_step(statement);

	             if ( res == SQLITE_ROW )
	             {
	                 int user = atoi((char*) sqlite3_column_text(statement, 0));
	                 int item = atoi((char*) sqlite3_column_text(statement, 1));
	                 double rate = atof((char*) sqlite3_column_text(statement, 2));

	                 numRates++;
	                 sum += rate;

	                 scaleDist.push_back(rate);

	                 int item_size = itemIds.size();
	                 int user_size = userIds.size();

	                 if (userIds.count(user) <= 0) userIds.insert(std::pair<int,int>(user,user_size));
	                 if (itemIds.count(item) <= 0) {itemIds.insert(std::pair<int,int>(item,item_size));}

	                 rating.item = item;
	                 rating.rate = rate;

	                 dataTable.insert(std::pair<int,rates>(user,rating));
	             }

	             globalMean = sum / numRates;

	             if ( res == SQLITE_DONE || res==SQLITE_ERROR) { break; }
	         }
	     }

	   sqlite3_close(dbfile);

       /*std::ifstream ifs;

		ifs.open("../chembl_19_mf1/netflix.csv", std::ifstream::in);

		while (std::getline(ifs, line))
		{
			split(data, line, is_any_of( "," ) );
			int user = atoi(data[0].c_str());
			int item = atoi(data[1].c_str());
			double rate = atof(data[2].c_str());

			numRates++;
			sum += rate;

			scaleDist.push_back(rate);

			int item_size = itemIds.size();
			int user_size = userIds.size();

			if (userIds.count(user) <= 0) userIds.insert(std::pair<int,int>(user,user_size));
			if (itemIds.count(item) <= 0) {itemIds.insert(std::pair<int,int>(item,item_size));}

			rating.item = item;
			rating.rate = rate;

			dataTable.insert(std::pair<int,rates>(user,rating));

			if (numRates > 50000 ) break;
	   }

		globalMean = sum / numRates;

		ifs.close();*/
    }

	numRows = userIds.size(), numCols = itemIds.size();
	numRates = scaleDist.size();

	std::sort(scaleDist.begin(), scaleDist.end());

	// if min-rate = 0.0, shift upper a scale
	double minRate = scaleDist[0];
	double epsilon = minRate == 0.0 ? scaleDist[1] - minRate : 0;

	if (epsilon > 0)
	{
		for (int i = 0, im = scaleDist.size(); i < im; i++)		scaleDist[i] += epsilon;

		for (int row = 0; row < numRows; row++)
			for (int col = 0; col < numCols; col++)
			 	if (dataTable.count(row) > 0 )
				{
					it = dataTable.find(row);

					if ((it->second.item) == col) //update the rating
						it->second.rate = it->second.rate  + epsilon;
				}
	}

	int nnz = dataTable.size();

	bool CRS = false; //si pas CRS c'st COO

	if (CRS) //CRS COMPRESSION FORMAT
	{
		rowPtr.resize(numRows + 1);
		colInd.resize(nnz);
		rowData.resize(nnz);

		int j = 0, cols = 0;

		//Read data in the CRS format

		std::pair <std::multimap<int, rates>::iterator, std::multimap<int, rates>::iterator> ret;

		for (int i = 1; i <= numRows; ++i)
		{
			ret = dataTable.equal_range(i - 1);

			for (std::multimap<int, rates>::iterator it = ret.first; it != ret.second; ++it)
			{
				colInd[j++] = it->second.item;
				cols++;

				if (cols >= numCols)
					printf("ERROR colInd[%d]=%d, which is not a valid column index \n", j, cols );

			}

			rowPtr[i] = cols;

			vector<int>::iterator iter = colInd.begin();

			std::sort(iter+rowPtr[i-1], iter+rowPtr[i]);
		}

		//TODO utiliser plutot le CCS parce que le nbre de colonnes est inférieur aux nbres de lignes donc autant avoir deux grandes matrices (inévitables) et avoir la troisième plus petite

		// set data

		for (it = dataTable.begin(); it != dataTable.end(); ++it)
		{
			int row = it->first;
			int col = it->second.item;
			double val = it->second.rate;
			vector<int>::iterator i = colInd.begin();

			i = find(i + rowPtr[row], i + rowPtr[row + 1], col);

			int nPosition = distance(colInd.begin(), i);

			if (nPosition >= 0 && colInd[nPosition] == col)
				  rowData[nPosition] = val;
			else
				printf("ERROR ENTRY %d,%d is not in the matrix structure \n", (row + 1) ,(col + 1));
		}
	}
	else //COO COMPRESSION FORMAT
	{
		rowPtr.resize(nnz);
		colInd.resize(nnz);
		rowData.resize(nnz);

		int j = 0, cols = 0;

		std::pair <std::multimap<int,rates>::iterator, std::multimap<int,rates>::iterator> ret;

		for (int i = 1; i <= numRows; ++i)
		{
			ret = dataTable.equal_range(i);

			for (std::multimap<int, rates>::iterator it = ret.first; it != ret.second; ++it)
			{
				rowPtr[j] = it->first;
				colInd[j]= it->second.item;
				rowData[j]= it->second.rate;

				j++;
			}
		}
	}

	// release memory of data table
	for (it=dataTable.begin(); it!=dataTable.end(); ++it)
	   	dataTable.erase(it);
}
