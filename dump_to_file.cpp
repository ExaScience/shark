/*
 * Copyright (c) 2010-2013, Vrije Universiteit Brussel.
 * Copyright (c) 2014-2015, imec
 * All rights reserved.
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include <getopt.h>
#include <shark.hpp>

using namespace std;
using namespace shark;
using namespace shark::types1d;


void test_dump() {
	int n = 511;

	Domain d(world(), {{n}});
	auto src = coord_val<0>(d);
	GlobalArrayI ga(d);
	ga = src;
	ga.dump(ga.region(), "ga.dump");
}

int main(int argc, char **argv) {
	Init(&argc, &argv);
	SetupThreads();

	test_dump();

	Finalize();

	return 0;
}
