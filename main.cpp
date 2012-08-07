
#include <iostream>
#include <shark/group.hpp>
#include <shark/globals.hpp>

using namespace std;
using namespace shark;

int main(int argc, char* argv[]) {
	Init(&argc, &argv);
	SetupThreads();

	cout << "Hello world from " << world().procid << " of " << world().nprocs << endl;

	Finalize();
}

