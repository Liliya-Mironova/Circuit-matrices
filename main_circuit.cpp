#include "matrix.hpp"

#define N 8 // number of input files

using namespace Parser;

int main() {

	const char* filename[N] = {"input1.txt", "input2.txt", "input3.txt", 
								"input4.txt", "input5.txt", "input6.txt",
								"input7.txt", "input8.txt"};// array of input files

	Circuit MyCircuit1 = parse(filename[0]);
	Circuit MyCircuit2 = parse(filename[1]);
	Circuit MyCircuit3 = parse(filename[2]);
	Circuit MyCircuit4 = parse(filename[3]);
	Circuit MyCircuit5 = parse(filename[4]);
	Circuit MyCircuit6 = parse(filename[5]);
	Circuit MyCircuit7 = parse(filename[6]);
	Circuit MyCircuit8 = parse(filename[7]);
	
	int i = 1;
	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit1.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit2.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit3.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit4.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit5.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit6.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit7.solve();

	printf("CIRCUIT %d:\n\n", i++);
	MyCircuit8.solve();
	MyCircuit8.delete0(); // delete zero connections
	printf("Zero connections deleted\n\n");
	MyCircuit8.solve();

	return 0;
}