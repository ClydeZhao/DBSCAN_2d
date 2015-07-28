#include "pregel_dbscan.h"

int main(int argc, char* argv[]){
	init_workers();
	pregel_dbscan("/2dtestPreOut", "/2dtestOut", 9, 5, 0.001);
	worker_finalize();
	return 0;
}