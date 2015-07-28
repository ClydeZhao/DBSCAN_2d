#include "pregel_pre_data.h"

int main(int argc, char* argv[]){
	init_workers();
	pregel_pre_data("/dbscanToy2", "/2dtestPreOut", 9);
	worker_finalize();
	return 0;
}
