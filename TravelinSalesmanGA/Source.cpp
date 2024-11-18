#include "crossover.h"

using namespace std;
using namespace std::chrono;

//Problem parameters
const int POP_SIZE = 16; //could change these to a rule based on how many cities there are in the path
const float CROSSOVER_PER = 0.5; //needs to be half or else we're fucked
const float MUTATION_PER = 0.5; //50% mutation rate
const int ELITISM = 2; //take the top 2 best solutions from each generations
const int MAX_GENERATIONS = 500; //how many generations we are doing, different gene pools

struct TSPProblemData {
	string name;
	string comment;
	string type;
	int dimension;
	string edge_weight_type;
	vector<City> initCities;
};

TSPProblemData readTSPFile(const string& filepath) {
	TSPProblemData data;
	ifstream infile(filepath);
	string line;

	if (!infile) {
		cerr << "Unable to open File" << filepath << endl;
		exit(1);
	}

	while (getline(infile, line)) {
		if (line.find("NAME") == 0) {
			data.name = line.substr(line.find(":") + 1);
		}
		else if (line.find("COMMENT") == 0) {
			data.comment = line.substr(line.find(":") + 1);
		}
		else if (line.find("TYPE") == 0) {
			data.type = line.substr(line.find(":") + 1);
		}
		else if (line.find("DIMENSION") == 0) {
			data.dimension = std::stoi(line.substr(line.find(":") + 1));
		}
		else if (line.find("EDGE_WEIGHT_TYPE") == 0) {
			data.edge_weight_type = line.substr(line.find(":") + 1);
		}
		else if (line == "NODE_COORD_SECTION") {
			string identifier;
			float x, y;
			while (infile >> identifier >> x >> y) {
				data.initCities.emplace_back(x, y, identifier);
			}
		}
	}
	return data;
}

//globally setting the crossover function before we start, saves time on branching in the simulation
using CrossoverFunc = void(*)(Trip&, Trip&, std::vector<Trip>&,int);
CrossoverFunc selectCrossoverFunction(const std::string& crossoverType) {
	if (crossoverType == "UX") {
		return uniformCrossover;
	}
	else if (crossoverType == "PMX") {
		return partiallyMappedCrossover;
	}
	else if (crossoverType == "SPX") {
		return spCrossover;
	}
	return nullptr; // Default case or handle error
}
//Note need to change num cities in crossover.h
int main(int argc, char* argv[]) {
	//Defaults
	string crossoverType = "SPX"; //SPX,PMX,UX...can optimize the branching with these
	string mutationType = "M"; //R (Scramble), S (Simple Swap), M (Moro Mutate)
	string selectionType = "newRWS"; //SUS (Stochastic Universal Sampling, RWS (Roulette Wheel Selection), LRS (Linear Rank Selection), newRWS
	string filePath = "./tsp/original10.tsp";
	int nn = 0;
	if (argc < 4) {
		cout << "No file arguemnts specificed" << endl;
		cout << "Arguments: filename -CrossoverType -MutationType -selectionType -NN (optional)";
		cout << "Defaults are SPX, M, and newRWS, and not using Nearest Neighbor" << endl;
	}
	else {//hack job, but seeing if it will work
		string filePath = "./tsp/" + string(argv[1]);
		string crossoverType = string(argv[2]);
		string mutationType = string(argv[3]);
		string selectionType = string(argv[4]);
		if (argc == 6) {
			if (argv[5] == "NN") {
				int nn = 1; //not sure about this
			}
		}
	}
	srand(time(NULL));
	int run = 1;
	CrossoverFunc crossoverFunction = selectCrossoverFunction(crossoverType);
	cout << "Using crossover function: " << crossoverType << endl;
	TSPProblemData data = readTSPFile(filePath);
	cout << data.name << endl;
	cout << data.comment << endl;
	cout << data.type << endl;
	cout << data.dimension << endl;
	int numCities = data.dimension;
	cout << "Number of cities: " << numCities << endl;
	cout << data.edge_weight_type << endl;
	vector<City> initCities = data.initCities;
	int mutationLength  = int(floor(0.2 * numCities));
	int numSwaps = int(floor(0.3*numCities)); //arbitrary rules, can change these

	if (run == 0) {
		Trip NNTrip = NearestNeighbor(initCities, numCities);
		NNTrip.printPath();
		NNTrip.printPathLength();
		return 0;
	}
	else {
		auto start = high_resolution_clock::now();

		//Add nearest Neighbor clause - need to write nearest neighbor
		vector<Trip> genePool;
		//nn == 1 means we are using nearest neighbor
		if (nn == 1) {
			auto startNN = high_resolution_clock::now();
			Trip NNTrip = NearestNeighbor(initCities, numCities);
			auto stopNN = high_resolution_clock::now();
			auto durationNN = duration_cast<microseconds>(stopNN - startNN);
			cout << endl << "Time taken by Nearest Neighbor function: "
				<< durationNN.count() << " microseconds" << endl;
			cout << "Original Nearest Neighbors" << endl;
			NNTrip.printPath();
			NNTrip.printPathLength();
			for (int i = 0; i < POP_SIZE; i++) {
				genePool.push_back(NNTrip);
			}
		}
		else {
			int picks[] = { 0 };
			for (int i = 0; i < numCities; i++) picks[i] = i; //populating picks

			for (int i = 0; i < POP_SIZE; i++) {
				Trip newTrip;
				vector<int> picks;
				for (int i = 0; i < numCities; i++) { picks.push_back(i); };
				vector<int> newPicks;
				newPicks.insert(newPicks.begin(), picks.begin(), picks.end());
				int n = numCities;

				for (int i = 0; i < numCities; i++) {
					int randIndex = rand() % n;
					int numToAdd = newPicks[randIndex];
					newPicks[randIndex] = newPicks[n - 1];
					n--;
					newTrip.addCity(initCities[numToAdd]);
				}
				newTrip.calcPathLength();
				genePool.push_back(newTrip);
			}
		}
		/*
		cout << "Genes: " << endl;
		for (auto& gene : genePool) {
			gene.printPath();
			gene.printPathLength();
			cout << endl;
		}
		*/
		int mutations = 0;//tracking mutations
		vector<Trip> newGen;
		for (int p = 0; p <= MAX_GENERATIONS; p++) {
			//roulette wheel probability
			vector<Trip> parents;
			if (selectionType == "RWS") {
				RWS(genePool, parents, POP_SIZE,CROSSOVER_PER);
			}
			else if(selectionType == "SUS") {
				SUSSelection(genePool, parents, POP_SIZE);
			}
			else if (selectionType == "newRWS") {
				newRWSSelection(genePool, parents, POP_SIZE);
			}
			else if (selectionType == "LRS") {
				linearRankSelection(genePool, parents, 2);
			}

			vector<Trip> children;
			for(size_t i = 0;i+1<parents.size();i+=2){
				crossoverFunction(parents[i], parents[i + 1],children,numCities);
			}
			
			//mutation  - swapping cities in a path - 20% mutation chance per gene in children pool - introducing new genes essentially
			//I could try and parallelize this
			for (size_t i = 0; i < children.size(); i++) {
				float mutateThreshold = genRandom();
				if (mutateThreshold > (1 - MUTATION_PER)) {
					mutations++;
					if (mutationType == "S") {
						mutate(children[i], numCities);
					}
					else if (mutationType == "R") {
						scrambleMutate(children[i], mutationLength, numCities);
					}
					else if (mutationType == "M") {
						moroMutate(children[i], numSwaps, numCities);
					}
				}
			}

			//now we add the children to the current population - sort - then move on
			//1. Grab top 2 from original gene pool
			//2. Put all parents and children in the same vector - add children to parent vector
			//3. sort the parent vector, grab enough to fill enough for a population
			sort(genePool.begin(), genePool.end(), comparePaths); //sort the original
			//1.
			for (auto i = genePool.begin(); i != genePool.begin() + ELITISM; ++i) {
				newGen.push_back(*i);
			}
			//2.
			parents.insert(parents.end(), children.begin(), children.end());
			sort(parents.begin(), parents.end(), comparePaths); //sort the original
			//3.
			for (auto i = parents.begin(); i != parents.end()-ELITISM; ++i) {
				newGen.push_back(*i);
			}

			genePool.clear();
			genePool = newGen;
			newGen.clear();
		}

		cout << endl << "Best Solution: " << endl;
		genePool[0].printPath();
		genePool[0].printPathLength();
		cout << endl;
		cout << "Num mutations: " << mutations << endl;
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<seconds>(stop - start);
		cout << endl << "Time taken by function: "
			<< duration.count() << " Seconds" << endl;
		return 0;
	}
	
}