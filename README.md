## The Genetic Algorithm and the Traveling Salesman Problem
#### How to Run on Linux

commands:
make
./main fileName Crossover Mutation Selection maxGenerations PopulationSize NN

- Crossover can be SPX (Single Point Crossover), UX (Uniform Crossover), or ERX (Edge Recombination Operator)

- Mutation can be S (Swap), R (Scramble), or M (Moro)

- Selection can be RWS or newRWS (Roulette Wheel), LRS (Linear Rank Selection), or  SUSN (Stochastic Universal Sampling)

- Max Generations dictates how many populations the algorithm will create

- Population Size dictates how many individuals are in a each generation

- NN can be left blank or if you put NN, it will run nearest neighbors before it starts with generations.

Note: you may need to remove main to recompile the program

#### Example of running:
./main att48.tsp SPX M LRS 500 16 NN
Here we are using att48.tsp as the problem, Single Point Crossover,Moro Mutate, Linear Rank Selection, 500 generations, 16 individuals, and asking to use Nearest Neighbor before we start

#### Another example is:
./main pr1002.tsp UX R SUSN 1000 32
Here are using pr1002.tsp as the problem, Universal Crossover, Scramble Mutate, and Stochastic Universal Sampling, 1000 genetations and 32 individuals
