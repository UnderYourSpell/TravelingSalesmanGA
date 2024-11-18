## The Genetic Algorithm and the Traveling Salesman Problem
#### How to Run on Linux
use the branch linux_base

commands:
make
./main fileName Crossover Mutation Selection NN

- Crossover can be either SPX (Single Point Crossover) or UX (Uniform Crossover)

- Mutation can be S (Swap), R (Scramble), or M (Moro)

- Selection can be RWS or newRWS (Roulette Wheel), LRS (Linear Rank Selection), or  SUS (Stochastic Universal Sampling)

- NN can be left blank or if you put NN, it will run nearest neighbors before it starts with generations.

#### Example of running:
./main att48.tsp SPX M LRS NN
Here we are using att48.tsp as the problem, Single Point Crossover,Moro Mutate, Linear Rank Selection, and asking to use Nearest Neighbor before we start

#### Another example is:
./main pr1002.tsp UX R SUS
Here are using pr1002.tsp as the problem, Universal Crossover, Scramble Mutate, and Stochastic Universal Sampling
