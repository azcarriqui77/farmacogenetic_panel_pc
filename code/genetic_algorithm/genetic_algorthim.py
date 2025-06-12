import numpy as np
import pandas as pd
import random
from fitness import evaluate_fitness

class GeneticAlgorithm:
    def __init__(self, population_size, generations, crossover_rate, mutation_rate):
        super().__init__()
        self.population_size = population_size
        self.generations = generations
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate

    # Generate an initial random population, where each individual's genetic code is a binary vector of size feature_count.  
    # feature_count is the number of features or variables in the data DataFrame.  
    # Default random initialization: 70% probability of being 0 (feature not selected) and 30% probability of being 1 (selected).  
    def initialize_population(self, feature_count, probability=[0.7, 0.3]):
        return [np.random.choice([0, 1], size=(feature_count,), p=probability) for _ in range(self.population_size)]

    # Crossover operator. It only occurs if random.random() is less than crossover_rate.  
    # A random point in the vector is chosen, and a part of each parent is combined.  
    def crossover(self, parent1, parent2):
        if random.random() < self.crossover_rate:
            point = random.randint(1, len(parent1) - 1)
            return np.concatenate((parent1[:point], parent2[point:]))
        else:
            return parent1

    # Mutation operator. Each gene of the individual is traversed, and with a probability of MUTATION_RATE, its value is flipped (0 to 1 and vice versa).  
    def mutate(self, individual):
        for i in range(len(individual)):
            if random.random() < self.mutation_rate:
                individual[i] = 1 - individual[i]
        return individual
    
    # Execution of the genetic algorithm on the dataset
    def run(self, data, target, n_trees, algorithm, probability=[0.7, 0.3]):

        df = pd.DataFrame(data)
        feature_count = df.shape[1]
        population = self.initialize_population(feature_count, probability)


        for generation in range(self.generations):
            fitness_scores = [(individual, evaluate_fitness(individual, data, target, n_trees, algorithm)) for individual in population]
            fitness_scores.sort(key=lambda x: x[1], reverse=True)
            new_population = [fitness_scores[0][0]]  # Elitism
            
            while len(new_population) < self.population_size:
                parent1, parent2 = random.sample([ind for ind, _ in fitness_scores[:2*self.population_size//3]], 2)
                offspring = self.mutate(self.crossover(parent1, parent2))
                new_population.append(offspring)
            
            population = new_population

        return (fitness_scores[0][0], fitness_scores[0][1])