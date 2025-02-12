import subprocess
from rpy2.robjects import r
from rc_data_class import RcData
from rc_folds_class  import rcFolds
from rc_pred_ann_model import PredAnnModel
from rc_individual_fold import RcFoldForANN
import numpy as np
import random
import pandas as pd
import random
import matplotlib.pyplot as plt
import pickle
from multiprocessing import Process
import time
import os

# Genetic Algorithm Parameters
POP_SIZE = 10       # Population size
N_GENERATIONS = 5    # Number of generations
TOURNAMENT_SIZE = 3  # Tournament selection (k=3)
CROSSOVER_RATE = 0.8 # Probability of crossover
MUTATION_RATE = 0.1 # Mutation probability per gene
ELITE_PERCENTAGE = 0.05 # Top 5% preserved

#Other parameters
p_value = 0.2
split_train = True
folds = 5
num_epochs = 50
dna_dict = {}  # Empty dictionary

# making a data frame to keep track of GA progress
column_names = [f'auc_{i+1}' for i in range(POP_SIZE)]
# Initialize an empty DataFrame with columns
ga_df = pd.DataFrame(columns=column_names)

# Global variables to store the results of each fold
first_fold_results = None
second_fold_results = None
third_fold_results = None
fourth_fold_results = None
fifth_fold_results = None

class TimerLogger:
    def __init__(self, log_file="execution_times.txt"):
        self.start_times = {}
        self.log_file = log_file
        
        # Clear the file at the start
        with open(self.log_file, "w") as f:
            f.write("Execution Timing Log\n")
            f.write("="*30 + "\n")

    def start(self, label):
        """Start a timer for a given label."""
        self.start_times[label] = time.time()

    def stop(self, label):
        """Stop the timer and log the elapsed time."""
        if label not in self.start_times:
            print(f"Warning: No timer started for '{label}'")
            return

        elapsed_time = time.time() - self.start_times.pop(label)
        log_message = f"{label}: {elapsed_time:.2f} seconds\n"

        # Print to console
        print(log_message.strip())

        # Write to file
        with open(self.log_file, "a") as f:
            f.write(log_message)

def get_genes_list(p_thresh, split_train):
    # Define the R script path
    r_script = "rc_get_diff_genes.r"
    
    # Build the command to run the R script
    command = ["Rscript", r_script, str(p_thresh), str(split_train)]
    
    result = subprocess.run(command, capture_output=True, text=True)
    
    # Check if the R script ran successfully
    if result.returncode == 0:
        print("R script executed successfully.")
    
        # Read the generated file

        # Get the directory of the current script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Define the file path relative to the script
        rds_path = os.path.join(script_dir, "ann_gene_set.rds")
        
        current_genes = r.readRDS(rds_path)
        print(len(current_genes))
        
    else:
        print("Error in R script execution:")
        print(result.stderr)

    return(current_genes.tolist())

def generate_individual():
    """Creates a binary chromosome for feature selection."""
    return [random.randint(0, 1) for _ in range(len(current_genes))]

def initiate_fold(current_folds,genes_list,fold,fold_name,results):
    genes_list = genes_list.tolist()
    current_fold = RcFoldForANN(current_folds,0)
    current_model = PredAnnModel(current_fold, genes_list, num_epochs=num_epochs)
    results = current_model.test_auc_list  # Store result

    # Write result to a file with a unique name for each fold
    result_filename = f"results_{fold_name}.txt"
    with open(result_filename, 'w') as f:
        f.write(str(results))  # Store the result (convert to string if necessary)

def evaluate_fitness(individual,gene_list,input_data,count):
    """Evaluates the fitness of an individual based on the average test auc value across folds."""
    start_time = time.time()
    # selected_features = [s for s, m in zip(gene_list, individual) if m]  #application of a binary mask to the genes list
    individual = np.array(individual, dtype=bool)  # Ensure it's a boolean array
    gene_list = np.array(gene_list)  # Convert to NumPy array if it's a list
    selected_features = gene_list[individual]
    
    if len(selected_features) == 0:
        return 0  # Prevent division by zero when no features are selected
    current_aucs = []
    current_folds = rcFolds(input_data,folds)
    current_member = count + 1
    print(f"Currently training, population member {current_member}")

    results = {}

    print(len(selected_features))
    if __name__ == '__main__':
        t1 = Process(target=initiate_fold, args=(current_folds, selected_features, 0, 'first', results))
        t2 = Process(target=initiate_fold, args=(current_folds, selected_features, 1, 'second', results))
        t3 = Process(target=initiate_fold, args=(current_folds, selected_features, 2, 'third', results))
        t4 = Process(target=initiate_fold, args=(current_folds, selected_features, 3, 'fourth', results))
        t5 = Process(target=initiate_fold, args=(current_folds, selected_features, 4, 'fifth', results))
        
        t1.start()
        t2.start()
        t3.start()
        t4.start()
        t5.start()
        
        t1.join()
        t2.join()
        t3.join()
        t4.join()
        t5.join()

     # Now read the results from the files after all processes are finished
    results = {}
    for fold_name in ['first', 'second', 'third', 'fourth', 'fifth']:
        result_filename = f"results_{fold_name}.txt"
        if os.path.exists(result_filename):
            with open(result_filename, 'r') as f:
                results[fold_name] = eval(f.read())  # Read and evaluate the result (could be a list or value)
        else:
            print(f"Warning: {result_filename} does not exist!")

    # Ensure all results are available before calculating score
    if any(val is None for val in results.values()):
        print("Warning: Some results are None. Check if subprocesses executed correctly.")
        return 0  # If any result is None, return 0 or handle it accordingly

    # Calculate the score using the results from all folds
    score = np.mean([max(results['first']), max(results['second']), max(results['third']),
                     max(results['fourth']), max(results['fifth'])])
    print(f"Score: {score}")

    end_time = time.time()
    print(f"Total time: {end_time - start_time} seconds")
    return score  # Higher auc average = better fitness

def tournament_selection(population, fitness_scores):
    """Selects a parent using tournament selection (k=3)."""
    competitors = random.sample(list(enumerate(fitness_scores)), TOURNAMENT_SIZE)
    best = max(competitors, key=lambda x: x[1])  # Select individual with best fitness
    return population[best[0]]


def uniform_crossover(parent1, parent2):
    """Performs uniform crossover (each gene has 50% chance of swapping)."""
    if random.random() < CROSSOVER_RATE:
        SWAP_PROBABILITY = 0.2
        child1 = [p1 if random.random() > SWAP_PROBABILITY else p2 for p1, p2 in zip(parent1, parent2)]
        child2 = [p2 if random.random() > SWAP_PROBABILITY else p1 for p1, p2 in zip(parent1, parent2)]
        return child1, child2
    return parent1[:], parent2[:]  # No crossover, children are copies

def mutate(individual):
    """Mutates an individual by flipping bits with a small probability."""
    return [1 - gene if random.random() < MUTATION_RATE else gene for gene in individual]

def plot_row_averages(df):
    """
    Plots the average of each row in the given DataFrame.
    
    Parameters:
    df (pd.DataFrame): Input DataFrame containing numerical values.
    """
    row_averages = df.mean(axis=1)  # Compute the average across each row
    
    plt.figure(figsize=(10, 5))
    plt.plot(row_averages, marker='o', linestyle='-', color='b', label='Row Averages')
    
    plt.xlabel("Generation")
    plt.ylabel("Average AUC Values")
    plt.title("Average of Each Row in DataFrame")
    plt.legend()
    plt.grid(True)
    
    plt.show()

timer = TimerLogger()
    
# establish the dataset object

timer.start("Step 1: Data preparation")
current_data = RcData()
timer.stop("Step 1: Data preparation")


# grab the intial feature set
timer.start("Step 2: Grab differentially expressed genes")
current_genes = get_genes_list(p_value, split_train)
timer.stop("Step 2: Grab differentially expressed genes")


timer.start("Step 3: Establish population")

# making a data frame to keep track of GA progress
column_names = [f'auc_{i+1}' for i in range(POP_SIZE)]

# Initialize an empty DataFrame with columns
ga_df = pd.DataFrame(columns=column_names)

# Initialize population
population = [generate_individual() for _ in range(POP_SIZE)]
dna_dict[1] = population
timer.stop("Step 3: Establish population")


# Evaluate initial fitness

timer.start("Step 4: Evaluate first generation")
fitness_scores = [evaluate_fitness(ind, current_genes, current_data, count) for count, ind in enumerate(population)]


best_fitness = max(fitness_scores)
print(f"Generation {1}, Best Accuracy: {best_fitness:.4f}")


ga_df.loc[len(ga_df)] = fitness_scores
timer.stop("Step 4: Evaluate first generation")



timer.start(f"Step 5: Evaluate {N_GENERATIONS} generations")


for gen in range(N_GENERATIONS):
    # Select the top individuals (elitism)
    elite_count = int(ELITE_PERCENTAGE * POP_SIZE)
    elites = [population[i] for i in np.argsort(fitness_scores)[-elite_count:]]  # Keep best individuals

    # Generate next generation
    new_population = elites[:]  # Start with elites

    while len(new_population) < POP_SIZE:
        # Select parents using tournament selection
        parent1 = tournament_selection(population, fitness_scores)
        parent2 = tournament_selection(population, fitness_scores)

        # Crossover to generate children
        child1, child2 = uniform_crossover(parent1, parent2)

        # Apply mutation
        child1 = mutate(child1)
        child2 = mutate(child2)

        # Add to new population (ensure we don't exceed population size)
        new_population.append(child1)
        if len(new_population) < POP_SIZE:
            new_population.append(child2)

    # Update population and fitness scores
    population = new_population
    current_generation = gen + 2
    dna_dict[current_generation] = new_population
    fitness_scores = [evaluate_fitness(ind, current_genes, current_data, count) for count, ind in enumerate(population)]
    ga_df.loc[len(ga_df)] = fitness_scores

    pickle.dump(dna_dict, open("dna_dict.pkl", "wb"))
    ga_df.to_pickle("ga_df.pkl")

    # Print best result every 10 generations
    if gen % 1 == 0:
        best_fitness = max(fitness_scores)
        print(f"Generation {current_generation}, Best Accuracy: {best_fitness:.4f}")
        
timer.stop(f"Step 5: Evaluate {N_GENERATIONS} generations")

