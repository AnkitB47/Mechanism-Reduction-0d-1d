import numpy as np
from dataclasses import dataclass
from typing import Callable, Iterable, List, Tuple

@dataclass
class GAOptions:
    population_size: int = 20
    generations: int = 10
    crossover_rate: float = 0.8
    mutation_rate: float = 0.1


def initialize_population(pop_size: int, genome_length: int) -> np.ndarray:
    return np.random.randint(0, 2, size=(pop_size, genome_length))


def fitness(population: np.ndarray, eval_fn: Callable[[np.ndarray], float]) -> np.ndarray:
    return np.array([eval_fn(ind) for ind in population])


def selection(population: np.ndarray, scores: np.ndarray) -> np.ndarray:
    probs = scores / scores.sum()
    idx = np.random.choice(len(population), size=len(population), p=probs)
    return population[idx]


def crossover(population: np.ndarray, rate: float) -> np.ndarray:
    next_pop = population.copy()
    for i in range(0, len(population), 2):
        if np.random.rand() < rate and i+1 < len(population):
            point = np.random.randint(1, population.shape[1])
            next_pop[i, point:], next_pop[i+1, point:] = population[i+1, point:].copy(), population[i, point:].copy()
    return next_pop


def mutate(population: np.ndarray, rate: float) -> np.ndarray:
    mutation = np.random.rand(*population.shape) < rate
    population[mutation] = 1 - population[mutation]
    return population


def run_ga(
    genome_length: int,
    eval_fn: Callable[[np.ndarray], float],
    options: GAOptions = GAOptions(),
    return_history: bool = False,
) -> np.ndarray | Tuple[np.ndarray, List[float]]:
    """Run a simple genetic algorithm.

    Parameters
    ----------
    genome_length:
        Length of the binary genome.
    eval_fn:
        Fitness evaluation function returning a scalar value.
    options:
        Control parameters for the GA.
    return_history:
        If ``True`` the function also returns a list of best scores per
        generation.

    Returns
    -------
    np.ndarray or Tuple[np.ndarray, List[float]]
        The best genome found and optionally the score history.
    """

    pop = initialize_population(options.population_size, genome_length)
    best = pop[0]
    best_score = -np.inf
    history: List[float] = []
    for _ in range(options.generations):
        scores = fitness(pop, eval_fn)
        if scores.max() > best_score:
            best_score = float(scores.max())
            best = pop[scores.argmax()]
        history.append(best_score)
        pop = selection(pop, scores)
        pop = crossover(pop, options.crossover_rate)
        pop = mutate(pop, options.mutation_rate)
    if return_history:
        return best, history
    return best
