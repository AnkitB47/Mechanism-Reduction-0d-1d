"""GA-GNN Mechanism Reduction Pipeline.

Integrates genetic algorithm search with graph neural network scoring
for automated mechanism reduction while maintaining accuracy.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass
import yaml
import random
import copy

from mechanism.loader import Mechanism
from reactor.plug_flow import PlugFlowReactor, ReactorGeometry, HeatTransferConfig, InletConditions
from validation.hp_pox_validator import HP_POXValidator


@dataclass
class GAParameters:
    """Genetic algorithm parameters."""
    
    population_size: int = 50
    generations: int = 100
    mutation_rate: float = 0.1
    crossover_rate: float = 0.8
    elite_size: int = 5
    tournament_size: int = 5


@dataclass
class GNNParameters:
    """Graph neural network parameters."""
    
    hidden_dim: int = 64
    num_layers: int = 3
    learning_rate: float = 0.001
    epochs: int = 100
    batch_size: int = 32


@dataclass
class ReductionTarget:
    """Mechanism reduction target."""
    
    max_species: int
    max_error_pct: float
    validation_cases: List[str]
    priority_species: List[str] = None  # Species to preserve


@dataclass
class Individual:
    """Individual in genetic algorithm population."""
    
    species_mask: np.ndarray  # Boolean array for species selection
    fitness: float = 0.0
    error_metrics: Dict[str, float] = None
    
    def __post_init__(self):
        if self.error_metrics is None:
            self.error_metrics = {}


class GAGNNPipeline:
    """GA-GNN mechanism reduction pipeline."""
    
    def __init__(
        self,
        mechanism_file: str,
        validation_cases: List[str],
        reduction_target: ReductionTarget,
        ga_params: GAParameters = None,
        gnn_params: GNNParameters = None
    ):
        """Initialize GA-GNN pipeline.
        
        Args:
            mechanism_file: Path to detailed mechanism file
            validation_cases: List of validation cases to use
            reduction_target: Target reduction parameters
            ga_params: Genetic algorithm parameters
            gnn_params: GNN parameters
        """
        self.mechanism_file = mechanism_file
        self.validation_cases = validation_cases
        self.reduction_target = reduction_target
        self.ga_params = ga_params or GAParameters()
        self.gnn_params = gnn_params or GNNParameters()
        
        # Load mechanism
        self.mechanism = Mechanism(mechanism_file)
        self.n_species = len(self.mechanism.species_names)
        
        # Initialize validator
        self.validator = HP_POXValidator()
        
        # Initialize GNN (placeholder - would use actual GNN implementation)
        self.gnn = None  # Placeholder for GNN model
        
        # Results storage
        self.population_history = []
        self.best_individuals = []
        self.convergence_data = []
        
    def run_reduction(self, output_dir: str = "results") -> Dict:
        """Run complete GA-GNN reduction pipeline.
        
        Args:
            output_dir: Output directory for results
            
        Returns:
            Reduction results
        """
        
        print("Starting GA-GNN mechanism reduction pipeline...")
        print(f"Target: {self.reduction_target.max_species} species, "
              f"max error: {self.reduction_target.max_error_pct}%")
        
        # Initialize population
        population = self._initialize_population()
        
        # Run genetic algorithm
        for generation in range(self.ga_params.generations):
            print(f"Generation {generation + 1}/{self.ga_params.generations}")
            
            # Evaluate fitness
            self._evaluate_population(population, generation)
            
            # Store best individuals
            best_individual = max(population, key=lambda x: x.fitness)
            self.best_individuals.append(copy.deepcopy(best_individual))
            
            # Store convergence data
            self.convergence_data.append({
                "generation": generation,
                "best_fitness": best_individual.fitness,
                "avg_fitness": np.mean([ind.fitness for ind in population]),
                "n_species": np.sum(best_individual.species_mask)
            })
            
            # Check convergence
            if self._check_convergence():
                print("Convergence achieved!")
                break
            
            # Create next generation
            population = self._create_next_generation(population)
        
        # Find best individual
        best_individual = max(self.best_individuals, key=lambda x: x.fitness)
        
        # Create reduced mechanism
        reduced_mechanism = self._create_reduced_mechanism(best_individual)
        
        # Validate reduced mechanism
        validation_results = self._validate_reduced_mechanism(reduced_mechanism)
        
        # Save results
        self._save_results(best_individual, reduced_mechanism, validation_results, output_dir)
        
        return {
            "best_individual": best_individual,
            "reduced_mechanism": reduced_mechanism,
            "validation_results": validation_results,
            "convergence_data": self.convergence_data
        }
    
    def _initialize_population(self) -> List[Individual]:
        """Initialize genetic algorithm population."""
        
        population = []
        
        for _ in range(self.ga_params.population_size):
            # Create random species mask
            species_mask = np.zeros(self.n_species, dtype=bool)
            
            # Ensure minimum number of species
            min_species = max(10, self.reduction_target.max_species // 2)
            n_selected = random.randint(min_species, self.reduction_target.max_species)
            
            # Select species randomly
            selected_indices = random.sample(range(self.n_species), n_selected)
            species_mask[selected_indices] = True
            
            # Ensure priority species are included
            if self.reduction_target.priority_species:
                for species in self.reduction_target.priority_species:
                    if species in self.mechanism.species_names:
                        idx = self.mechanism.species_names.index(species)
                        species_mask[idx] = True
            
            population.append(Individual(species_mask=species_mask))
        
        return population
    
    def _evaluate_population(self, population: List[Individual], generation: int):
        """Evaluate fitness of population."""
        
        for i, individual in enumerate(population):
            print(f"  Evaluating individual {i+1}/{len(population)}")
            
            try:
                # Create reduced mechanism
                reduced_mechanism = self._create_reduced_mechanism(individual)
                
                # Evaluate fitness
                fitness, error_metrics = self._evaluate_fitness(reduced_mechanism)
                
                individual.fitness = fitness
                individual.error_metrics = error_metrics
                
            except Exception as e:
                print(f"    Error evaluating individual: {e}")
                individual.fitness = 0.0
                individual.error_metrics = {"error": str(e)}
    
    def _create_reduced_mechanism(self, individual: Individual) -> Mechanism:
        """Create reduced mechanism from individual."""
        
        # Get species to remove
        species_to_remove = []
        for i, keep in enumerate(individual.species_mask):
            if not keep:
                species_to_remove.append(self.mechanism.species_names[i])
        
        # Create copy of mechanism
        reduced_mechanism = copy.deepcopy(self.mechanism)
        
        # Remove species
        reduced_mechanism.remove_species(species_to_remove)
        
        return reduced_mechanism
    
    def _evaluate_fitness(self, reduced_mechanism: Mechanism) -> Tuple[float, Dict[str, float]]:
        """Evaluate fitness of reduced mechanism.
        
        Args:
            reduced_mechanism: Reduced mechanism to evaluate
            
        Returns:
            Tuple of (fitness_score, error_metrics)
        """
        
        # Save reduced mechanism temporarily
        temp_file = "temp_reduced_mechanism.yaml"
        reduced_mechanism.save(temp_file)
        
        try:
            # Validate against test cases
            validation_results = self.validator.validate_all_cases(
                mechanism_file=temp_file,
                cases=self.validation_cases,
                output_dir="temp_validation"
            )
            
            # Calculate fitness
            fitness = self._calculate_fitness(validation_results, reduced_mechanism)
            
            # Calculate error metrics
            error_metrics = self._calculate_error_metrics(validation_results)
            
            return fitness, error_metrics
            
        finally:
            # Clean up temporary files
            Path(temp_file).unlink(missing_ok=True)
            import shutil
            shutil.rmtree("temp_validation", ignore_errors=True)
    
    def _calculate_fitness(self, validation_results: Dict, reduced_mechanism: Mechanism) -> float:
        """Calculate fitness score for reduced mechanism."""
        
        # Base fitness from validation results
        if validation_results["overall_passed"]:
            base_fitness = 100.0
        else:
            base_fitness = 50.0
        
        # Penalty for number of species
        n_species = len(reduced_mechanism.species_names)
        species_penalty = max(0, n_species - self.reduction_target.max_species) * 2.0
        
        # Bonus for meeting target
        if n_species <= self.reduction_target.max_species:
            species_bonus = 20.0
        else:
            species_bonus = 0.0
        
        # Calculate total fitness
        fitness = base_fitness - species_penalty + species_bonus
        
        return max(0.0, fitness)
    
    def _calculate_error_metrics(self, validation_results: Dict) -> Dict[str, float]:
        """Calculate error metrics from validation results."""
        
        error_metrics = {}
        
        for case, results in validation_results["case_results"].items():
            if "error" not in results:
                error_metrics[f"{case}_ch4_conversion"] = results.get("ch4_conversion", 0.0)
                error_metrics[f"{case}_h2_co_ratio"] = results.get("h2_co_ratio", 0.0)
                error_metrics[f"{case}_passed"] = 1.0 if results.get("passed", False) else 0.0
            else:
                error_metrics[f"{case}_error"] = 1.0
        
        return error_metrics
    
    def _check_convergence(self) -> bool:
        """Check if algorithm has converged."""
        
        if len(self.convergence_data) < 10:
            return False
        
        # Check if fitness has plateaued
        recent_fitness = [data["best_fitness"] for data in self.convergence_data[-10:]]
        fitness_std = np.std(recent_fitness)
        
        return fitness_std < 0.1  # Converged if fitness variation is small
    
    def _create_next_generation(self, population: List[Individual]) -> List[Individual]:
        """Create next generation using genetic operators."""
        
        new_population = []
        
        # Elitism: keep best individuals
        sorted_population = sorted(population, key=lambda x: x.fitness, reverse=True)
        for i in range(self.ga_params.elite_size):
            new_population.append(copy.deepcopy(sorted_population[i]))
        
        # Generate offspring
        while len(new_population) < self.ga_params.population_size:
            # Selection
            parent1 = self._tournament_selection(population)
            parent2 = self._tournament_selection(population)
            
            # Crossover
            if random.random() < self.ga_params.crossover_rate:
                child1, child2 = self._crossover(parent1, parent2)
            else:
                child1, child2 = copy.deepcopy(parent1), copy.deepcopy(parent2)
            
            # Mutation
            if random.random() < self.ga_params.mutation_rate:
                child1 = self._mutate(child1)
            if random.random() < self.ga_params.mutation_rate:
                child2 = self._mutate(child2)
            
            new_population.extend([child1, child2])
        
        return new_population[:self.ga_params.population_size]
    
    def _tournament_selection(self, population: List[Individual]) -> Individual:
        """Tournament selection."""
        
        tournament = random.sample(population, self.ga_params.tournament_size)
        return max(tournament, key=lambda x: x.fitness)
    
    def _crossover(self, parent1: Individual, parent2: Individual) -> Tuple[Individual, Individual]:
        """Single-point crossover."""
        
        # Choose crossover point
        crossover_point = random.randint(1, self.n_species - 1)
        
        # Create children
        child1_mask = np.concatenate([
            parent1.species_mask[:crossover_point],
            parent2.species_mask[crossover_point:]
        ])
        child2_mask = np.concatenate([
            parent2.species_mask[:crossover_point],
            parent1.species_mask[crossover_point:]
        ])
        
        child1 = Individual(species_mask=child1_mask)
        child2 = Individual(species_mask=child2_mask)
        
        return child1, child2
    
    def _mutate(self, individual: Individual) -> Individual:
        """Mutation operator."""
        
        mutated_mask = individual.species_mask.copy()
        
        # Randomly flip some bits
        for i in range(self.n_species):
            if random.random() < 0.1:  # 10% chance per bit
                mutated_mask[i] = not mutated_mask[i]
        
        # Ensure minimum number of species
        if np.sum(mutated_mask) < 5:
            # Add some random species
            available_indices = np.where(~mutated_mask)[0]
            n_to_add = min(5, len(available_indices))
            selected_indices = random.sample(available_indices.tolist(), n_to_add)
            mutated_mask[selected_indices] = True
        
        return Individual(species_mask=mutated_mask)
    
    def _validate_reduced_mechanism(self, reduced_mechanism: Mechanism) -> Dict:
        """Validate reduced mechanism against all test cases."""
        
        # Save reduced mechanism
        reduced_file = "reduced_mechanism.yaml"
        reduced_mechanism.save(reduced_file)
        
        try:
            # Validate against all cases
            validation_results = self.validator.validate_all_cases(
                mechanism_file=reduced_file,
                cases=self.validation_cases,
                output_dir="reduced_validation"
            )
            
            return validation_results
            
        finally:
            # Clean up
            Path(reduced_file).unlink(missing_ok=True)
    
    def _save_results(
        self,
        best_individual: Individual,
        reduced_mechanism: Mechanism,
        validation_results: Dict,
        output_dir: str
    ):
        """Save reduction results."""
        
        output_path = Path(output_dir) / "mechanism_reduction"
        output_path.mkdir(exist_ok=True)
        
        # Save reduced mechanism
        reduced_mechanism.save(str(output_path / "reduced_mechanism.yaml"))
        
        # Save convergence data
        convergence_df = pd.DataFrame(self.convergence_data)
        convergence_df.to_csv(output_path / "convergence.csv", index=False)
        
        # Save best individual
        best_individual_data = {
            "fitness": best_individual.fitness,
            "n_species": np.sum(best_individual.species_mask),
            "species_mask": best_individual.species_mask.tolist(),
            "error_metrics": best_individual.error_metrics
        }
        
        with open(output_path / "best_individual.yaml", 'w') as f:
            yaml.dump(best_individual_data, f, default_flow_style=False)
        
        # Save validation results
        with open(output_path / "validation_results.yaml", 'w') as f:
            yaml.dump(validation_results, f, default_flow_style=False)
        
        # Generate plots
        self._generate_reduction_plots(output_path)
    
    def _generate_reduction_plots(self, output_path: Path):
        """Generate mechanism reduction plots."""
        
        # Convergence plot
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle("GA-GNN Mechanism Reduction Results", fontsize=16)
        
        # Fitness convergence
        ax1 = axes[0, 0]
        generations = [data["generation"] for data in self.convergence_data]
        best_fitness = [data["best_fitness"] for data in self.convergence_data]
        avg_fitness = [data["avg_fitness"] for data in self.convergence_data]
        
        ax1.plot(generations, best_fitness, 'b-', linewidth=2, label='Best')
        ax1.plot(generations, avg_fitness, 'r--', linewidth=2, label='Average')
        ax1.set_xlabel("Generation")
        ax1.set_ylabel("Fitness")
        ax1.set_title("Fitness Convergence")
        ax1.legend()
        ax1.grid(True)
        
        # Species count convergence
        ax2 = axes[0, 1]
        n_species = [data["n_species"] for data in self.convergence_data]
        ax2.plot(generations, n_species, 'g-', linewidth=2)
        ax2.axhline(y=self.reduction_target.max_species, color='r', linestyle='--', 
                   label=f'Target: {self.reduction_target.max_species}')
        ax2.set_xlabel("Generation")
        ax2.set_ylabel("Number of Species")
        ax2.set_title("Species Count Evolution")
        ax2.legend()
        ax2.grid(True)
        
        # Species selection heatmap
        ax3 = axes[1, 0]
        species_mask = self.best_individuals[-1].species_mask
        species_names = self.mechanism.species_names
        
        # Create heatmap data
        heatmap_data = species_mask.reshape(1, -1)
        im = ax3.imshow(heatmap_data, cmap='RdYlGn', aspect='auto')
        ax3.set_xlabel("Species Index")
        ax3.set_title("Species Selection (Green=Selected)")
        ax3.set_yticks([])
        
        # Add species names as x-tick labels (every 10th)
        tick_positions = range(0, len(species_names), 10)
        ax3.set_xticks(tick_positions)
        ax3.set_xticklabels([species_names[i] for i in tick_positions], rotation=45)
        
        # Fitness distribution
        ax4 = axes[1, 1]
        final_fitness = [ind.fitness for ind in self.best_individuals[-10:]]
        ax4.hist(final_fitness, bins=10, alpha=0.7, edgecolor='black')
        ax4.set_xlabel("Fitness")
        ax4.set_ylabel("Frequency")
        ax4.set_title("Final Generation Fitness Distribution")
        ax4.grid(True)
        
        plt.tight_layout()
        plt.savefig(output_path / "reduction_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Species importance plot
        self._plot_species_importance(output_path)
    
    def _plot_species_importance(self, output_path: Path):
        """Plot species importance based on selection frequency."""
        
        # Calculate selection frequency across all generations
        selection_frequency = np.zeros(self.n_species)
        
        for individual in self.best_individuals:
            selection_frequency += individual.species_mask
        
        selection_frequency /= len(self.best_individuals)
        
        # Sort by frequency
        sorted_indices = np.argsort(selection_frequency)[::-1]
        sorted_frequencies = selection_frequency[sorted_indices]
        sorted_names = [self.mechanism.species_names[i] for i in sorted_indices]
        
        # Plot top 20 species
        top_n = min(20, len(sorted_names))
        fig, ax = plt.subplots(figsize=(12, 8))
        
        bars = ax.bar(range(top_n), sorted_frequencies[:top_n])
        ax.set_xlabel("Species")
        ax.set_ylabel("Selection Frequency")
        ax.set_title("Species Importance (Selection Frequency)")
        ax.set_xticks(range(top_n))
        ax.set_xticklabels(sorted_names[:top_n], rotation=45, ha='right')
        
        # Color bars by frequency
        for i, bar in enumerate(bars):
            bar.set_color(plt.cm.viridis(sorted_frequencies[i]))
        
        plt.tight_layout()
        plt.savefig(output_path / "species_importance.png", dpi=300, bbox_inches='tight')
        plt.close()
