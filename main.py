#!/usr/bin/env python3
"""HP-POX Benchmark Framework - Main CLI Interface.

Industry-grade 1-D reactor framework for HP-POX validation and plant applications.
Supports mechanism reduction, plasma surrogates, and comprehensive validation.
"""

import argparse
import sys
import yaml
from pathlib import Path
from typing import List, Optional

from validation.hp_pox_validator import HP_POXValidator
from applications.plant_application import PlantApplication, OperatingEnvelope, PlantValidationData
from surrogates.plasma_models import PlasmaSurrogateManager, PlasmaConfig
from reduction.ga_gnn_pipeline import GAGNNPipeline, ReductionTarget, GAParameters, GNNParameters


def main():
    """Main CLI entry point."""
    
    parser = argparse.ArgumentParser(
        description="HP-POX Benchmark Framework - Industry-grade 1-D reactor modeling",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate HP-POX Cases 1 & 4
  python main.py validate --case case_1,case_4 --mechanism data/gri30.yaml

  # Run plant application sweep
  python main.py plant --geometry plant_a --mechanism data/gri30.yaml --sweep

  # Run mechanism reduction
  python main.py reduce --mechanism data/gri30.yaml --target-species 30 --cases case_1,case_4

  # Test plasma surrogates
  python main.py plasma --mechanism data/gri30.yaml --type thermal --power-range 50,100,150,200

  # Full validation pipeline
  python main.py pipeline --mechanism data/gri30.yaml --output results/full_validation
        """
    )
    
    # Global arguments
    parser.add_argument("--mechanism", required=True, help="Path to mechanism file (.yaml/.cti)")
    parser.add_argument("--output", default="results", help="Output directory")
    parser.add_argument("--config", default="config", help="Configuration directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    # Subcommands
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Validation command
    validate_parser = subparsers.add_parser("validate", help="Validate against HP-POX benchmark")
    validate_parser.add_argument("--case", default="case_1,case_4", 
                               help="Comma-separated list of cases to validate")
    validate_parser.add_argument("--n-points", type=int, default=1000, 
                               help="Number of axial discretization points")
    
    # Plant application command
    plant_parser = subparsers.add_parser("plant", help="Plant application analysis")
    plant_parser.add_argument("--geometry", choices=["plant_a", "plant_b"], required=True,
                            help="Plant geometry to use")
    plant_parser.add_argument("--sweep", action="store_true", 
                            help="Run operating envelope sweep")
    plant_parser.add_argument("--validate", help="Path to plant validation data file")
    plant_parser.add_argument("--n-points", type=int, default=1000,
                            help="Number of axial discretization points")
    
    # Mechanism reduction command
    reduce_parser = subparsers.add_parser("reduce", help="Mechanism reduction using GA-GNN")
    reduce_parser.add_argument("--target-species", type=int, required=True,
                             help="Target number of species")
    reduce_parser.add_argument("--cases", default="case_1,case_4",
                             help="Validation cases to use")
    reduce_parser.add_argument("--max-error", type=float, default=5.0,
                             help="Maximum allowed error percentage")
    reduce_parser.add_argument("--generations", type=int, default=100,
                             help="Number of GA generations")
    reduce_parser.add_argument("--population-size", type=int, default=50,
                             help="GA population size")
    
    # Plasma surrogate command
    plasma_parser = subparsers.add_parser("plasma", help="Plasma surrogate analysis")
    plasma_parser.add_argument("--type", choices=["thermal", "radical"], required=True,
                             help="Plasma surrogate type")
    plasma_parser.add_argument("--power-range", help="Comma-separated power range (kW)")
    plasma_parser.add_argument("--injection-location", type=float, default=0.5,
                             help="Injection location (m from inlet)")
    plasma_parser.add_argument("--case", default="case_1", help="Base case for plasma analysis")
    
    # Full pipeline command
    pipeline_parser = subparsers.add_parser("pipeline", help="Run complete validation pipeline")
    pipeline_parser.add_argument("--skip-reduction", action="store_true",
                               help="Skip mechanism reduction step")
    pipeline_parser.add_argument("--skip-plasma", action="store_true",
                               help="Skip plasma analysis step")
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Run command
    try:
        if args.command == "validate":
            run_validation(args)
        elif args.command == "plant":
            run_plant_application(args)
        elif args.command == "reduce":
            run_mechanism_reduction(args)
        elif args.command == "plasma":
            run_plasma_analysis(args)
        elif args.command == "pipeline":
            run_full_pipeline(args)
        else:
            print(f"Unknown command: {args.command}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def run_validation(args):
    """Run HP-POX validation."""
    
    print("Running HP-POX validation...")
    
    # Parse cases
    cases = [case.strip() for case in args.case.split(",")]
    
    # Initialize validator
    validator = HP_POXValidator()
    
    # Run validation
    results = validator.validate_all_cases(
        mechanism_file=args.mechanism,
        cases=cases,
        output_dir=args.output
    )
    
    # Print summary
    print("\nValidation Summary:")
    print(f"Mechanism: {args.mechanism}")
    print(f"Cases validated: {len(cases)}")
    print(f"Overall passed: {results['overall_passed']}")
    
    for case, result in results["case_results"].items():
        status = "PASSED" if result.get("passed", False) else "FAILED"
        print(f"  {case}: {status}")
        if "error" in result:
            print(f"    Error: {result['error']}")
    
    print(f"\nResults saved to: {args.output}")


def run_plant_application(args):
    """Run plant application analysis."""
    
    print(f"Running plant application for {args.geometry}...")
    
    # Initialize plant application
    plant_app = PlantApplication()
    
    if args.sweep:
        # Run operating envelope sweep
        envelope = OperatingEnvelope(
            pressure_range_bar_g=[30, 80],
            temperature_range_C=[1000, 1500],
            flow_range_kgph=[500, 2000],
            o2_c_ratio_range=[0.5, 0.7],
            diluent_fractions={"H2O": 0.1, "N2": 0.05}
        )
        
        results = plant_app.run_operating_envelope_sweep(
            geometry_name=args.geometry,
            mechanism_file=args.mechanism,
            envelope=envelope,
            output_dir=args.output
        )
        
        print(f"\nSweep Results:")
        print(f"Total simulations: {results['n_simulations']}")
        print(f"Successful: {results['successful']}")
        
    elif args.validate:
        # Validate against plant data
        plant_data = load_plant_validation_data(args.validate)
        
        results = plant_app.validate_against_plant_data(
            geometry_name=args.geometry,
            mechanism_file=args.mechanism,
            plant_data=plant_data,
            output_dir=args.output
        )
        
        print(f"\nPlant Validation Results:")
        print(f"Cases validated: {results['n_cases']}")
        
        for result in results["results"]:
            status = "PASSED" if result.get("passed", False) else "FAILED"
            print(f"  {result['case']}: {status}")
    
    print(f"\nResults saved to: {args.output}")


def run_mechanism_reduction(args):
    """Run mechanism reduction."""
    
    print("Running GA-GNN mechanism reduction...")
    
    # Parse cases
    cases = [case.strip() for case in args.cases.split(",")]
    
    # Create reduction target
    reduction_target = ReductionTarget(
        max_species=args.target_species,
        max_error_pct=args.max_error,
        validation_cases=cases,
        priority_species=["CH4", "O2", "H2", "CO", "CO2", "H2O", "N2"]
    )
    
    # Create GA parameters
    ga_params = GAParameters(
        generations=args.generations,
        population_size=args.population_size
    )
    
    # Initialize pipeline
    pipeline = GAGNNPipeline(
        mechanism_file=args.mechanism,
        validation_cases=cases,
        reduction_target=reduction_target,
        ga_params=ga_params
    )
    
    # Run reduction
    results = pipeline.run_reduction(output_dir=args.output)
    
    # Print summary
    best_individual = results["best_individual"]
    print(f"\nReduction Results:")
    print(f"Target species: {args.target_species}")
    print(f"Final species: {np.sum(best_individual.species_mask)}")
    print(f"Best fitness: {best_individual.fitness:.2f}")
    print(f"Generations: {len(results['convergence_data'])}")
    
    print(f"\nResults saved to: {args.output}/mechanism_reduction")


def run_plasma_analysis(args):
    """Run plasma surrogate analysis."""
    
    print(f"Running plasma {args.type} analysis...")
    
    # Parse power range
    if args.power_range:
        power_range = [float(p.strip()) for p in args.power_range.split(",")]
    else:
        power_range = [50, 100, 150, 200, 250, 300]
    
    # Initialize plasma manager
    plasma_manager = PlasmaSurrogateManager()
    
    # Create base reactor for analysis
    from applications.plant_application import PlantApplication
    plant_app = PlantApplication()
    reactor = plant_app.create_reactor("plant_a", args.mechanism)
    
    # Set up base case
    from validation.hp_pox_validator import HP_POXValidator
    validator = HP_POXValidator()
    
    # Get base case data
    case_data = validator.config["cases"][args.case]
    inlet = plant_app._create_inlet_conditions(case_data)
    reactor.set_inlet_conditions(inlet)
    
    # Run base simulation
    base_result = reactor.solve()
    
    # Run plasma sweep
    results = plasma_manager.run_plasma_sweep(
        reactor=reactor,
        base_result=base_result,
        plasma_type=args.type,
        power_range=power_range,
        output_dir=args.output
    )
    
    print(f"\nPlasma Analysis Results:")
    print(f"Type: {args.type}")
    print(f"Power range: {power_range}")
    print(f"Tests run: {results['n_tests']}")
    print(f"Successful: {results['successful']}")
    
    print(f"\nResults saved to: {args.output}/plasma_sweeps")


def run_full_pipeline(args):
    """Run complete validation pipeline."""
    
    print("Running complete validation pipeline...")
    
    # Step 1: HP-POX validation
    print("\n1. HP-POX Validation")
    run_validation(args)
    
    # Step 2: Mechanism reduction (if not skipped)
    if not args.skip_reduction:
        print("\n2. Mechanism Reduction")
        reduction_args = argparse.Namespace(
            mechanism=args.mechanism,
            target_species=30,
            cases="case_1,case_4",
            max_error=5.0,
            generations=50,
            population_size=30,
            output=args.output
        )
        run_mechanism_reduction(reduction_args)
    
    # Step 3: Plant application
    print("\n3. Plant Application")
    plant_args = argparse.Namespace(
        geometry="plant_a",
        sweep=True,
        validate=None,
        n_points=1000,
        mechanism=args.mechanism,
        output=args.output,
        verbose=args.verbose
    )
    run_plant_application(plant_args)
    
    # Step 4: Plasma analysis (if not skipped)
    if not args.skip_plasma:
        print("\n4. Plasma Analysis")
        plasma_args = argparse.Namespace(
            type="thermal",
            power_range="50,100,150,200",
            injection_location=0.5,
            case="case_1",
            mechanism=args.mechanism,
            output=args.output,
            verbose=args.verbose
        )
        run_plasma_analysis(plasma_args)
    
    print(f"\nComplete pipeline finished. Results saved to: {args.output}")


def load_plant_validation_data(file_path: str) -> List[PlantValidationData]:
    """Load plant validation data from file."""
    
    with open(file_path, 'r') as f:
        data = yaml.safe_load(f)
    
    plant_data = []
    for case_data in data["plant_data"]:
        plant_data.append(PlantValidationData(
            case_name=case_data["case_name"],
            pressure_bar_g=case_data["pressure_bar_g"],
            temperature_C=case_data["temperature_C"],
            flow_rate_kgph=case_data["flow_rate_kgph"],
            o2_c_ratio=case_data["o2_c_ratio"],
            ch4_conversion_pct=case_data["ch4_conversion_pct"],
            co2_conversion_pct=case_data["co2_conversion_pct"],
            h2_co_ratio=case_data["h2_co_ratio"],
            outlet_composition=case_data["outlet_composition"],
            pressure_drop_bar=case_data["pressure_drop_bar"],
            ignition_length_m=case_data["ignition_length_m"],
            ng_composition=case_data["ng_composition"],
            inlet_temperatures=case_data["inlet_temperatures"]
        ))
    
    return plant_data


if __name__ == "__main__":
    main()
