/*!
 * \file keywords.h
 *
 * Deon Garrett
 * University of Memphis
 * jdgarrett@gmail.com
 *
 * Defines the known keywords in the system.  Provides a singleton
 * instance which is used to access the keywords.
 *
 */

#ifndef _KEYWORDS_H_
#define _KEYWORDS_H_

#include <string>

using namespace std;

/*!
 * \class store the available keywords/values for sls configuration
 *
 * \author deong
 * \date 05/09/2007
 */
class keywords
{
public:
    
    // generic parameters
    static const string TRIALS;
    static const string RANDOM_SEED;

    // problem parameters
    static const string PROBLEM;
    static const string PROBLEM_DATA;
    static const string DIMENSIONS;
    static const string OBJECTIVES;
    static const string TRUE_PARETO_FRONT;
    static const string PENALIZATION_FACTOR;
    static const string UNASSIGNED_PENALTY;
    static const string CAPACITY_PENALTY;
    
    // algorithm parameters
    static const string ALGORITHM;
    static const string GENITOR;
    static const string SIMPLE_GA;
    static const string NSGA2;
    static const string SPEA2;
    static const string EPSILON_MOEA;
    static const string CHC;
    static const string RANDOM_RESTART_LS;
    static const string MULTIOBJECTIVE_RRLS;
    static const string SAMPLE_AND_MUTATE;
    static const string TWO_PHASE_LS;
    static const string TABU_SEARCH;
    static const string SIMULATED_ANNEALING;
    static const string VARIABLE_DEPTH_SEARCH;
    static const string LANDSCAPE;
    
    // hill climbing parameters
    static const string LOCAL_SEARCH;
    static const string DEBUG_LS_GENERATIONS;
    static const string NEXT_DESCENT;
    static const string STEEPEST_DESCENT;
    static const string LOCAL_SEARCH_ITERATIONS;
    static const string SECOND_PHASE_ITERATIONS;
    static const string NONDOMINATED_POINTS;
    static const string DIVERSIFICATION;
    static const string NUM_WEIGHTS;
    static const string PHASE_TWO_ITERATIONS;
    static const string PHASE_TWO_RUNS;
    static const string INITIAL_50_PERCENTILE;
    static const string FINAL_50_PERCENTILE;
    static const string COOLING_STEPS;
    
    // encoding parameters
    static const string ENCODING;
    static const string BOOLEAN_ENCODING;
    static const string BINARY_ENCODING;
    static const string REAL_ENCODING;
    static const string PERMUTATION_ENCODING;
    static const string INTEGER_ENCODING;
    static const string GAP_ENCODING;
    static const string GSAP_ENCODING;
        
    // evoutionary algorithm parameters
    static const string POPULATION_SIZE;
    
    static const string SELECTION_SCHEME;
    static const string TOURNAMENT_SELECTION;
    static const string RANKING_SELECTION;
    static const string RANDOM_SELECTION;
    static const string RANKING_BIAS;
    
    static const string REPLACEMENT_SCHEME;
    static const string GENERATIONAL_REPLACEMENT;
    static const string ELITIST_REPLACEMENT;
    static const string ELITES;
    static const string TRUNCATION_REPLACEMENT;
    static const string REPLACE_WORST_REPLACEMENT;

    static const string CROSSOVER_OPERATOR;
    static const string CROSSOVER_RATE;
    static const string UNIFORM_CROSSOVER;
    static const string ONE_POINT_CROSSOVER;
    static const string TWO_POINT_CROSSOVER;
    static const string HUX_CROSSOVER;
    static const string ORDER_CROSSOVER;
    static const string CYCLE_CROSSOVER;
    static const string SBX_CROSSOVER;
    static const string SBX_ETA;
    
    static const string MUTATION_OPERATOR;
    static const string MUTATION_RATE;
    static const string BITWISE_MUTATION;
    static const string SWAP_MUTATION;
    static const string GAUSSIAN_MUTATION;
    static const string GAUSSIAN_MUTATION_MU;
    static const string GAUSSIAN_MUTATION_SIGMA;
    static const string POLYNOMIAL_MUTATION;
    static const string POLYNOMIAL_ETA;
    static const string SHIFT_MUTATION;
    static const string SSS_MUTATION;

    static const string REPAIR_OPERATOR;
    static const string GAP_REPAIR_ND;
    static const string GAP_REPAIR_SD;
    static const string GSAP_REPAIR;
    
    static const string ARCHIVE_SIZE;
    static const string SPEA2_K;
    
    // local search operators
    static const string NEIGHBORHOOD;
    static const string SWAP_NEIGHBORHOOD;
    static const string HAMMING_NEIGHBORHOOD;
    static const string SHIFT_NEIGHBORHOOD;
    static const string SSS_NEIGHBORHOOD;
    
    // local search strategies
    static const string HC_STRATEGY;
    static const string STRATEGY_ALL;
    static const string STRATEGY_BEST;
    static const string STRATEGY_WORST;
    static const string STRATEGY_NONE;
    static const string STRATEGY_RANDOM;
    static const string HC_RATE;
    
    // tabu search parameters
    static const string MIN_TABU_TENURE;
    static const string MAX_TABU_TENURE;
    static const string TABU_COMPARATOR;
    static const string TABU_COMPARATOR_ANY;
    static const string TABU_COMPARATOR_ALL;
    
    // performance metrics
    static const string METRIC;
    static const string EVALUATION_COUNTER;
    static const string GENERATION_COUNTER;
    static const string BEST_SOLUTION;
	static const string REPORT_ALL_BEST;
    static const string HYPERVOLUME;
    static const string REFERENCE_POINT;
    static const string POPULATION_ENTROPY;
	static const string PRINT_EVERY_GENERATION;
	
    // termination criteria
    static const string TERMINATOR;
    static const string EVALUATION_LIMIT;
    static const string GENERATION_LIMIT;
    static const string NULL_TERMINATOR;
    static const string MAX_EVALUATIONS;
    static const string MAX_GENERATIONS;

    // landscape analysis
    static const string LANDSCAPE_TOOL;
    static const string RANDOM_WALK;
    static const string RANDOM_WALK_BETWEEN_OPTIMA;
    static const string RUGGEDNESS;
    static const string RUGGEDNESS_BETWEEN_OPTIMA;
    static const string FITNESS_DISTANCE_CORRELATION;
    static const string PERTURBATION_SEARCH;
    static const string INFEASIBILITY_REGION;
    static const string NUM_LOCAL_OPTIMA;
    static const string RANDOM_WALK_LENGTH;
    static const string NUM_RANDOM_WALKS;
    static const string PARETO_PLATEAUS;
    static const string GRAPHVIZ_FILE;
    
    // comparison methods for chromosomes
    static const string COMPARATOR;
    static const string FITNESS_COMPARATOR;
    static const string PARETO_DOMINANCE_COMPARATOR;
    static const string WEAK_DOMINANCE_COMPARATOR;
    static const string STRONG_DOMINANCE_COMPARATOR;
    static const string EPSILON_DOMINANCE_COMPARATOR;
    static const string SCALARIZING_COMPARATOR;
    static const string WEIGHT_VECTOR;
    static const string EPSILON;
};

#endif

