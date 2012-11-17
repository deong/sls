/*!
 * \file keywords.cpp
 *
 * Defines the known keywords in the system.  Provides a singleton
 * instance which is used to access the keywords.
 *
 * Deon Garrett
 * University of Memphis
 * deong@acm.org
 *
 */

#include "keywords.h"

const string keywords::TRIALS = "trials";
const string keywords::RANDOM_SEED = "random_seed";

const string keywords::PROBLEM = "problem";
const string keywords::PROBLEM_DATA = "problem_data";
const string keywords::DIMENSIONS = "dimensions";
const string keywords::OBJECTIVES = "objectives";
const string keywords::TRUE_PARETO_FRONT = "true_pareto_front";
const string keywords::PENALIZATION_FACTOR = "penalization_factor";
const string keywords::UNASSIGNED_PENALTY = "unassigned_penalty";
const string keywords::CAPACITY_PENALTY = "capacity_penalty";

const string keywords::ALGORITHM = "algorithm";
const string keywords::GENITOR = "genitor";
const string keywords::SIMPLE_GA = "simple_ga";
const string keywords::NSGA2 = "nsga2";
const string keywords::SPEA2 = "spea2";
const string keywords::EPSILON_MOEA = "emoea";
const string keywords::CHC = "chc";
const string keywords::RANDOM_RESTART_LS = "rrls";
const string keywords::MULTIOBJECTIVE_RRLS = "morrls";
const string keywords::SAMPLE_AND_MUTATE = "s_and_m";
const string keywords::TWO_PHASE_LS = "two_phase_ls";
const string keywords::TABU_SEARCH = "tabu_search";
const string keywords::SIMULATED_ANNEALING = "simulated_annealing";
const string keywords::VARIABLE_DEPTH_SEARCH = "variable_depth_search";
const string keywords::LANDSCAPE = "landscape";

const string keywords::LOCAL_SEARCH = "local_search";
const string keywords::DEBUG_LS_GENERATIONS = "debug_ls_generations";
const string keywords::NEXT_DESCENT = "next_descent";
const string keywords::STEEPEST_DESCENT = "steepest_descent";
const string keywords::LOCAL_SEARCH_ITERATIONS = "local_search_iterations";
const string keywords::SECOND_PHASE_ITERATIONS = "second_phase_iterations";
const string keywords::NONDOMINATED_POINTS = "nondominated_points";
const string keywords::DIVERSIFICATION = "diversification";
const string keywords::NUM_WEIGHTS = "num_weights";
const string keywords::PHASE_TWO_ITERATIONS = "phase_two_iterations";
const string keywords::PHASE_TWO_RUNS = "phase_two_runs";
const string keywords::INITIAL_50_PERCENTILE = "initial_50_percentile";
const string keywords::FINAL_50_PERCENTILE = "final_50_percentile";
const string keywords::COOLING_STEPS = "cooling_steps";

const string keywords::ENCODING = "encoding";
const string keywords::PERMUTATION_ENCODING = "permutation";
const string keywords::BOOLEAN_ENCODING = "boolean";
const string keywords::BINARY_ENCODING = "binary";
const string keywords::REAL_ENCODING = "real";
const string keywords::INTEGER_ENCODING = "integer";
const string keywords::GAP_ENCODING = "gap";
const string keywords::GSAP_ENCODING = "gsap";

const string keywords::POPULATION_SIZE = "population_size";

const string keywords::SELECTION_SCHEME = "selection";
const string keywords::TOURNAMENT_SELECTION = "tournament";
const string keywords::RANKING_SELECTION = "ranking";
const string keywords::RANDOM_SELECTION = "random";
const string keywords::RANKING_BIAS = "ranking_bias";

const string keywords::REPLACEMENT_SCHEME = "replacement";
const string keywords::GENERATIONAL_REPLACEMENT = "generational";
const string keywords::ELITIST_REPLACEMENT = "elitist";
const string keywords::ELITES = "elites";
const string keywords::TRUNCATION_REPLACEMENT = "truncation";
const string keywords::REPLACE_WORST_REPLACEMENT = "replace_worst";

const string keywords::CROSSOVER_OPERATOR = "crossover";
const string keywords::CROSSOVER_RATE = "crossover_rate";
const string keywords::UNIFORM_CROSSOVER = "uniform";
const string keywords::ONE_POINT_CROSSOVER = "one_point";
const string keywords::TWO_POINT_CROSSOVER = "two_point";
const string keywords::HUX_CROSSOVER = "hux";
const string keywords::ORDER_CROSSOVER = "order";
const string keywords::CYCLE_CROSSOVER = "cycle";
const string keywords::SBX_CROSSOVER = "sbx";
const string keywords::SBX_ETA = "sbx_eta";

const string keywords::MUTATION_OPERATOR = "mutation";
const string keywords::MUTATION_RATE = "mutation_rate";
const string keywords::BITWISE_MUTATION = "bitwise";
const string keywords::SWAP_MUTATION = "swap";
const string keywords::GAUSSIAN_MUTATION = "gaussian";
const string keywords::GAUSSIAN_MUTATION_MU = "mu";
const string keywords::GAUSSIAN_MUTATION_SIGMA = "sigma";
const string keywords::POLYNOMIAL_MUTATION = "polynomial";
const string keywords::POLYNOMIAL_ETA = "polynomial_eta";
const string keywords::SHIFT_MUTATION = "shift";
const string keywords::SSS_MUTATION = "sss";

const string keywords::REPAIR_OPERATOR="repair_operator";
const string keywords::GAP_REPAIR_ND="gap_repair_nd";
const string keywords::GAP_REPAIR_SD="gap_repair_sd";
const string keywords::GSAP_REPAIR="gsap_repair";

const string keywords::ARCHIVE_SIZE = "archive_size";
const string keywords::SPEA2_K = "spea2_k";

const string keywords::NEIGHBORHOOD = "neighborhood";
const string keywords::SWAP_NEIGHBORHOOD = "swap";
const string keywords::HAMMING_NEIGHBORHOOD = "hamming";
const string keywords::SHIFT_NEIGHBORHOOD = "shift";
const string keywords::SSS_NEIGHBORHOOD = "sss";

const string keywords::HC_STRATEGY = "strategy";
const string keywords::STRATEGY_ALL = "all";
const string keywords::STRATEGY_BEST = "best";
const string keywords::STRATEGY_WORST = "worst";
const string keywords::STRATEGY_NONE = "none";
const string keywords::STRATEGY_RANDOM = "random";
const string keywords::HC_RATE = "hc_rate";

const string keywords::MIN_TABU_TENURE = "min_tabu_tenure";
const string keywords::MAX_TABU_TENURE = "max_tabu_tenure";
const string keywords::TABU_COMPARATOR = "tabu_comparator";
const string keywords::TABU_COMPARATOR_ANY = "any";
const string keywords::TABU_COMPARATOR_ALL = "all";

const string keywords::METRIC = "metric";
const string keywords::EVALUATION_COUNTER = "evaluation_count";
const string keywords::GENERATION_COUNTER = "generation_count";
const string keywords::BEST_SOLUTION = "best_solution";
const string keywords::REPORT_ALL_BEST = "report_all_best";
const string keywords::HYPERVOLUME = "hypervolume";
const string keywords::REFERENCE_POINT = "reference_point";
const string keywords::POPULATION_ENTROPY = "population_entropy";
const string keywords::PRINT_EVERY_GENERATION = "print_every_generation";

const string keywords::TERMINATOR = "terminator";
const string keywords::EVALUATION_LIMIT = "evaluation_limit";
const string keywords::GENERATION_LIMIT = "generation_limit";
const string keywords::NULL_TERMINATOR = "null";
const string keywords::MAX_EVALUATIONS = "max_evaluations";
const string keywords::MAX_GENERATIONS = "max_generations";

const string keywords::LANDSCAPE_TOOL = "landscape_tool";
const string keywords::PERTURBATION_SEARCH = "perturbation_search";
const string keywords::FITNESS_DISTANCE_CORRELATION = "fitness_distance_correlation";
const string keywords::RANDOM_WALK = "random_walk";
const string keywords::RANDOM_WALK_BETWEEN_OPTIMA = "random_walk_between_optima";
const string keywords::RUGGEDNESS = "ruggedness";
const string keywords::RUGGEDNESS_BETWEEN_OPTIMA = "ruggedness_between_optima";
const string keywords::INFEASIBILITY_REGION = "infeasibility_region";
const string keywords::NUM_LOCAL_OPTIMA = "num_local_optima";
const string keywords::RANDOM_WALK_LENGTH = "random_walk_length";
const string keywords::NUM_RANDOM_WALKS = "num_random_walks";
const string keywords::PARETO_PLATEAUS = "pareto_plateaus";
const string keywords::GRAPHVIZ_FILE = "graphviz_file";

const string keywords::COMPARATOR = "comparator";
const string keywords::FITNESS_COMPARATOR = "fitness";
const string keywords::PARETO_DOMINANCE_COMPARATOR = "pareto_dominance";
const string keywords::WEAK_DOMINANCE_COMPARATOR = "weak_dominance";
const string keywords::STRONG_DOMINANCE_COMPARATOR = "strong_dominance";
const string keywords::EPSILON_DOMINANCE_COMPARATOR = "epsilon_dominance";
const string keywords::SCALARIZING_COMPARATOR = "scalarizing";
const string keywords::WEIGHT_VECTOR = "weight_vector";
const string keywords::EPSILON = "epsilon";
