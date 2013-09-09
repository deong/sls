/*!
 * \file slsmain.cpp
 *
 * driver program for the stochastic local search framework
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include <cstdlib>
#include <string>
#include "sls.h"
#include "genitor.h"
#include "simplega.h"
#include "chc.h"
#include "rrls.h"
#include "morrls.h"
#include "twophasels.h"
#include "nsga2.h"
#include "spea2.h"
#include "emoea.h"
#include "landscape.h"
#include "sandm.h"
#include "chromosome.h"
#include "encoding.h"
#include "comparator.h"
#include "problems.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

/*!
 * \mainpage
 *
 * 
 * This code formed the experimental basis for my dissertation work on
 * multiobjective landscape analysis. It is a fairly general platform
 * for implementing multiobjective metaheuristics for numeric or
 * combinatorial optimization. The focus of the system is on making it
 * fairly straightforward to describe an algorithm, and this comes
 * somewhat at the expense of the clarity of the internal code. I'm
 * not too happy with the system as it stands today, as I think the
 * complexity has outpaced the benefits. I'm currently working on a
 * ground-up rewrite of the package that will remove much of the
 * incidental complexity, but as of now, the newer package is not yet
 * ready to use (and as of mid-2011, I'm no longer devoting much time
 * to this line of research, so the future of that code is unknown).
 * 
 * The code is all C++, and relies rather heavily on features of
 * modern C++. Most of the code gets pulled into a single giant
 * template instantiation, which has a number of bad consequences. The
 * most obvious is that compile times are quite long, and very little
 * can be done in the way of separate or incremental compilation.
 */

/*!
 * \brief main driver for sls application
 *
 * \author deong
 * \date 05/12/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
int main(int argc, char** argv)
{
    if(argc < 2)
    {
        cerr << "usage: sls [config_file]+" << endl;
        exit(1);
    }

    // read any configuration files from the command line
    for(int i=1; i<argc; i++)
    {
        configuration::read_configuration_file(argv[i]);
    }

    // figure out how many trials to run
    int trials = 1;
    configuration::integer_parameter(keywords::TRIALS, trials);

    // what type of algorithm are we using
    string alg;
    configuration::string_parameter(keywords::ALGORITHM, alg, true);

    // what type of encoding are we using
    string enc;
    configuration::string_parameter(keywords::ENCODING, enc, true);

    // run the trials
    for(int i=1; i<=trials; i++)
    {
        cout << "trial " << i << endl;

        // initialize the random number generator
        mtrandom mt;
        cout << "seed: " << mtrandom::seed() << endl;
    

        if(alg == keywords::SIMPLE_GA)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                simple_ga<chromosome,binary_encoding> sga;
                sga.initialize();
                sga.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                simple_ga<chromosome,boolean_encoding> sga;
                sga.initialize();
                sga.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                simple_ga<chromosome,real_encoding> sga;
                sga.initialize();
                sga.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                simple_ga<chromosome,permutation_encoding> sga;
                sga.initialize();
                sga.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                simple_ga<chromosome,integer_encoding> sga;
                sga.initialize();
                sga.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                simple_ga<chromosome,gap_encoding> sga;
                sga.initialize();
                sga.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                simple_ga<chromosome,gsap_encoding> sga;
                sga.initialize();
                sga.run();
            }
        }
        else if(alg == keywords::GENITOR)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                genitor<chromosome,binary_encoding> gen;
                gen.initialize();
                gen.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                genitor<chromosome,boolean_encoding> gen;
                gen.initialize();
                gen.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                genitor<chromosome,real_encoding> gen;
                gen.initialize();
                gen.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                genitor<chromosome,permutation_encoding> gen;
                gen.initialize();
                gen.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                genitor<chromosome,integer_encoding> gen;
                gen.initialize();
                gen.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                genitor<chromosome,gap_encoding> gen;
                gen.initialize();
                gen.run();
            }
            else if(enc == keywords::GSAP_ENCODING)
            {
                genitor<chromosome,gsap_encoding> gen;
                gen.initialize();
                gen.run();
            }
        }
        else if(alg == keywords::CHC)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                chc<chromosome,binary_encoding> chc;
                chc.initialize();
                chc.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                chc<chromosome,boolean_encoding> chc;
                chc.initialize();
                chc.run();
            }
        }
        else if(alg == keywords::RANDOM_RESTART_LS)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                rrls<chromosome,binary_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                rrls<chromosome,boolean_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                rrls<chromosome,real_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                rrls<chromosome,permutation_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                rrls<chromosome,integer_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                rrls<chromosome,gap_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                rrls<chromosome,gsap_encoding> ls;
                ls.initialize();
                ls.run();
            }
        }
        else if(alg == keywords::MULTIOBJECTIVE_RRLS)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                morrls<chromosome,binary_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                morrls<chromosome,boolean_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                morrls<chromosome,real_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                morrls<chromosome,permutation_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                morrls<chromosome,integer_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                morrls<chromosome,gap_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                morrls<chromosome,gsap_encoding> ls;
                ls.initialize();
                ls.run();
            }
        }
        else if(alg == keywords::TWO_PHASE_LS)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                twophasels<chromosome,binary_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                twophasels<chromosome,boolean_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                twophasels<chromosome,real_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                twophasels<chromosome,permutation_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                twophasels<chromosome,integer_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                twophasels<chromosome,gap_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else 
                                if(enc == keywords::GSAP_ENCODING)
            {
                twophasels<chromosome,gsap_encoding> ls;
                ls.initialize();
                ls.run();
            }
        }
        else if(alg == keywords::NSGA2)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                nsga2<binary_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                nsga2<boolean_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                nsga2<real_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                nsga2<permutation_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                nsga2<integer_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                nsga2<gap_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                nsga2<gsap_encoding> moea;
                moea.initialize();
                moea.run();
            }
        }
        else if(alg == keywords::SPEA2)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                spea2<binary_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                spea2<boolean_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                spea2<real_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                spea2<permutation_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                spea2<integer_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                spea2<gap_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                spea2<gsap_encoding> moea;
                moea.initialize();
                moea.run();
            }
        }
        else if(alg == keywords::EPSILON_MOEA)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                emoea<binary_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                emoea<boolean_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                emoea<real_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                emoea<permutation_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                emoea<integer_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                emoea<gap_encoding> moea;
                moea.initialize();
                moea.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                emoea<gsap_encoding> moea;
                moea.initialize();
                moea.run();
            }
        }
        else if(alg == keywords::LANDSCAPE)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                landscape<chromosome,binary_encoding>* lndscp =
                    landscape_factory<chromosome,binary_encoding>::construct();
                lndscp->run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                landscape<chromosome,boolean_encoding>* lndscp =
                    landscape_factory<chromosome,boolean_encoding>::construct();
                lndscp->run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                landscape<chromosome,real_encoding>* lndscp =
                    landscape_factory<chromosome,real_encoding>::construct();
                lndscp->run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                landscape<chromosome,permutation_encoding>* lndscp =
                    landscape_factory<chromosome,permutation_encoding>::construct();
                lndscp->run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                landscape<chromosome,integer_encoding>* lndscp =
                    landscape_factory<chromosome,integer_encoding>::construct();
                lndscp->run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                landscape<chromosome,gap_encoding>* lndscp =
                    landscape_factory<chromosome,gap_encoding>::construct();
                lndscp->run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                landscape<chromosome,gsap_encoding>* lndscp =
                    landscape_factory<chromosome,gsap_encoding>::construct();
                lndscp->run();
            }
        }
        else if(alg == keywords::SAMPLE_AND_MUTATE)
        {
            if(enc == keywords::BINARY_ENCODING)
            {
                sandm<chromosome,binary_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::BOOLEAN_ENCODING)
            {
                sandm<chromosome,boolean_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::REAL_ENCODING)
            {
                sandm<chromosome,real_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::PERMUTATION_ENCODING)
            {
                sandm<chromosome,permutation_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::INTEGER_ENCODING)
            {
                sandm<chromosome,integer_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else if(enc == keywords::GAP_ENCODING)
            {
                sandm<chromosome,gap_encoding> ls;
                ls.initialize();
                ls.run();
            }
            else
                                if(enc == keywords::GSAP_ENCODING)
            {
                sandm<chromosome,gsap_encoding> ls;
                ls.initialize();
                ls.run();
            }
        }
        else
        {
            cerr << "illegal algorithm: " << alg << " specified" << endl;
            exit(1);
        }
        
        cout << endl;
    }
    
    return 0;
}

