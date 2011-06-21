/*!
 * \file landscape.cpp
 *
 * Deon Garrett
 * deong@acm.org
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <cmath>
#include <ios>
#include <iomanip>
#include "landscape.h"
#include "sls.h"
#include "pfront.h"
#include "comparator.h"
#include "localsearch.h"
#include "utilities.h"

using namespace std;

/*!
 * \brief constructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
landscape<Chromosome,Encoding>::landscape() :
    m_hc(0)
{
}

/*!
 * \brief destructor
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
landscape<Chromosome,Encoding>::~landscape()
{
    delete m_hc;
}

/*!
 * \brief initialize the landscape components
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void landscape<Chromosome,Encoding>::initialize()
{
    sls<Chromosome,Encoding>::initialize();

    // construct the hill climber
    local_search_factory<Chromosome,Encoding> lsf;
    lsf.set_prefix("ls_");
    m_hc = lsf.construct();

    // how many tabu search iterations
    configuration::unsigned_integer_parameter(keywords::LOCAL_SEARCH_ITERATIONS, m_lsiter, true);
}

/*!
 * \brief search neighborhood of nondominated solutions
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void perturbation_search<Chromosome,Encoding>::run()
{
    mtrandom mt;

    // read the maximum number of nondominated neighbors to try
    unsigned int num_nondom = 100;
    configuration::unsigned_integer_parameter(keywords::NONDOMINATED_POINTS, num_nondom, false);

    pareto_dominance_comparator<Chromosome,Encoding> pdc;
    scalarizing_comparator<Chromosome,Encoding> comp;

    unsigned int diff = 0;
    unsigned int atts = 0;
    unsigned int iter = 0;
    while(iter++ < num_nondom)
    {
        unsigned int mu = 50;
        unsigned int index = mt.random(0,mu);
    
        comp.randomize_weights(this->m_fitfunc->objectives());
        comp.weights[0] = double(mu-index-1)/(mu-1);
        comp.weights[1] = 1.0 - comp.weights[0];
    
        Chromosome<Encoding> chr(this->m_fitfunc);
        chr.randomize();
        chr.evaluate(this->m_fitfunc);
        this->m_hc->improve(chr, &comp, this->m_fitfunc);

        // at this point, chr is on the front (we assume)
        // time to diversify about its location
        cout << "point located on front: " << endl;
        cout << chr << endl << endl;
        Chromosome<Encoding> popt = chr;

        chr = popt;
        unsigned int indexm1 = index-1;
        if(indexm1 > 0)
        {
            comp.weights[0] = max(0.0,double(mu-indexm1-1)/(mu-1));
            comp.weights[1] = 1.0 - comp.weights[0];
            this->m_hc->improve(chr, &comp, this->m_fitfunc);
            if(chr!=popt)
            {
                diff++;
            }
            atts++;
                
            cout << "point located on front: " << endl;
            cout << chr << endl << endl;
        }

        chr = popt;
        unsigned int indexp1 = index+1;
        if(indexp1 < mu)
        {
            comp.weights[0] = min(1.0,double(mu-indexp1-1)/(mu-1));
            comp.weights[1] = 1.0 - comp.weights[0];
            this->m_hc->improve(chr, &comp, this->m_fitfunc);
            if(chr!=popt)
            {
                diff++;
            }
            atts++;

            cout << "point located on front: " << endl;
            cout << chr << endl << endl << endl;
        }
    }
    cout << double(diff)/atts;
}

/*!
 * \brief compute fitness distance correlation 
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void fitness_distance_correlation<Chromosome,Encoding>::run()
{
    // construct the true pareto front
    list<string> truepf;
    configuration::list_parameter(keywords::TRUE_PARETO_FRONT,truepf,true);
    pareto_front<Chromosome,Encoding> pf;
    for(list<string>::iterator it=truepf.begin(); it!=truepf.end(); it++)
    {
        pf.construct_front(*it,this->m_fitfunc);
    }
    
    // construct the comparator object to use for the local search
    scalarizing_comparator<Chromosome,Encoding> comp;

    // generate several local optima and compare their distance to the
    // nearest point on the true pareto front
    unsigned int max_points;
    configuration::unsigned_integer_parameter(keywords::NUM_LOCAL_OPTIMA, max_points, true);
    unsigned int count = 0;
    while(count < max_points)
    {
        Chromosome<Encoding> chr(this->m_fitfunc);
        chr.randomize();
        chr.evaluate(this->m_fitfunc);
        comp.randomize_weights(this->m_fitfunc->objectives());
        this->m_hc->improve(chr, &comp, this->m_fitfunc);

        // now compute the nearest point on the pareto front
        if(chr.feasible())
        {
            neighborhood<Chromosome,Encoding>* nf = this->m_hc->get_neighborhood();
            typename Encoding::Genotype mindist = nf->distance_between(chr, pf[0]);
            unsigned int minidx = 0;
            for(unsigned int i=1; i<pf.size(); i++)
            {
                typename Encoding::Genotype currdist = nf->distance_between(chr, pf[i]);
                if(currdist < mindist)
                {
                    mindist = currdist;
                    minidx = i;
                }
            }

            cout << mindist << " " << chr.fitness_distance(pf[minidx]) << endl;
            count++;
        }
    }
}

/*!
 * \brief perform a random walk from nondominated solutions
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void random_walk<Chromosome,Encoding>::run()
{
    mtrandom mt;
    
    // construct the true pareto front
    list<string> truepf;
    configuration::list_parameter(keywords::TRUE_PARETO_FRONT,truepf,true);
    pareto_front<Chromosome,Encoding> pf;
    for(list<string>::iterator it=truepf.begin(); it!=truepf.end(); it++)
    {
        pf.construct_front(*it,this->m_fitfunc);
    }
    
    // how many local optima to generate
    unsigned int max_points;
    configuration::unsigned_integer_parameter(keywords::NUM_LOCAL_OPTIMA, max_points, true);

    // how many steps to take on the random walk
    unsigned int num_steps;
    configuration::unsigned_integer_parameter(keywords::RANDOM_WALK_LENGTH, num_steps, true);

    unsigned int count = 0;
    while(count < max_points)
    {
        Chromosome<Encoding> chr = pf[mt.random(pf.size())];

        cout << "0 0 " << chr.fitness << endl;
        for(unsigned int j=0; j<num_steps; j++)
        {
            int p1 = mt.random(chr.length());
            int p2 = mt.random(chr.length());
            swap(chr[p1],chr[p2]);
            chr.evaluate(this->m_fitfunc);

            // now compute the nearest point on the pareto front
            double minfdist = chr.fitness_distance(pf[0]);
            unsigned int minfidx = 0;
            typename Encoding::Genotype mingdist = this->m_hc->get_neighborhood()->distance_between(chr, pf[0]);
            unsigned int mingidx = 0;
            for(unsigned int k=1; k<pf.size(); k++)
            {
                double currfdist = chr.fitness_distance(pf[k]);
                typename Encoding::Genotype currgdist = this->m_hc->get_neighborhood()->distance_between(chr, pf[k]);
                if(currfdist < minfdist)
                {
                    minfdist = currfdist;
                    minfidx = k;
                }
                if(currgdist < mingdist)
                {
                    mingdist = currgdist;
                    mingidx = k;
                }
            }
            cout << mingdist << " " << minfdist << " " << chr.fitness << endl;
        }
        cout << endl;
        count++;
    }
}

/*!
 * \brief perform a random walk from one nondominated solution to another
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void random_walk_between_optima<Chromosome,Encoding>::run()
{
    mtrandom mt;
    
    // construct the true pareto front
    list<string> truepf;
    configuration::list_parameter(keywords::TRUE_PARETO_FRONT,truepf,true);
    pareto_front<Chromosome,Encoding> pf;
    for(list<string>::iterator it=truepf.begin(); it!=truepf.end(); it++)
    {
        pf.construct_front(*it,this->m_fitfunc);
    }
    
    // how many local optima to generate
    unsigned int max_points;
    configuration::unsigned_integer_parameter(keywords::NUM_LOCAL_OPTIMA, max_points, true);

    unsigned int count = 0;
    while(count < max_points)
    {
        // select two points from the pareto front
        Chromosome<Encoding> c1 = pf[mt.random(pf.size())];
        Chromosome<Encoding> c2 = pf[mt.random(pf.size())];
        typename Encoding::Genotype dist = this->m_hc->get_neighborhood()->distance_between(c1, c2);
        
        // take random walks from c1 to c2
        cout << "0 0 " << c2.fitness << endl;
        vector<int> order = mt.permutation(c1.length());
        for(unsigned int i=0; i<c1.length(); i++)
        {
            if(c1[order[i]] != c2[order[i]])
            {
                // find which element of c2 to swap
                unsigned int idx = 0;
                for(idx=0; idx<c2.length(); idx++)
                {
                    if(c2[order[idx]] == c1[order[i]])
                    {
                        swap(c2[order[idx]], c2[order[i]]);
                        dist--;
                        break;
                    }
                }
                c2.evaluate(this->m_fitfunc);
            }
            
            // now compute the nearest point on the pareto front
            double minfdist = c2.fitness_distance(pf[0]);
            unsigned int minfidx = 0;
            typename Encoding::Genotype mingdist = this->m_hc->get_neighborhood()->distance_between(c2, pf[0]);
            unsigned int mingidx = 0;
            for(unsigned int j=1; j<pf.size(); j++)
            {
                double currfdist = c2.fitness_distance(pf[j]);
                typename Encoding::Genotype currgdist = this->m_hc->get_neighborhood()->distance_between(c2, pf[j]);
                if(currfdist < minfdist)
                {
                    minfdist = currfdist;
                    minfidx = j;
                }
                if(currgdist < mingdist)
                {
                    mingdist = currgdist;
                    mingidx = j;
                }
            }
            cout << mingdist << " " << minfdist << " " << c2.fitness << endl;
        }

        cout << endl;
        count++;
    }
}

/*!
 * \brief compute autocorrelation along random walk 
 * 
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void ruggedness<Chromosome,Encoding>::run()
{
    mtrandom mt;
    
    // construct the true pareto front
    list<string> truepf;
    configuration::list_parameter(keywords::TRUE_PARETO_FRONT,truepf,true);
    pareto_front<Chromosome,Encoding> pf;
    for(list<string>::iterator it=truepf.begin(); it!=truepf.end(); it++)
    {
        pf.construct_front(*it,this->m_fitfunc);
    }
    
    // how many local optima to generate
    unsigned int max_points;
    configuration::unsigned_integer_parameter(keywords::NUM_LOCAL_OPTIMA, max_points, true);

    // how many steps to take on the random walk
    unsigned int num_steps;
    configuration::unsigned_integer_parameter(keywords::RANDOM_WALK_LENGTH, num_steps, true);

    unsigned int count = 0;
    while(count < max_points)
    {
        Chromosome<Encoding> chr = pf[mt.random(pf.size())];
        for(unsigned int i=0; i<chr.fitness.size(); i++)
        {
            cout << chr.fitness[i] << " ";
        }
        cout << "\n";
        
        for(unsigned int j=0; j<num_steps; j++)
        {
            int p1 = mt.random(chr.length());
            int p2 = mt.random(chr.length());
            swap(chr[p1],chr[p2]);
            chr.evaluate(this->m_fitfunc);
            for(unsigned int i=0; i<chr.fitness.size(); i++)
            {
                cout << chr.fitness[i] << " ";
            }
            cout << "\n";
        }
        cout << endl;
        count++;
    }
}

/*!
 * \brief compute ruggedness along path between pareto optima
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void ruggedness_between_optima<Chromosome,Encoding>::run()
{
    mtrandom mt;
    
    // construct the true pareto front
    list<string> truepf;
    configuration::list_parameter(keywords::TRUE_PARETO_FRONT,truepf,true);
    pareto_front<Chromosome,Encoding> pf;
    for(list<string>::iterator it=truepf.begin(); it!=truepf.end(); it++)
    {
        pf.construct_front(*it,this->m_fitfunc);
    }
    
    // how many local optima to generate
    unsigned int max_points;
    configuration::unsigned_integer_parameter(keywords::NUM_LOCAL_OPTIMA, max_points, true);

    unsigned int count = 0;
    while(count < max_points)
    {
        // select two points from the pareto front
        Chromosome<Encoding> c1 = pf[mt.random(pf.size())];
        Chromosome<Encoding> c2 = pf[mt.random(pf.size())];
        typename Encoding::Genotype dist = this->m_hc->get_neighborhood()->distance_between(c1, c2);
        
        // take random walks from c1 to c2
        vector<int> order = mt.permutation(c1.length());
        for(unsigned int i=0; i<c1.length(); i++)
        {
            if(c1[order[i]] != c2[order[i]])
            {
                // find which element of c2 to swap
                unsigned int idx = 0;
                for(idx=0; idx<c2.length(); idx++)
                {
                    if(c2[order[idx]] == c1[order[i]])
                    {
                        swap(c2[order[idx]], c2[order[i]]);
                        dist--;
                        break;
                    }
                }
                c2.evaluate(this->m_fitfunc);
            }

            // now print the new objective vector
            for(unsigned int j=0; j<c2.fitness.size(); j++)
            {
                cout << c2.fitness[j] << " ";
            }
            cout << "\n";
        }

        cout << endl;
        count++;
    }
}

/*!
 * \brief compute fraction of paths around local optima that are infeasible
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void infeasibility_region<Chromosome,Encoding>::run()
{
    // construct the true pareto front
    list<string> truepf;
    configuration::list_parameter(keywords::TRUE_PARETO_FRONT,truepf,true);
    pareto_front<Chromosome,Encoding> pf;
    for(list<string>::iterator it=truepf.begin(); it!=truepf.end(); it++)
    {
        pf.construct_front(*it,this->m_fitfunc);
    }
    
    // construct the comparator object to use for the local search
    scalarizing_comparator<Chromosome,Encoding> comp;

    // generate several local optima and compare their distance to the
    // nearest point on the true pareto front
    unsigned int max_points;
    configuration::unsigned_integer_parameter(keywords::NUM_LOCAL_OPTIMA, max_points, true);
    unsigned int count = 0;
    while(count < max_points)
    {
        Chromosome<Encoding> chr(this->m_fitfunc);
        chr.randomize();
        chr.evaluate(this->m_fitfunc);
        comp.randomize_weights(this->m_fitfunc->objectives());
        this->m_hc->improve(chr, &comp, this->m_fitfunc);

        // now compute the nearest point on the pareto front
        if(chr.feasible())
        {
            neighborhood<Chromosome,Encoding>* nf = this->m_hc->get_neighborhood();
            typename Encoding::Genotype mindist = nf->distance_between(chr, pf[0]);
            unsigned int minidx = 0;
            for(unsigned int i=1; i<pf.size(); i++)
            {
                typename Encoding::Genotype currdist = nf->distance_between(chr, pf[i]);
                if(currdist < mindist)
                {
                    mindist = currdist;
                    minidx = i;
                }
            }

            typename Encoding::Genotype totaldist;
            typename Encoding::Genotype infeasdist;
            nf->feasible_distance_between(chr, pf[minidx], totaldist, infeasdist, this->m_fitfunc);
            cout << mindist << " " << infeasdist << " " << chr.fitness_distance(pf[minidx]) << endl;
            count++;
        }
    }
}

/*!
 * \brief override to prevent initialization of useless landscape components
 *
 * \author deong
 * \date 06/06/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void hypervolume_analysis<Chromosome,Encoding>::initialize()
{
    // read in the objective function information
    this->m_fitfunc = Encoding::ProblemFactoryType::construct();
    
    // initialize the encoding parameters
    Encoding::initialize_parameters(this->m_fitfunc);
}

/*!
 * \brief compute the hypervolume of a set of pareto fronts after the fact
 *
 * \author deong
 * \date 06/06/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
void hypervolume_analysis<Chromosome,Encoding>::run()
{
    // fill the pareto front
    string runfile;
    configuration::string_parameter(keywords::TRUE_PARETO_FRONT,runfile,true);
    pareto_front<Chromosome,Encoding> pf;
    pf.construct_front(runfile, this->m_fitfunc);
    
    // compute the hypervolume for the front
    hypervolume<Chromosome,Encoding> hv;
    hv.initialize();
    hv.generation_completed(pf);
    hv.compute();
    hv.report(cout);
}

/*!
 * \brief initialize the plateau tool
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <template <typename> class Chromosome, typename Encoding>
void pareto_plateaus<Chromosome,Encoding>::initialize()
{
}

/*!
 * \brief sort the population into nonoverlapping fronts
 *
 * \author deong@acm.org
 * \date 10/15/2008
 */
template <template <typename> class Chromosome, typename Encoding>
void pareto_plateaus<Chromosome,Encoding>::nondominated_sort(population<Chromosome,Encoding>& pop,
                                                             vector<pareto_front<Chromosome,Encoding> >& fronts)
{
    // clear any previous sorting
    fronts.clear();

    // sort the population into separate fronts
    unsigned int fnum = 0;
    while(pop.size() > 0)
    {
        pareto_front<Chromosome,Encoding> curr;
        for(unsigned int i=0; i<pop.size(); ++i)
        {
            curr.add(pop[i]);
        }
        fronts.push_back(curr);
        pop -= curr;
        fnum++;
    }
}

/*!
 * \brief run the plateaus tool
 *
 * \author deong@acm.org
 * \date 10/14/2008
 */
template <template <typename> class Chromosome, typename Encoding>
void pareto_plateaus<Chromosome,Encoding>::run()
{
    string infile;
    configuration::string_parameter(keywords::TRUE_PARETO_FRONT,infile,true);
    ifstream in(infile.c_str());
    if(!in)
    {
        cout << "error opening dump file" << endl;
        return;
    }

    population<Chromosome,Encoding> pop;
    in >> pop;
    
    //! do a nondominated sort on the population
    population<Chromosome,Encoding> saved_pop=pop;
    vector<pareto_front<Chromosome,Encoding> > fronts;
    nondominated_sort(pop,fronts);
    
    //! and label each solution with its rank
    map<Chromosome<Encoding>,int> ranks;
    for(unsigned int i=0; i<fronts.size(); ++i)
    {
        cout << fronts[i] << endl << endl;
        for(unsigned int j=0; j<fronts[i].size(); ++j)
        {
            ranks[fronts[i][j]]=i;
        }
    }

    //! keep a record of how many pathways go between each pair of ranks
    //! (the "+1" is to make a special rank for infeasible solutions)
    const int infeasible=fronts.size();
    vector<vector<int> > tableaux(fronts.size()+1);
    for(unsigned int i=0; i<fronts.size()+1; ++i)
        tableaux[i].resize(fronts.size()+1);
    
    for(unsigned int i=0; i<saved_pop.size(); ++i)
    {
        Chromosome<Encoding> sol=saved_pop[i];
        neighborhood_factory<Chromosome,Encoding> nf;
        neighborhood<Chromosome,Encoding>* hood=nf.construct();
        hood->initialize(sol);
        while(hood->has_more_neighbors())
        {
            int src=ranks[sol];
            int dest;
            
            //! compute the next neighbor of the current point and figure out
            //! what rank it is
            move<Chromosome,Encoding> m=hood->next_neighbor();
            Chromosome<Encoding> tmp=sol;
            m.apply(tmp);
            if(tmp==sol)        // don't consider moves that did nothing
                continue;
            typename map<Chromosome<Encoding>,int>::iterator it=ranks.find(tmp);
            if(it!=ranks.end())
                dest=it->second;
            else
                dest=infeasible;
            
            tableaux[src][dest]++;
        }
    }
    
    // write out the transition matrix
    for(unsigned int i=0; i<tableaux.size(); ++i)
    {
        for(unsigned int j=0; j<tableaux.size(); ++j)
        {
            cout << setw(6) << tableaux[i][j];
        }
        cout << endl;
    }

    //! write a graphviz file to generate the graphs
    string graphviz_file;
    configuration::string_parameter(keywords::GRAPHVIZ_FILE,graphviz_file,true);
    ofstream out(graphviz_file.c_str());
    if(!out)
    {
        cout << "could not open graphviz file" << endl;
        return;
    }
    out << "digraph plateau {" << endl;
    for(unsigned int i=0; i<tableaux.size(); ++i)
    {
        int total=0;
        for(unsigned int j=0; j<tableaux.size(); ++j)
        {
            total+=tableaux[i][j];
        }
        
        for(unsigned int j=0; j<tableaux.size(); ++j)
        {
            if(tableaux[i][j]>0)
            {
                out << "  R" << i << " -> R" << j << " [ label = \"";
                out << setprecision(2) << static_cast<double>(tableaux[i][j])/total;
                out << "\" ];" << endl;
            }
        }
    }
    out << "}" << endl;
}

/*!
 * \brief create and return initialized landscape operators
 *
 * \author deong
 * \date 05/09/2007
 *
 * \code
 * Modification History
 * MM/DD/YYYY   DESCRIPTION
 * \endcode
 */
template <template <typename> class Chromosome, typename Encoding>
landscape<Chromosome,Encoding>* landscape_factory<Chromosome,Encoding>::construct()
{
    // what tool to utilize
    string toolname;
    configuration::string_parameter(keywords::LANDSCAPE_TOOL, toolname, true);
    if(toolname == keywords::PERTURBATION_SEARCH)
    {
        landscape<Chromosome,Encoding>* lnd = new perturbation_search<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname == keywords::FITNESS_DISTANCE_CORRELATION)
    {
        landscape<Chromosome,Encoding>* lnd = new fitness_distance_correlation<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname == keywords::RANDOM_WALK)
    {
        landscape<Chromosome,Encoding>* lnd = new random_walk<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname == keywords::RANDOM_WALK_BETWEEN_OPTIMA)
    {
        landscape<Chromosome,Encoding>* lnd = new random_walk_between_optima<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname == keywords::RUGGEDNESS)
    {
        landscape<Chromosome,Encoding>* lnd = new ruggedness<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname == keywords::RUGGEDNESS_BETWEEN_OPTIMA)
    {
        landscape<Chromosome,Encoding>* lnd = new ruggedness_between_optima<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname == keywords::INFEASIBILITY_REGION)
    {
        landscape<Chromosome,Encoding>* lnd = new infeasibility_region<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname==keywords::HYPERVOLUME)
    {
        landscape<Chromosome,Encoding>* lnd=new hypervolume_analysis<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else if(toolname==keywords::PARETO_PLATEAUS)
    {
        landscape<Chromosome,Encoding>* lnd=new pareto_plateaus<Chromosome,Encoding>;
        lnd->initialize();
        return lnd;
    }
    else
    {
        cerr << "illegal landscape tool selected: " << toolname << endl;
        exit(1);
    }
}
