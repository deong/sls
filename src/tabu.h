/*!
 * \file tabu.h
 *
 * tabu search algorithms
 *
 * Deon Garrett
 * jdgarrett@gmail.com
 */

#ifndef _TABU_H_
#define _TABU_H_

#include <map>
#include <vector>
#include <utility>
#include "sls.h"
#include "chromosome.h"
#include "encoding.h"
#include "problems.h"
#include "comparator.h"
#include "localsearch.h"
#include "move.h"

using namespace std;

/*!
 * \class tabu_list
 */
template <template <typename> class Chromosome, typename Encoding>
class tabu_list
{
public:
    typedef pair<move<Chromosome,Encoding>,unsigned int> tlist_item;

public:
    tabu_list();
    ~tabu_list();
    void set_prefix(const string& prefix);
    void initialize();
    void clear();
    void accept_move(const Chromosome<Encoding>& chr, const move<Chromosome,Encoding>& m, unsigned int iter);
    bool tabu(const move<Chromosome,Encoding>& m, unsigned int iter) const;

private:
    unsigned int _ttmin;
    unsigned int _ttmax;
    deque<tlist_item> _internal;
    string m_prefix;
};

/*!
 * \class tabu_search
 */
template <template <typename> class Chromosome, typename Encoding>
class tabu_search : public local_search<Chromosome,Encoding>
{
public:
    tabu_search();
    virtual ~tabu_search();
    virtual void set_prefix(const string& prefix);
    virtual void initialize();
    virtual void reset();
    virtual void improve(Chromosome<Encoding>& chr,
                         comparator<Chromosome,Encoding>* comp,
                         const typename Encoding::ProblemType* prob);

protected:
    tabu_list<Chromosome,Encoding>* _tlist;
};

#include "tabu.cpp"

#endif
