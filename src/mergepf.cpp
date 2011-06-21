#include <iostream>
#include "chromosome.h"
#include "encoding.h"
#include "problems.h"
#include "pfront.h"
#include "configuration.h"
#include "keywords.h"

using namespace std;

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        cerr << "usage: mergepf <config>" << endl;
        exit(EXIT_FAILURE);
    }

    configuration::read_configuration_file(argv[1]);
    permutation_encoding::ProblemType* prob = permutation_encoding::ProblemFactoryType::construct();

    string outname;
    configuration::string_parameter(keywords::TRUE_PARETO_FRONT,outname,true);
    pareto_front<chromosome,permutation_encoding> pf;
    pf.construct_front(outname, prob);

    for(unsigned int i=0; i<pf.size(); i++)
    {
        cout << pf[i].fitness << endl;
    }

    return 0;
}

