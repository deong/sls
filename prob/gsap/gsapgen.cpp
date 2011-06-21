#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ios>
#ifndef WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif
#include <boost/random.hpp>
#include <boost/program_options.hpp>

using namespace std;
namespace po=boost::program_options;

int main(int argc, char** argv)
{
    boost::uint32_t seed;
    int agents=0;
    int tasks=0;
    int obj=0;
    int days=0;
    int slots=0;
    double resource_mu=0;
    double resource_sigma=0;
    double cost_mu=0;
    double cost_sigma=0;
    double tasks_per_slot_mu=0;
    double tasks_per_slot_sigma=0;
    int min_tasks_per_slot=2;
    double agent_task_prob=0;
    double slack=0;
    boost::mt19937 rng;

    po::options_description opt("Allowed options");
    opt.add_options()
        ("help,h", "display program usage information")
        ("agents,n", po::value<int>(&agents), "number of agents")
        ("tasks,m", po::value<int>(&tasks), "number of task classes")
        ("objectives,k", po::value<int>(&obj), "number of objectives")
        ("days,d", po::value<int>(&days), "number of days to schedule")
        ("slots,s", po::value<int>(&slots), "number of time slots per day")
        ("resource_mu", po::value<double>(&resource_mu)->default_value(50),
         "mean resources for agent/task pairs")
        ("resource_sigma", po::value<double>(&resource_sigma)->default_value(20),
         "standard deviation of agent/task resources")
        ("cost_mu", po::value<double>(&cost_mu)->default_value(50),
         "mean cost for agent/task pairs")
        ("cost_sigma", po::value<double>(&cost_sigma)->default_value(20),
         "standard deviation of agent/task costs")
        ("tasks_per_slot_mu", po::value<double>(&tasks_per_slot_mu)->default_value(5),
         "mean number of tasks per time slot")
        ("tasks_per_slot_sigma", po::value<double>(&tasks_per_slot_sigma)->default_value(2),
         "standard deviation of tasks per time slot")
        ("agent_task_prob", po::value<double>(&agent_task_prob)->default_value(0.5),
         "probability agent/task capability")
        ("slack", po::value<double>(&slack)->default_value(0.25),"looseness of constraints");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opt), vm);
    po::notify(vm);

    if(vm.count("help") || agents==0 || tasks==0 || obj==0 || days==0 || slots==0)
    {
        cout << opt << endl;
        return EXIT_FAILURE;
    }

#ifndef WIN32    
    // try to read from /dev/random; if that fails, use the system time
    int fd=open("/dev/random",O_RDONLY|O_NONBLOCK);
    if((fd!=-1) && (read(fd,&seed,sizeof(unsigned int))==sizeof(unsigned int)))
    {
        seed=static_cast<boost::uint32_t>(seed);
    }
    else
    {
        cerr << "failed to read random bytes from /dev/random...falling back to system clock" << endl;
        seed = static_cast<boost::uint32_t>(time(NULL));
    }
#else
    seed = static_cast<boost::uint32_t>(time(NULL));
#endif  
    rng.seed(seed);
    
    cout << ";; agents: " << agents << ";   tasks: " << tasks << ";   objectives: " <<
        obj << ";   days: " << days << ";   slots: " << slots << ";   resource_mu: " <<
        resource_mu << ";   resource_sigma: " << resource_sigma << ";   cost_mu: " <<
        cost_mu << ";   cost_sigma: " << cost_sigma << ";   tasks_per_slot_mu: " <<
        tasks_per_slot_mu << ";   tasks_per_slot_sigma: " << tasks_per_slot_sigma <<
        ";   agent_task_prob: " << agent_task_prob << ";   slack: " << slack << "\n\n";

    // write first few lines of problem information
    cout << agents << "\n" << tasks << "\n" << obj << "\n" << days << "\n" << slots << "\n\n";

    // compute appropriate capacity constraint distribution and generate random values
    double mean_cap=(days*slots*tasks_per_slot_mu*resource_mu)/(agents*agent_task_prob)*(1+slack);
    cerr << "mean_cap = " << mean_cap << endl;
    boost::normal_distribution<> capdist(mean_cap, 0.1*mean_cap);
    boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > capgen(rng,capdist);
    for(int i=0; i<agents; ++i)
        cout << static_cast<int>(capgen()+0.5) << " ";
    cout << "\n\n";

    // generate random tasks for each time slot
    boost::normal_distribution<> numtasksdist(tasks_per_slot_mu, tasks_per_slot_sigma);
    boost::uniform_int<> tasksdist(0, tasks-1);
    boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > numtasksgen(rng, numtasksdist);
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > tasksgen(rng, tasksdist);
    for(int i=0; i<days; ++i)
    {
        for(int j=0; j<slots; ++j)
        {
            int num_tasks = max(min_tasks_per_slot, static_cast<int>(numtasksgen()+0.5));
            for(int k=0; k<num_tasks; ++k)
                cout << tasksgen() << " ";
            cout << "\n";
        }
        cout << "\n";
    }
    
    // generate random resource allocation matrix
    // note that I'm generating the matrix in memory rather than writing it out directly
    // because I need to generate a certain number of zeros to indicate agents who cannot
    // perform certain task classes, and I want to match the zeros in the cost matrices
    // up with the ones in the resource matrix.
    boost::uniform_real<> availdist(0.0, 1.0);
    boost::normal_distribution<> resdist(resource_mu, resource_sigma);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > availgen(rng, availdist);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > resgen(rng, resdist);
    vector<vector<int> > resources;
    resources.resize(agents);
    for(int i=0; i<agents; ++i)
    {
        resources[i].resize(tasks);
        for(int j=0; j<tasks; ++j)
        {
            if(availgen()<agent_task_prob)
                resources[i][j] = max(1, static_cast<int>(resgen()+0.5));
            else
                resources[i][j] = 0;
            cout << setw(5) << resources[i][j];
        }
        cout << "\n";
    }
    cout << "\n";

    // generate random cost matrices
    boost::normal_distribution<> costdist(cost_mu, cost_sigma);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > costgen(rng, costdist);
    vector<vector<vector<int> > > costs;
    costs.resize(obj);
    for(int k=0; k<obj; ++k)
    {
        costs[k].resize(agents);
        for(int i=0; i<agents; ++i)
        {
            costs[k][i].resize(tasks);
            for(int j=0; j<tasks; ++j)
            {
                if(resources[i][j]!=0)
                    costs[k][i][j] = max(1, static_cast<int>(costgen()+0.5));
                else
                    costs[k][i][j] = 0;
                cout << setw(5) << costs[k][i][j];
            }
            cout << "\n";
        }
        cout << "\n";
    }
    cout << endl;
    
    return EXIT_SUCCESS;
}
