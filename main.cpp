#include <iostream>
//#include "myCommandManagement.h"

#include "qmodel.h"
#include "mygrammar.h"
#include "mycompilation.h"
#include "mynewcommandmanager.h"
#include "myscriptmanager.h"
#include "Experiment.h"
#include "Markov.h"
#include "simulation.h"
#include "commands.h"
int main(int argc, char **argv)
{
    std::cerr<<argv[0]<<"\n";
    std::cerr<<argv[1]<<"\n";

    std::map<double, std::string> p{{0.5,"gsdgfs"},{9,"gsdhse"}};
    std::cout<<p;
   // std::cin>>p;
    std::cout<<p;

    std::cout<<std::tuple(1.4566,"qr");

    auto A=Allosteric_Model
            (
    {"R","L","R","L","R","L"},//const std::vector<std::string>& conformational_changes,
    {{{"R",true},"beta"},{{"R",false},"alpha"},{{"L",false},"koff"},{{"L",true},"kon"}},//const std::map<std::pair<std::string,bool>,std::string> conformational_changes_names,
    {1,3,5},//const std::set<std::size_t>& agonist_changes,
    {0,2,4},//const std::set<std::size_t>& conductance_changes,
    {{0u,"g_0"},{1u,"g_1"},{2u,"g_2"},{3u,"g_3"}},//const std::map<std::size_t, std::string>& conductance_names,

    {
                    {0,{{1},{"RL","RL_R"}}},
                    {0,{{2},{"RR","RR_I"}}},
                    {0,{{4},{"RR","RR_D"}}},
                    {0,{{5},{"LR","LR_R"}}},
                    {0,{{1,2},{"RLR","RLR_I"}}},
                    {0,{{2,4},{"RRR","RRR_M"}}},
                    {0,{{4,5},{"RLR","RLR_D"}}},

                    {1,{{0},{"RL","RL_L"}}},
                    {1,{{2},{"LR","LR_L"}}},
                    {1,{{0,2},{"RLR","RLR_L"}}},

                }//const std::multimap<std::size_t, std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>& conformational_interactions,
                );

    std::map<std::string, double> P
        {
            {"LR",300}	,{"LR_L",0},	{"LR_R",1},	{"RL",30},	{"RLR",300},	{"RLR_D",1},	{"RLR_I",1},	{"RLR_L",0},	{"RL_L",0},	{"RL_R",1},	{"RR",3},	{"RRR",1},	{"RRR_M",0.5},	{"RR_D",0.5},	{"RR_I",0.5},	{"alpha",1e7}, 	{"beta",10},	{"g_0",0.0},{"g_1",-0.001},	{"g_2",-0.1}	,{"g_3",-1},	{"koff",1e7},	{"kon",500}
        };
    SingleLigandModel SM(A.Qs(P),A.g(P));

    write(std::cout,A);
    std::ofstream f;
    f.open("output.txt");
    write(f,A);
    f.close();
    std::ifstream fi;
    fi.open("output.txt");
    Allosteric_Model B{};
    read(fi,B);
    write(std::cout, B);
    std::ifstream fe;
    fe.open("/home/lmoffatt/Code/macro-dr/Data/Moffatt_Hume_2007_ATP.txt");
    io::myDataFrame<double> da;
    da.read(fe);
    da.write(std::cout);



    auto e=experiment::DataFrame_to_Experiment(da,"t","ns","xATP","yCurrent", 50E3);
    Markov_Model_calculations<SingleLigandModel,experiment::basic_Experiment<experiment::point<double,double>>,double> MC(SM,1000,e);


    auto s=simulate<decltype (e),Allosteric_Model>()(0,e,A,P,1000,10);
    auto ds=experiment::Experiment_to_DataFrame(s);
    ds.write(std::cout);



   std::cout<<std::endl;
    std::map<std::string, double> beta
    {
        {"LR",3}	,{"LR_L",0.5},	{"LR_R",0.5},	{"RL",3},	{"RLR",3},	{"RLR_D",0.5},	{"RLR_I",0.5},	{"RLR_L",0.5},	{"RL_L",0.5},	{"RL_R",0.5},	{"RR",3},	{"RRR",3},	{"RRR_M",0.5},	{"RR_D",0.5},	{"RR_I",0.5},	{"alpha",1e3},	{"beta",1},	{"g_1",0.01},	{"g_2",0.3}	,{"g_3",1},	{"koff",1e5},	{"kon",100}
    };


    auto Q=A.Qs(P);


    typedef experiment::basic_Experiment<experiment::point<double,double>> singleLigandExperiment;



    typedef grammar::CommandManager<typename Objects::types,typename Objects::commands> CM;
//    typedef grammar::CommandManager<> CM;

    CM cm;

 //   cm.insert_function(C<Allosteric_Model>{},Allosteric_Model::name.c_str(),Constructor<Allosteric_Model>{},Allosteric_Model::get_constructor_fields());

    /* cmd_["readIt"]=new readCommand(this);
      cmd_["align"]=new AlignCommand(this);
      cmd_["writeIt"]=new writeCommand(this);
      cmd_["merge"]=new MergeCommand(this);
      cmd_["distances"]=new DistancesCommand(this);
      cmd_["histogram"]=new HistogramCommand(this);
      cmd_["simulate"]=new SimulateCommand(this);
      cmd_["experiment"]=new ExperimentCommand(this);
      cmd_["likelihood"]=new LikelihoodCommand(this);
      cmd_["optimize"]=new OptimizeCommand(this);
      cmd_["evidence"]=new EvidenceCommand(this);
      cmd_["tempering"]=new TemperingCommand(this);
    */
    /*
    cm.push_command("readIt",
                    make_Command
                    (C<CM>(),C<void>(),Co<Cls>(),readIt<CM>(),
                     std::pair<CM*,std::string>{&cm,"CommandManager"},
                     std::pair<std::string,std::string>{"","filename"},
                     std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"}));

    cm.push_command("writeIt",
                    make_Command
                    (C<CM>(),C<bool>(),Co<Cls>(),writeIt<CM>(),
                     std::pair<CM*,std::string>{&cm,"CommandManager"},
                     std::pair<std::string,std::string>{"","id"},
                     std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"},
                     std::pair<std::string,std::string>{"","filename"},
                     std::pair<bool,std::string>{false,"append"}
                     ));
    cm.push_command("dataFrame",
                    make_Command
                    (C<CM>(),C<bool>(),Co<Cls>(),dataFrame<CM>(),
                     std::pair<CM*,std::string>{&cm,"CommandManager"},
                     std::pair<std::string,std::string>{"","id"},
                     std::pair<std::ostream*,std::string>{&std::cerr,"log_stream"},
                     std::pair<std::string,std::string>{"","filename"}
                     ));


*/


    switch (argc) {
    case 1:

        break;
    case 2:
    {
        myScript<CM> s(&cm);
        return s.run(argv[1], std::cerr);
    }
        break;
    case 4:
    {
        myScript<CM> s(&cm);
        return s.runDefine(argv[1],{argv[2]},{argv[3]}, std::cerr);
    }
    break;

    default:
    {
        std::size_t nlabels=argc/2-1;
        std::vector<std::string> labels(nlabels);
        std::vector<std::string> valuesInPlace(nlabels);
        for (std::size_t i=0; i<nlabels; ++i)
        {
            labels[i]=argv[2*(i+1)];
            valuesInPlace[i]=argv[2*(i+1)+1];
        }
        myScript<CM> s(&cm);
        s.runDefine(argv[1],labels,valuesInPlace, std::cerr);


    }
        break;




        return 0;
}

}
