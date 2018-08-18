#include <iostream>
//#include "myCommandManagement.h"

#include "qmodel.h"
#include "mygrammar.h"
#include "mycompilation.h"
#include "mynewcommandmanager.h"
#include "myscriptmanager.h"
#include "Experiment.h"
#include "commands.h"
#include "myevidence.h"
int main(int argc, char **argv)
{
    std::cerr<<argv[0]<<"\n";
    std::cerr<<argv[1]<<"\n";

    std::map<double, std::string> p{{0.5,"gsdgfs"},{9,"gsdhse"}};
    std::cout<<p;
   // std::cin>>p;
    std::cout<<p;

    std::cout<<std::tuple(1.4566,"qr");

    Conformational_change_label L("L");
    Conformational_change_label R("R");

    Conformational_change binding(1, 0, Conformational_change_label("L"),  rate_Parameter_label("kon"),  rate_Parameter_label("koff"));
    Conformational_change rocking(0, 1, Conformational_change_label("R"),  rate_Parameter_label("beta"),  rate_Parameter_label("alpha"));


    Conformational_interaction RL({R, L}, Coupling_factor_Parameter_label( "RL"),
    {Coupling_coefficient_Parameter_label("RL_R"),Coupling_coefficient_Parameter_label("RL_L")});
    Conformational_interaction LR({L, R}, Coupling_factor_Parameter_label( "LR"),
    {Coupling_coefficient_Parameter_label("LR_L"),Coupling_coefficient_Parameter_label("LR_R")});
    Conformational_interaction RR({R, R}, Coupling_factor_Parameter_label( "RR"),
    {Coupling_coefficient_Parameter_label("RR_1"),Coupling_coefficient_Parameter_label("RR_2")});
     Conformational_interaction RLR({R, L, R}, Coupling_factor_Parameter_label( "RLR"),
    {Coupling_coefficient_Parameter_label("RLR_1"),Coupling_coefficient_Parameter_label("RLR_3"),Coupling_coefficient_Parameter_label("RLR_L")});
    Conformational_interaction RRR({R, R, R}, Coupling_factor_Parameter_label( "RRR"),
    {Coupling_coefficient_Parameter_label("RRR_R"),Coupling_coefficient_Parameter_label("RRR_R"),Coupling_coefficient_Parameter_label("RRR_R")});

    Allosteric_Model A(3, {{L, binding},{R, rocking}}, {R,L}, {RL,LR,RR,RLR,RRR},
    {{0,Conductance_Parameter_label("g_0")},{1,Conductance_Parameter_label("g_1")},{2,Conductance_Parameter_label("g_2")},{3, Conductance_Parameter_label("g_3")}});


 /*  auto A=Allosteric_Model
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
*/

    std::map<std::string, double> par=
            {
                {"LR",300}	,{"LR_L",0},	{"LR_R",1},	{"RL",30},	{"RLR",300},	{"RLR_1",1},	{"RLR_3",1},	{"RLR_L",0},	{"RL_L",0},	{"RL_R",1},	{"RR",3},	{"RRR",1},	{"RRR_R",0.5},	{"RR_1",0.5},	{"RR_2",0.5},	{"alpha",1e7}, 	{"beta",10},	{"g_0",0.0},{"g_1",-0.001},	{"g_2",-0.1}	,{"g_3",-1},	{"koff",1e7},	{"kon",500}, {"Number_of_Channels" , 100} , {"gaussian_noise" , 1.0 }};

    Parameters_values<Allosteric_Model> P(par);
    SingleLigandModel SM(A.Qs(P),A.g(P),1000,0);

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
    Markov_Model_calculations<Markov_Transition_step,SingleLigandModel,experiment::basic_Experiment<experiment::point<double,double>,measure_just_y<double>>,double> MC(SM,e);


    auto s=simulate<decltype (e),Allosteric_Model>::run(0,e,A,P,10);
    auto ds=experiment::Experiment_to_DataFrame(s);
    ds.write(std::cout);
    write(std::cout,s);


   std::cout<<std::endl;
    std::map<std::string, double> beta
    {
        {"LR",3}	,{"LR_L",0.5},	{"LR_R",0.5},	{"RL",3},	{"RLR",3},	{"RLR_D",0.5},	{"RLR_I",0.5},	{"RLR_L",0.5},	{"RL_L",0.5},	{"RL_R",0.5},	{"RR",3},	{"RRR",3},	{"RRR_M",0.5},	{"RR_D",0.5},	{"RR_I",0.5},	{"alpha",1e3},	{"beta",1},	{"g_1",0.01},	{"g_2",0.3}	,{"g_3",1},	{"koff",1e5},	{"kon",100},
      {"Number_of_Channels" , 100} , {"gaussian_noise" , 1.0 }};


    auto Q=A.Qs(P);


    //typedef experiment::basic_Experiment<experiment::point<double,double>> singleLigandExperiment;





    typedef grammar::CommandManager<typename Objects::types,typename Objects::commands, typename Objects::templateCommands> CM;
//    typedef grammar::CommandManager<> CM;
    static_assert (!has_base_type<grammar::Compiled_Statement<CM>>::value,"" );

    CM cm;

    //typedef typename optional_tag_t<grammar::Compiled_Statement<CM>*>::test test;

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
