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


    typedef grammar::CommandManager<Objects> CM;
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
