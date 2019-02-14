#include "mynewcommandmanager.h"


template<class Objects>
grammar::CommandManager<Objects>::CommandManager():d_{}
{
    insert_constructor(myTypes());
    insert_valuer(myTypes());
    insert_loader(myTypes());
    insert_Literal(myTypes());
    insert_command(myCommands());
}


#include "commands.h"

template class grammar::CommandManager<Objects>;
