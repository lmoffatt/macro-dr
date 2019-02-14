#include "mynewcommandmanager.h"


template<class Objects>
template<class... type>
void grammar::CommandManager<Objects>::insert_command(Cs<type...>)
{
    (insert_command<type>(),...);
}



#include "commands.h"

template  void grammar::CommandManager<Objects>::insert_command<>(typename grammar::CommandManager<Objects>::myCommands);

