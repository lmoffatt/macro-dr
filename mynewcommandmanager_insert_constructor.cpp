#include "mynewcommandmanager.h"


template<class Objects>
template<class... type>
void grammar::CommandManager<Objects>::insert_constructor(Cs<type...>)
{
    (insert_constructor<type>(),...);
}

template<class Objects>
template<class... type>
void grammar::CommandManager<Objects>::insert_valuer(Cs<type...>)
{
    (insert_valuer<type>(),...);
}

template<class Objects>
template<class... type>
void grammar::CommandManager<Objects>::insert_Literal(Cs<type...>)
{
    (insert_Literal<type>(),...);
}


template<class Objects>
template<class... type>
void grammar::CommandManager<Objects>::insert_loader(Cs<type...>)
{
    (insert_loader<type>(),...);
}



#include "commands.h"


template  void grammar::CommandManager<Objects>::insert_constructor<>(typename grammar::CommandManager<Objects>::myTypes);

template  void grammar::CommandManager<Objects>::insert_loader<>(typename grammar::CommandManager<Objects>::myTypes);
template  void grammar::CommandManager<Objects>::insert_valuer<>(typename grammar::CommandManager<Objects>::myTypes);
template  void grammar::CommandManager<Objects>::insert_Literal<>(typename grammar::CommandManager<Objects>::myTypes);

