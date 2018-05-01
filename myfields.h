#ifndef MYFIELDS_H
#define MYFIELDS_H

#include "mytypetraits.h"

#include <string>
#include <optional>
#include <type_traits>
namespace grammar {



template <typename T>
struct argument
{
    typedef T  argument_type;
    typedef std::optional<std::decay_t<argument_type>> default_type;
    std::string idField;
    default_type default_value;
    argument(C<T>,const char* id):idField{id}, default_value{}{}
    argument(C<T>,const char* id, T val):idField{id}, default_value{std::move(val)}{}
 };


}
template <auto>
struct function_trait{};





namespace grammar {


template <class Object,class Method>
struct field
{
    typedef
    std::invoke_result_t<Method,Object> return_type;

    typedef Method member_type;
    typedef Object object_type;
    typedef  std::decay_t<return_type> result_type;
     typedef std::optional<std::decay_t<return_type>> default_type;

     std::string idField;
     member_type access_method;
     default_type default_value;
     field(C<Object>,const char* id, member_type get):idField{id},access_method{get}, default_value{}{}
};



template <class Object,class... Method>
bool has_all(const std::tuple<field<Object,Method>...>& fs)
{
    return std::apply([](auto&...x){return (x.default_value.has_value()&&...);},fs);
}









} // namespace grammar



#endif // MYFIELDS_H
