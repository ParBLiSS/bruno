/**
 * @file    function_tratis.hpp
 * @ingroup utils
 * @author  Tony Pan <tpan7@gatech.edu>
 * @brief   Implements helper functions for fucntors.
 *
 * Copyright (c) 2014 Georgia Institute of Technology
 *
 * TODO add Licence
 */

#ifndef FUNCTION_TRAITS_HPP_
#define FUNCTION_TRAITS_HPP_

namespace bliss
{
  namespace functional
  {

    //==== functor and function pointer traits support.
    /**
     * @class   function_traits
     * @brief function traits to get return type.
     * Supports:
     *   overloaded operator() and functions : via parameter types.
     *     for functors, fortunately we can construct the list from the Base_Iterator (our use case).
     *     for function pointers, (member, nonmember) these are already explicitly specified by the user.
     *   static function and member function pointers.
     *   for overloaded function, need to construct the function pointer, which requires specifying the exact function input parameters for decltype.
     *        use static_cast<RET (*)(ARGS)>
     *        else create a typedef for the function pointer type, assign the function to it (compiler resolve it) and then do decltype on that.
     *   works with lambda functions (except that class type may be main function.  don't know how this affects evaluation.
     *   functor input const/ref:  user/iterator specified arg types, which may not match correctly.  when calling, need to supply
     *       references instead, then let the function call cast to const, ref, const ref, or just make copy.
     *       result_of appears to be smart enough to deal with const/ref during decltype step.
     *  function returning const type:  this is considered not useful and obselete with c++11.  not tested here.
     *
     * Possible Issues:
     *   if operator() overloading is based on input parameter const/ref, then decltype inside result_of will not be able to resolve.  does not happen with func ptrs.
     *
     * Use:  'dereference' function is implemented in the func_traits class as static functions.
     *
     */
    template<typename F, typename ... Args>
    struct function_traits
    {
        typedef decltype(std::declval<F>()(std::declval<typename std::add_lvalue_reference<Args>::type>()...)) return_type;
    };

  } /* namespace functional */
} /* namespace bliss */
#endif /* FUNCTION_TRAITS_HPP_ */
