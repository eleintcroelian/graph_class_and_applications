#include "MinStack.hpp"
#include <cassert>
#include <iostream>
#include "reshape.hpp"
#include "symmetric.hpp"
#include "transpose.hpp"
#include "intersection.hpp"

int main()
{
    MinStack<int> s;
    s.push(5);
    s.push(2);
    s.push(100);
    assert(s.size() == 3);
    assert(s.min() == 2);
    assert(s.top() == 100);
    s.pop();
    // std::cout<<s.size()<<std::endl;
    assert(s.size() == 2);
    assert(s.top() == 2);
    s.pop();
    assert(s.size() == 1);
    // std::cout << s.min() << std::endl;
    assert(s.min() == 5);

    MinStack<int> k(s);
    assert(k.min() == 5);
    assert(k.size() == 1);
}
