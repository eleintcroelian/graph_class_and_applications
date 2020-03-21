#include <cassert>
#include <algorithm>
#include <iostream>
// Copy and move assignments and copy and move constructors were influenced by the lecture notes 10

/** Dynamically memory allocating array class. dynamically allocated, built-in arrays that grow
* in size based on user's needs.
* push() inserts a new element onto the stack with O(1) complexity.
* pop() removes the most recently inserted element onto the stack, if it is not empty, in O(1)
* top() returns the value of the most recently inserted element onto the stack, if it is not empty,
* in O(1)
* size() returns the current number of elements O(1)
*
* MinStack keeps track of the minimum element as the array is being built.
* min () returns the minimum value, if it is not empty, in O(1) time.
*/
template <class T>
class MinStack
{
private:
    T *array_pointer; // Pointer to the array on stack
    int array_length; // Stack length of array
    int elem_num;     // Number of elements stored in array
    T current_min;    // Current minimum element
public:
    // Default Constructor
    MinStack()
    {
        array_pointer = new T[5];
        for (int i = 0; i < 5; i++)
            array_pointer[i] = 0;

        array_length = 5;
        elem_num = 0;
        current_min = std::numeric_limits<T>::max();
    }

    // Copy Constructor
    MinStack(const MinStack &m) : array_pointer{new T[m.size()]},
                                  array_length{m.array_length},
                                  elem_num{m.elem_num},
                                  current_min{m.current_min}

    {
        for (int i = 0; i < elem_num; i++)
        {
            array_pointer[i] = m.array_pointer[i];
        }
    }

    // Move Constructor
    MinStack(MinStack &&a) // Move constructor .
        : array_length{a.array_length},
          array_pointer{a.array_pointer},
          elem_num{a.elem_num},
          current_min{a.current_min}
    {
        a.elem_num = 0;
        a.array_length = 0;
        a.array_pointer = nullptr;
    }
    // Destructor

    ~MinStack()
    {
        delete[] array_pointer;
    }
    // Copy assignment

    MinStack &operator=(const MinStack &m)
    {
        assert(m.size() == this->size());
        std::copy(m.array_pointer, m.array_pointer + m.size(), array_pointer); // Copy elements .
        current_min(m.current_min);
    }

    // Move assignment
    MinStack<T> &operator=(MinStack &&a)
    {
        swap(current_min, a.current_min);
        swap(elem_num, a.elem_num);
        swap(array_pointer, a.array_pointer);
        return *this;
    }

    void push(T val)
    {
        if (elem_num == array_length)
        {
            T *temp;
            array_length = array_length * 2;
            temp = new T[array_length];
            for (int i = 0; i < array_length; i++)
                temp[i] = (i < elem_num) ? array_pointer[i] : 0;
            delete[] array_pointer;
            array_pointer = temp;
        }
        if (val < current_min)
        {
            current_min = val;
        }
        array_pointer[elem_num] = val;
        elem_num++;
    }

    int size() const
    {
        return elem_num;
    }

    T min()
    {
        assert(elem_num != 0);
        return current_min;
    }

    T top()
    {
        assert(elem_num != 0);
        return array_pointer[elem_num - 1];
    }

    void pop()
    {
        assert(elem_num != 0);
        bool reassign_min = false; // If the minimum element in the array is to be deleted,
                                   // minimum element will be found afterwards.

        if (array_pointer[elem_num - 1] == current_min)
        {
            reassign_min = true;
        };

        array_pointer[elem_num - 1] = 0;
        --elem_num;

        if (reassign_min)
        {
            current_min = array_pointer[0];
            for (int i = 1; i < elem_num; i++)
            {
                if (array_pointer[i] < current_min)
                    current_min = array_pointer[i];
            }
        }
    }
};
