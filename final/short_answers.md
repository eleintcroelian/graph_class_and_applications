1.1
printed value from line 2:
1
printed value from line 3:
2

1.2
First difference between a struct and a class is that members of struct
are public by default and members of class can be set to public, private
or protected.
Second difference is that classes deallocate memory before being destructed
but we cannot assume structs will do the same.
reference: https://www.learncpp.com/cpp-tutorial/82-classes-and-class-members/

1.3
We need to create a default constructor, otherwise it is a compile error. Compilers don't create default constructors automatically.

1.4
For example initializing an integer variable with an double value is valid
in initializing with "=" (although it may raise a warning).

However the same way of initialization is not valid for brackets and it
causes an error by the compiler. Initialization with brackets is safer because
it prevents data loss. In addition, "{}" leads to performance gains over "="
therefore, it is preferable.

See the code snippet below:
"a" holds the value "1" as a result of "=" initialization, decimal part is
discarded.
"{}" initialization throws an error, therefore saves us from data loss. 
	
int main()
{
  int a=1.1;
  int b {1.1};
  
  return 0;
}

1.5
The memory allocated on the stack is deallocated after it is popped
off the stack. Therefore it may not have the value we expect and it is
not safe to use it.
https://www.learncpp.com/cpp-tutorial/79-the-stack-and-the-heap/

1.6

RAII keeps memory management with constructors and destructors. When we hide a new in a class constructor, and delete in a destructor, we automatically call new and delete as we get in and out of a scope using this proxy classes.

1.7
Lifetime of an object is the time between its initialization and the call
of destructor (deallocation of the memory). We need to consider if an object is valid for example when we want to access data attributes or if it is a pointer, when we are dereferencing it.

1.8
The compilation is held in two phases, where it first allocates the templates and then when the method being called, filling in the template. Thus, we need to place template definitions inside the header files.