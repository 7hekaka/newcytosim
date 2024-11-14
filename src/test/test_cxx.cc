// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/*
 This is a test of C++ extensions
 defined by ISOC++11
 
 http://en.wikipedia.org/wiki/C%2B%2B11
 
 it should be compiled as: c++ -std=c++17
 
 C++ version:
 199711L
 201103L
 201300L
 201402L
 201703L
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <new>


class Bot
{
protected:
    int x;
public:
    Bot() { printf(" Bot\n"); }
    Bot(int) { x=0; printf("Bot(0)\n"); }
    ~Bot() { printf("~Bot\n"); }
    virtual int val() { return x; }
    virtual void func() { ++x; printf("Bot::inc\n"); }
};

class Top: public Bot
{
public:
    Top() { printf(" Top\n"); }
    ~Top() { printf("~Top\n"); }
    virtual void func() { --x; printf("Top::dec\n");  }
};


/*
 This test may fail with optimizations level that include de-virtualization:
 hence compile with: -fno-devirtualize
 */
void test_placement()
{
    printf("sizeof(Bot) = %lu ", sizeof(Bot));
    printf("sizeof(Top) = %lu\n", sizeof(Top));
    
    Bot B(0);
    (&B)->func();
    printf("value = %i\n", B.val());

    // construct Obj in place of B
    Top* O = new (&B) Top;
    
    if ( O != &B )
        printf("pointers: %p != %p\n", &B, O);

    B.func();      // not-virtual: Base::func called
    (&B)->func();  // virtual: Obj::func is called
    
    if ( O->val() == 1 )
        printf("PASSED: value == 1\n");
    else
        printf("FAILED: value %i != 1\n", O->val());

    // destructor must be called explicitly
    O->~Top();
}


void test_vector()
{
    std::vector<double> vec(3);
    vec[0] = 1;
    vec[1] = 2;
    vec[2] = 3;
    std::clog << "sizeof std::vector<>(3) = " << sizeof(vec) << '\n';
    std::clog << "  vec    @ " << &vec << '\n';
    std::clog << "  vec[0] @ " << &vec[0] << '\n';
    std::clog << "  vec[1] @ " << &vec[1] << '\n';
    std::clog << "  vec[2] @ " << &vec[2] << '\n';

    std::vector<double> vec6(6);
    vec6[0] = 1;
    vec6[1] = 2;
    vec6[2] = 3;
    std::clog << "sizeof std::vector<>(6) = " << sizeof(vec6) << '\n';
}

void test_array()
{
    int array[] = { 1, 2, 3, 5, 7, 11, 13, 17, 23 };
    for (int &x : array) {
        std::cout << " " << x;
    }
}

void test_constants()
{
    // hexadecimal floating point literal is a C++17 feature
    constexpr double CONSTANT = 0x1p-31;
    
    // constexpr initialization is a GNU extension:
    // constexpr double SQ3 = std::sqrt(3);
    const double SQ3 = std::sqrt(3);
    
    std::clog << "0x1p-31 = " << CONSTANT << '\n';
    std::clog << "std::sqrt(3) = " << SQ3 << '\n';
}

int main ()
{
    std::clog << "C++ version " << __cplusplus << '\n';
    std::clog << __VERSION__ << '\n';

    test_placement();
    std::cout << std::endl;
}
