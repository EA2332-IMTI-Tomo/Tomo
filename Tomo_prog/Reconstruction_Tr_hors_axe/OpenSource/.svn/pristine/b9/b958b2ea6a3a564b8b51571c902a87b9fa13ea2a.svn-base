#include "vChrono.h"
//#include <boost.h> testé sur boost 1.4.8
#include <iostream>


using namespace std;


int main()
{
  vChrono<boost::chrono::system_clock> t1;
  vChrono<boost::chrono::steady_clock> t2;
  vChrono<boost::chrono::high_resolution_clock> t3;

  std::cout << "Type the Enter key: ";
  std::cin.get();

  std::cout << std::fixed; //<< std::setprecision(9);
  std::cout << "system_clock-----------: "
            << t1.seconds() << " seconds\n";
  std::cout << "steady_clock--------: "
            << t2.seconds() << " seconds\n";
  std::cout << "high_resolution_clock--: "
            << t3.seconds() << " seconds\n and" << t3.milliseconds() << " ms\n" ;


  
  t1.reset();
  std::cout << "resetted. Type the Enter key: ";
  std::cin.get();
  
  std::cout << "system_clock-----------: "
            << t1.seconds() << " seconds\n";


  return 0;
}



// g++ vChrono.cc  -o clock -L /usr/local/phd/boost/lib/ -lboost_system -lboost_chrono
