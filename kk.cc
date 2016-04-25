#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;

int main ( int argc, char* argv[] )
{
    if ( argc != 2 ) // argc should be 2 for correct execution
    {
        // We print argv[0] assuming it is the program name
        cout <<"usage: ./kk  <inputfile>" << endl;

        // End program
        return -1;
    }

    // Call python file
    std::string input = argv[1];

    std::string command = "python assignment.py " + input;
    system(command.c_str());
}