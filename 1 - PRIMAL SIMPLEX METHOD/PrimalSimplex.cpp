/* Program to Solve a Linear Programming Problem using the Primal Simplex Method */

#include<iostream>
#include<vector>

#include "SimplexSolver.h"

int main(){
    // Variable to store the number of Decision Variables in the LPP
    int numVar;
    // Variable to store the number of Constraints in the LPP
    int numConstraints;

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                                                                          |" << std::endl;
    std::cout << "|                                      PRIMAL SIMPLEX METHOD                               |" << std::endl;
    std::cout << "|                                                                                          |" << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;

    std::cout << "| ENTER THE NUMBER OF DECISION VARIABLES IN THE LPP : ";
    std::cin >> numVar;

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;    

    std::cout << "| ENTER THE NUMBER OF CONTRAINTS IN THE LPP : ";
    std::cin >> numConstraints;

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;

    // We need to keep track of Basic and Non Basic Variables
    // these are indicator vectors for the Basic and Non Basic Variables
    std::vector<int> BasicVariables(numConstraints);
    std::vector<int> NonBasicVariables(numVar);

    // SLACK VARIABLES :- 0, 1, 2, ..., numConstraints-1
    for(int i=0 ; i<numConstraints ; i++){
        BasicVariables[i] = i;
    }

    // ORIGINAL VARIABLES :- numConstraints+0, numConstraints+1, ..., numConstraints+numVar-1
    for(int i=numConstraints ; i<numConstraints+numVar ; i++){
        NonBasicVariables[i-numConstraints] = i;
    }

    // Vector to store the coefficients in the Objective Function
    std::vector<double> ObjectiveFunction(numVar+1);

    // 2D vector to store the Initial Simplex Table and then use it for Iterations
    std::vector< std::vector<double> > SimplexTable( numConstraints+1 , std::vector<double>(numVar+1) );

    std::cout << "| ENTER THE COEEFFICIENTS OF THE OBJECTIVE FUNCTION (ALONG WITH THE CONSTANT)              |" << std::endl;

    std::cout << "| -> ";
    
    for(int i=0 ; i<=numVar ; i++){
        std::cin >> ObjectiveFunction[i];
    }

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;

    std::cout << "| ENTER THE COEFFICIENTS IN THE LHS OF THE LPP CONSTRAINTS  ALONG WITH CONSTANTS(A|b)       |" << std::endl;

    std::cout << "+-------------------------------------------------------------------------------------------+" << std::endl;

    for(int i=0 ; i<numConstraints ; i++){
        std::cout << "| CONSTRAINT " << i+1 << " : ";
        for(int j=0 ; j<=numVar ; j++){
            std::cin >> SimplexTable[i][j];
        }
    }

    // Las Row of the Initial Simplex Table is the Objective Function Coefficients themselves with a '-' sign
    for(int i=0 ; i<=numVar ; i++){
        SimplexTable[numConstraints][i] = (-1.0)*(ObjectiveFunction[i]);
    }

    // Now we have the Initial Simplex Table and we can proceed for the Simplex Iterations

    // Printing the Initial Simplex Table before the Iterations start
    _LogSimplexTable_(SimplexTable, BasicVariables, NonBasicVariables);

    // First we need to find the Pivot Element in the Initial Simplex Table
    int PivotRow;
    int PivotColumn;

    _FindPivot_(SimplexTable, &PivotRow, &PivotColumn);

    // Now we need to swap the Entering and Departing Variable Values in their indiator vectors
    double bufferVar;

    bufferVar = BasicVariables[ PivotRow ];
    BasicVariables[ PivotRow ] = NonBasicVariables[ PivotColumn ];
    NonBasicVariables[ PivotColumn ] = bufferVar;

    // now we are ready to start the Primal Simplex Iterations
    _SimplexIterations_(SimplexTable, BasicVariables, NonBasicVariables, &PivotRow, &PivotColumn);

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                  FINAL OPTIMAL SOLUTION                                  |" << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| OPTIMAL VALUE FOR THE LPP IS : Z = " << SimplexTable[numConstraints][numVar] << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| THIS OPTIMAL SOLUTION IS ACHIEVED WHEN :-                                                |" << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    for(int i=0 ; i<BasicVariables.size() ; i++){
        if( BasicVariables[i] >= numConstraints ){
            std::cout << "| X" << BasicVariables[i]-numConstraints+1 << " = " << SimplexTable[i][numVar] << std::endl;
        }
    }
    for(int i=0 ; i<NonBasicVariables.size() ; i++){
        if(NonBasicVariables[i] >= numConstraints){
            std::cout << "| X" << NonBasicVariables[i] - numConstraints + 1 << " = " << 0 << std::endl;
        }
    }

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;

    // Also we need to check for alternate optimal solutions if any indicator variable is 0
    int checkAlternateOptimal = _CheckAlternateOptimal_(SimplexTable);

    if( checkAlternateOptimal == 1 ){
        // Alternate optimal solution exists
        _FindAlternateOptimal_(SimplexTable, BasicVariables, NonBasicVariables);
    
        _LogSimplexTable_(SimplexTable, BasicVariables, NonBasicVariables);
    }

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                              ALTERNATE OPTIMAL SOLUTION                                  |" << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| OPTIMAL VALUE FOR THE LPP IS : Z = " << SimplexTable[numConstraints][numVar] << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "| THIS ALTERNATE OPTIMAL SOLUTION IS ACHIEVED WHEN :-                                      |" << std::endl;
    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;
    for(int i=0 ; i<BasicVariables.size() ; i++){
        if( BasicVariables[i] >= numConstraints ){
            std::cout << "| X" << BasicVariables[i]-numConstraints+1 << " = " << SimplexTable[i][numVar] << std::endl;
        }
    }
    for(int i=0 ; i<NonBasicVariables.size() ; i++){
        if(NonBasicVariables[i] >= numConstraints){
            std::cout << "| X" << NonBasicVariables[i] - numConstraints + 1 << " = " << 0 << std::endl;
        }
    }

    std::cout << "+------------------------------------------------------------------------------------------+" << std::endl;


    return 0;
}