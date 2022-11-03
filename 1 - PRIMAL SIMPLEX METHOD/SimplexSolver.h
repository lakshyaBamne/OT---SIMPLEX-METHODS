/* HEADER FILE FOR METHODS USED IN PRIMAL SIMPLEX ITERATIONS */

#include<iostream>
#include<vector>

// Funtion to Print a Simplex Table to the Console
void _LogSimplexTable_(std::vector< std::vector<double> > SimplexTable, std::vector<int> BasicVariables, std::vector<int> NonBasicVariables){
    // Temporary variables for the number of rows and columns in the Simplex Table
    int numConstraints = BasicVariables.size();
    int numvar = NonBasicVariables.size();

    std::cout << "+-------------------------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                     SIMPLEX TABLE                                         |" << std::endl;
    std::cout << "+-------------------------------------------------------------------------------------------+" << std::endl;

    // First we need to Print the Labels for the COLUMNS
    for(int i=0 ; i<numvar ; i++){
        std::cout << "|     ";

        if( NonBasicVariables[i] < numConstraints ){
            // SLACK
            std::cout << "S" << NonBasicVariables[i] + 1 << "     ";
        }
        else{
            // ORIGINAL 
            std::cout << "X" << NonBasicVariables[i] - numConstraints + 1 << "     ";
        }
    }
    std::cout << "|      1      |" << std::endl;

    // Now we can print the values in the Simplex Table Row wise
    for( int i=0 ; i<numConstraints ; i++ ){
        // Actual values in the Rows of the Simplex Table
        for(int j=0 ; j<=numvar ; j++){
            std::cout << "|      " << SimplexTable[i][j] << "     ";
        }

        // at the end of each row we have to print the corresponding Basic Variable
        if( BasicVariables[i] < numConstraints ){
            // SLACK VARIABLE
            std::cout << "|     =S" << BasicVariables[i] + 1;
        }
        else{
            // ORIGINAL VARIABLE
            std::cout << "|     =X" << BasicVariables[i] - numConstraints + 1;
        }

        std::cout << "\n";
    }
    std::cout << "+-------------------------------------------------------------------------------------------+" << std::endl;
    for(int i=0 ; i<=numvar ; i++){
        std::cout << "|     " << SimplexTable[numConstraints][i] << "     ";
    }
    std::cout << "|     =Z" << std::endl;
    std::cout << "+-------------------------------------------------------------------------------------------+" << std::endl;

}

// Function used to find the Pivot Element in a given simplex table
void _FindPivot_(std::vector< std::vector<double> > SimplexTable, int* PivotRow, int* PivotColumn){
    int IndicatorRowIndex = SimplexTable.size()-1;
    int numColumns = SimplexTable[0].size();

    double minIndicator = 0;

    int ColumnIndex;
    int RowIndex;

    // Loop to find the Column in which the Pivot element is present
    for(int i=0 ; i<numColumns-1 ; i++){
        if( SimplexTable[IndicatorRowIndex][i] < minIndicator ){
            minIndicator = SimplexTable[IndicatorRowIndex][i];
            ColumnIndex = i;
        }
    }

    double MinRatio = 1000000000.0;

    double Numerator, Denominator;
    double IndividualRatio; 

    for(int i=0 ; i<SimplexTable.size()-1 ; i++){
        Numerator = SimplexTable[ i ][ numColumns-1 ];
        Denominator = SimplexTable[ i ][ ColumnIndex ];

        IndividualRatio = Numerator / Denominator;

        // The LEAST POSITIVE ratio must only be selected
        if( (IndividualRatio < MinRatio) && ( IndividualRatio > 0 ) ){
            MinRatio = IndividualRatio;
            RowIndex = i;
        }
    }

    std::cout << "Pivot Row is    : " << RowIndex+1 << std::endl; 
    std::cout << "Pivot Column is : " << ColumnIndex+1 << std::endl;
    
    *PivotRow = RowIndex;
    *PivotColumn = ColumnIndex;

}

// Function to carry out Simplex Iterations after Initial Table has been formulated
void _SimplexIterations_(std::vector< std::vector<double> > &SimplexTable, std::vector<int> &BasicVariables, std::vector<int> &NonBasicVariables, int* PivotRow, int* PivotColumn){
    // We need an indicator variable to determine when to stop the Simplex Iterations
    int IndicatorCount = 1;

    double PivotValue = SimplexTable[*PivotRow][*PivotColumn];

    int numRow = SimplexTable.size();
    int numColumn = SimplexTable[0].size();

    std::vector< std::vector<double> > SecondaryTable( numRow, std::vector<double>(numColumn) );

    double bufferVar;

    double IndividualValue, Numerator, Denominator;

    while( IndicatorCount != 0 ){
        // Initializing the indicator variable for this particular iteration        
        IndicatorCount = 0;

        // First we need to calculate new simplex table
        SecondaryTable[*PivotRow][*PivotColumn] = 1.0 / SimplexTable[*PivotRow][*PivotColumn];

        for( int i=0 ; i<numColumn ; i++ ){
            if( i != *PivotColumn ){
                SecondaryTable[*PivotRow][i] = (SimplexTable[*PivotRow][i]) / (SimplexTable[*PivotRow][*PivotColumn]);
            }
        }

        for( int i=0 ; i<numRow ; i++ ){
            if( i != *PivotRow ){
                SecondaryTable[i][*PivotColumn] = (-1.0) * ( (SimplexTable[i][*PivotColumn]) / (SimplexTable[*PivotRow][*PivotColumn]) );
            }
        }

        for(int i=0 ; i<numRow ; i++){
            for(int j=0 ; j<numColumn ; j++){
                if( (i != *PivotRow) && (j != *PivotColumn) ){
                    Numerator = ( SimplexTable[*PivotRow][*PivotColumn] * SimplexTable[i][j] ) - ( SimplexTable[i][*PivotColumn] * SimplexTable[*PivotRow][j] );
                    Denominator = SimplexTable[*PivotRow][*PivotColumn];
                
                    IndividualValue = Numerator / Denominator;

                    SecondaryTable[i][j] = IndividualValue;
                }
            }
        }

        // Now we need to copy the elements from secondary Simplex Table to Primary Simplex Table
        for(int i=0 ; i<numRow ; i++){
            for(int j=0 ; j<numColumn ; j++){
                SimplexTable[i][j] = SecondaryTable[i][j];
            }
        }

        // We can print the new simplex table after recalculation
        _LogSimplexTable_(SimplexTable, BasicVariables, NonBasicVariables);

        // Now we need to check if there are any negative elements in the Indicator row
        for(int i=0 ; i<numColumn-1 ; i++){
            if( SimplexTable[numRow-1][i] < 0 ){
                IndicatorCount++;
            }
        }

        if(IndicatorCount == 0){
            break;
        }

        // If loop reaches here, it means there are negative elements in the Indicator variable
        // So we need to find the new pivot for the next iteration
        _FindPivot_(SimplexTable, PivotRow, PivotColumn);

        // now we have the pivot so we need to swap the Indicators for the Entering and Departing Variables
        bufferVar = BasicVariables[*PivotRow];
        BasicVariables[*PivotRow] = NonBasicVariables[*PivotColumn];
        NonBasicVariables[*PivotColumn] = bufferVar;
    }
}

// Function to check if the LPP has any alternate optimal solution
// Returns 1 if alternate optimal exists otherwise 0
int _CheckAlternateOptimal_(std::vector< std::vector<double> > SimplexTable){
    int numRow = SimplexTable.size();
    int numColumn = SimplexTable[0].size();
    
    // variable to give solution back
    int output = 0;

    // we need to find if a zero exists in the Indicator Row
    int countZero = 0;

    for(int i=0 ; i<numColumn-1 ; i++){
        if( SimplexTable[numRow-1][i] == 0 ){
            // zero found
            output = 1;
        }
    }

    return output;

}

// Function to find alternate optimal solution for an LPP given it has an alternate optimal solution
void _FindAlternateOptimal_(std::vector< std::vector<double> > &SimplexTable, std::vector<int> &BasicVariables, std::vector<int> &NonBasicVariables){
    // First we need to find the pivot corresponding to the Alternate optimal solution
    int PivotRow;
    int PivotColumn;

    // variables to store the number of rows and columns in the Simplex Table
    int numRow = SimplexTable.size();
    int numColumn = SimplexTable[0].size();

    for( int i=0 ; i<numColumn-1 ; i++ ){
        if( SimplexTable[numRow-1][i] == 0 ){
            PivotColumn = i;
            break;
        }
    }

    // now we need to find the minimum positive ration to find the Pivot
    double IndividualRatio;

    double minRatio = 1000000000.0;

    for(int i=0 ; i<numRow-1 ; i++){
        IndividualRatio = SimplexTable[i][numColumn-1] / SimplexTable[i][PivotColumn];

        if((IndividualRatio < minRatio) && (IndividualRatio > 0)){
            minRatio = IndividualRatio;
            PivotRow = i;
        }
    }

    std::cout << "Pivot Position is : ( " << PivotRow << " , " << PivotColumn << " )" << std::endl;

    // now we have the position of the pivot and we can calculate the simplex table
    std::vector< std::vector<double> > SecondaryTable( numRow, std::vector<double>(numColumn) );

    SecondaryTable[PivotRow][PivotColumn] =  1.0 / SimplexTable[PivotRow][PivotColumn];

    for(int i=0 ; i<numColumn ; i++){
        if( i != PivotColumn ){
            SecondaryTable[PivotRow][i] = (SimplexTable[PivotRow][i]) / (SimplexTable[PivotRow][PivotColumn]);
        }
    }

    for(int i=0 ; i<numRow ; i++){
        if( i != PivotRow ){
            SecondaryTable[i][PivotColumn] = (-1)*( (SimplexTable[i][PivotColumn]) / (SimplexTable[PivotRow][PivotColumn]) );
        }
    }

    double IndividualValue, Numerator, Denominator;

    for(int i=0 ; i<numRow ; i++){
        for(int j=0 ; j<numColumn ; j++){
            if( (i != PivotRow) && (j != PivotColumn ) ){
                Numerator = ( SimplexTable[PivotRow][PivotColumn] * SimplexTable[i][j] ) - ( SimplexTable[i][PivotColumn] * SimplexTable[PivotRow][j] );
                Denominator = SimplexTable[PivotRow][PivotColumn];
            
                IndividualValue = Numerator / Denominator;

                SecondaryTable[i][j] = IndividualValue;
            }    
        }
    }

    // Also we need to swap the Indicator variables for the Entering and Departing variables
    double bufferVar;

    bufferVar = BasicVariables[PivotColumn];
    BasicVariables[PivotColumn] = NonBasicVariables[PivotRow];
    NonBasicVariables[PivotRow] = bufferVar;

    // we need to change the Primary vector which has been passed as a Reference
    for(int i=0 ; i<numRow ; i++){
        for(int j=0 ; j<numColumn ; j++){
            SimplexTable[i][j] = SecondaryTable[i][j];
        }
    }

}

