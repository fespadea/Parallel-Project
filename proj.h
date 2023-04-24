#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double normFro2(double** A, int n, int m){
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            sum += A[i][j] * A[i][j];
        }
    }
    return sum;
}

double norm1(double** A, int n, int m){
    double sumMax = 0;
    for(int j = 0; j < m; j++){
        double sum = 0;
        for(int i = 0; i < n; i++){
            sum += fabs(A[i][j]);
        }
        if(sum > sumMax){
            sumMax = sum;
        }
    }
    return sumMax;
}

/* function for exchanging two rows of
   a matrix */
void swap(double ** mat, int R, int C, int row1, int row2,
          int col)
{
    for (int i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}

/* function for finding rank of matrix */
int rankOfMatrix(double ** mat, int R, int C)
{
    int rank = C;
 
    for (int row = 0; row < rank; row++)
    {
        // Before we visit current row 'row', we make
        // sure that mat[row][0],....mat[row][row-1]
        // are 0.
 
        // Diagonal element is not zero
        if (mat[row][row])
        {
           for (int col = 0; col < R; col++)
           {
               if (col != row)
               {
                 // This makes all entries of current
                 // column as 0 except entry 'mat[row][row]'
                 double mult = (double)mat[col][row] /
                                       mat[row][row];
                 for (int i = 0; i < rank; i++)
                   mat[col][i] -= mult * mat[row][i];
              }
           }
        }
 
        // Diagonal element is already zero. Two cases
        // arise:
        // 1) If there is a row below it with non-zero
        //    entry, then swap this row with that row
        //    and process that row
        // 2) If all elements in current column below
        //    mat[r][row] are 0, then remove this column
        //    by swapping it with last column and
        //    reducing number of columns by 1.
        else
        {
            int reduce = 1;
 
            /* Find the non-zero element in current
                column  */
            for (int i = row + 1; i < R;  i++)
            {
                // Swap the row with non-zero element
                // with this row.
                if (mat[i][row])
                {
                    swap(mat, R, C, row, i, rank);
                    reduce = 0;
                    break;
                }
            }
 
            // If we did not find any row with non-zero
            // element in current column, then all
            // values in this column are 0.
            if (reduce)
            {
                // Reduce number of columns
                rank--;
 
                // Copy the last column here
                for (int i = 0; i < R; i ++)
                    mat[i][row] = mat[i][rank];
            }
 
            // Process this row again
            row--;
        }
 
       // Uncomment these lines to see intermediate results
       // display(mat, R, C);
       // printf("\n");
    }
    return rank;
}

double ** matrixSparsification(double ** A, int n, int m, double epsilon, double delta, int sMult, double alpha){
    // get paramters ready
    double AF2 = normFro2(A, n, m);
    double A1 = norm1(A, n, m);
    int k = rankOfMatrix(A, n, m);
    int s = sMult * k * (n+m);

    printf("%i\n", s);

    // calculate probabilities for each value of being chosen
    int totalLength = n*m;
    double * probabilities = (double *)malloc(sizeof(double)*(totalLength));
    double sum = 0;
    for(int i = 0; i < totalLength; i++){
        double Aij = A[(int)(i / m)][i % m];
        probabilities[i] = alpha * fabs(Aij) / A1 + (1 - alpha) * (Aij * Aij) / AF2;
        sum += probabilities[i];
    }
    printf("probabilities calculated\n");

    // choose the indexes using the probabilities
    int * choices = (int *)malloc(sizeof(int)*s);
    for(int i = 0; i < s; i++){
        double prob = ((double)rand() / (double)RAND_MAX) * sum;
        double probSum = 0;
        for(int j = 0; probSum < prob; j++){
            choices[i] = j;
            probSum += probabilities[j];
        }
        printf("%i\n", i);
    }
    printf("choices chosen\n");

    // combine the chosen values into a sparse matrix
    double ** ATilde = (double**)calloc(n * sizeof(double*), sizeof(double*));
    for(int i = 0; i < n; i++){
        ATilde[i] = (double*)calloc(m * sizeof(double), sizeof(double*));
    }
    for(int k = 0; k < s; k++){
        int choice = choices[k];
        int i = (int)(choice / m);
        int j = choice % m;
        ATilde[i][j] += A[i][j] / probabilities[choice] / s;
    }
    printf("sparsified\n");
    
    free(probabilities);
    free(choices);

    return ATilde;
}

double error(double ** A, double ** ATilde, int n, int m){
    double ** ADiff = (double**)malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++){
        ADiff[i] = (double*)malloc(m * sizeof(double));
        for(int j = 0; j < m; j++){
            ADiff[i][j] = A[i][j] - ATilde[i][j];
        }
    }
    double error = norm1(ADiff, n, m) / norm1(A, n, m);
    free(ADiff);
    return log(error) / log(2);
}