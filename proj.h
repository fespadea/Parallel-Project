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

// Code from geeksforgeeks: https://www.geeksforgeeks.org/program-for-rank-of-matrix/
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

// Code from geeksforgeeks: https://www.geeksforgeeks.org/program-for-rank-of-matrix/
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

// Code from geeksforgeeks: https://www.geeksforgeeks.org/merge-sort/
void merge(double arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
 
    /* create temp arrays */
    double * L = (double *)malloc(sizeof(double) * n1);
    double * R = (double *)malloc(sizeof(double) * n2);
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
 
// Code from geeksforgeeks: https://www.geeksforgeeks.org/merge-sort/
void mergeSort(double arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
 
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
 
        merge(arr, l, m, r);
    }
}

double ** matrixSparsification(double ** A, int n, int m, double epsilon, double delta, int sMult, double alpha){
    // get paramters ready
    double AF2 = normFro2(A, n, m);
    double A1 = norm1(A, n, m);
    int k = rankOfMatrix(A, n, m);
    int s = sMult * k * (n+m);
    printf("s: %i\n", s);

    // calculate probabilities for each value of being chosen
    int totalLength = n*m;
    double * probabilities = (double *)malloc(sizeof(double)*(totalLength));
    double sum = 0;
    for(int i = 0; i < totalLength; i++){
        double Aij = A[(int)(i / m)][i % m];
        probabilities[i] = alpha * fabs(Aij) / A1 + (1 - alpha) * (Aij * Aij) / AF2;
        sum += probabilities[i];
    }

    // choose the indexes using the probabilities
    double * probs = (double *)malloc(sizeof(double)*s);
    for(int i = 0; i < s; i++){
        int ra = rand();
        probs[i] = ((double)ra / (double)RAND_MAX) * sum;
    }
    mergeSort(probs, 0, s-1);
    int * choices = (int *)malloc(sizeof(int)*s);
    double probSum = 0;
    int p = 0;
    for(int j = 0; p < s; j++){
        probSum += probabilities[j];
        while(probSum >= probs[p] && p < s){
            choices[p] = j;
            p++;
        }
    }

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
    
    free(probabilities);
    free(choices);
    free(probs);

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