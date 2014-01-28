/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.survival.cox.stats;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class Cholesky2 {
   
       /* $Id: cholesky2.c 11357 2009-09-04 15:22:46Z therneau $
     **
     ** subroutine to do Cholesky decompostion on a matrix: C = FDF'
     **   where F is lower triangular with 1's on the diagonal, and D is diagonal
     **
     ** arguments are:
     **     n         the size of the matrix to be factored
     **     **matrix  a ragged array containing an n by n submatrix to be factored
     **     toler     the threshold value for detecting "singularity"
     **
     **  The factorization is returned in the lower triangle, D occupies the
     **    diagonal and the upper triangle is left undisturbed.
     **    The lower triangle need not be filled in at the start.
     **
     **  Return value:  the rank of the matrix (non-negative definite), or -rank
     **     it not SPD or NND
     **
     **  If a column is deemed to be redundant, then that diagonal is set to zero.
     **
     **   Terry Therneau
     */
    /**
     *
     * @param matrix
     * @param n
     * @param toler
     * @return
     */
    public static int process(double[][] matrix, int n, double toler) {
        double temp;
        int i, j, k;
        double eps, pivot;
        int rank;
        int nonneg;

        nonneg = 1;
        eps = 0;
        for (i = 0; i < n; i++) {
            if (matrix[i][i] > eps) {
                eps = matrix[i][i];
            }
            for (j = (i + 1); j < n; j++) {
                matrix[j][i] = matrix[i][j];
            }
        }
        eps *= toler;

        rank = 0;
        for (i = 0; i < n; i++) {
            pivot = matrix[i][i];
            if (pivot < eps) {
                matrix[i][i] = 0;
                if (pivot < -8 * eps) {
                    nonneg = -1;
                }
            } else {
                rank++;
                for (j = (i + 1); j < n; j++) {
                    temp = matrix[j][i] / pivot;
                    matrix[j][i] = temp;
                    matrix[j][j] -= temp * temp * pivot;
                    for (k = (j + 1); k < n; k++) {
                        matrix[k][j] -= temp * matrix[k][i];
                    }
                }
            }
        }
        return (rank * nonneg);
    }
    
}
