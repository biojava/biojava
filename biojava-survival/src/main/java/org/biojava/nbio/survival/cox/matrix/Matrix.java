/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */
package org.biojava.nbio.survival.cox.matrix;


/**
 *
 * http://introcs.cs.princeton.edu/java/22library/Matrix.java.html
 */
/**
 * ***********************************************************************
 * Compilation: javac Matrix.java Execution: java Matrix
 *
 * A bare-bones collection of static methods for manipulating matrices.
 *
 ************************************************************************
 */
public class Matrix {



	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[][] sqrt(double[][] A) {
		double[][] d = new double[A.length][A[0].length];
		for (int i = 0; i < d.length; i++) {
			for (int j = 0; j < d[0].length; j++) {
				d[i][j] = Math.sqrt(A[i][j]);
			}
		}
		return d;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[] sqrt(double[] A) {
		double[] d = new double[A.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = Math.sqrt(A[i]);
		}
		return d;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[][] oneDivide(double[][] A) {
		double[][] d = new double[A.length][A[0].length];
		for (int i = 0; i < d.length; i++) {
			for (int j = 0; j < d[0].length; j++) {
				d[i][j] = 1.0 / A[i][j];
			}
		}
		return d;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[] oneDivide(double[] A) {
		double[] d = new double[A.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = 1.0 / A[i];
		}
		return d;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[] diag(double[][] A) {
		double[] d = new double[A.length];
		for (int i = 0; i < d.length; i++) {
			d[i] = A[i][i];
		}
		return d;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[][] diag(double[] A) {
		double[][] d = new double[A.length][A.length];
		for (int i = 0; i < d.length; i++) {
			d[i][i] = A[i];
		}
		return d;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[][] abs(double[][] A) {
		int m = A.length;
		int n = A[0].length;
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = Math.abs(A[i][j]);
			}
		}
		return C;
	}

	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[] abs(double[] A) {
		int m = A.length;

		double[] C = new double[m];
		for (int i = 0; i < m; i++) {
			C[i] = Math.abs(A[i]);
		}
		return C;
	}

	// return a random m-by-n matrix with values between 0 and 1
	/**
	 *
	 * @param m
	 * @param n
	 * @return
	 */
	public static double[][] random(int m, int n) {
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = Math.random();
			}
		}
		return C;
	}

	// return n-by-n identity matrix I
	/**
	 *
	 * @param n
	 * @return
	 */
	public static double[][] identity(int n) {
		double[][] I = new double[n][n];
		for (int i = 0; i < n; i++) {
			I[i][i] = 1;
		}
		return I;
	}

	// return x^T y
	/**
	 *
	 * @param x
	 * @param y
	 * @return
	 */
	public static double dot(double[] x, double[] y) {
		if (x.length != y.length) {
			throw new RuntimeException("Illegal vector dimensions.");
		}
		double sum = 0.0;
		for (int i = 0; i < x.length; i++) {
			sum += x[i] * y[i];
		}
		return sum;
	}

	// return C = A^T
	/**
	 *
	 * @param A
	 * @return
	 */
	public static double[][] transpose(double[][] A) {
		int m = A.length;
		int n = A[0].length;
		double[][] C = new double[n][m];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[j][i] = A[i][j];
			}
		}
		return C;
	}

	// return C = A + B
	/**
	 *
	 * @param A
	 * @param B
	 * @return
	 */
	public static double[][] add(double[][] A, double[][] B) {
		int m = A.length;
		int n = A[0].length;
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] + B[i][j];
			}
		}
		return C;
	}

	// return C = A - B
	/**
	 *
	 * @param A
	 * @param B
	 * @return
	 */
	public static double[][] subtract(double[][] A, double[][] B) {
		int m = A.length;
		int n = A[0].length;
		double[][] C = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				C[i][j] = A[i][j] - B[i][j];
			}
		}
		return C;
	}

	// return C = A * B
	/**
	 *
	 * @param A
	 * @param B
	 * @return
	 */
	public static double[][] multiply(double[][] A, double[][] B) {
		int mA = A.length;
		int nA = A[0].length;
		int mB = B.length;
		int nB = B[0].length;
		if (nA != mB) {
			throw new RuntimeException("Illegal matrix dimensions.");
		}
		double[][] C = new double[mA][nB];
		for (int i = 0; i < mA; i++) {
			for (int j = 0; j < nB; j++) {
				for (int k = 0; k < nA; k++) {
					C[i][j] += (A[i][k] * B[k][j]);
				}
			}
		}
		return C;
	}

	// matrix-vector multiplication (y = A * x)
	/**
	 *
	 * @param A
	 * @param x
	 * @return
	 */
	public static double[] multiply(double[][] A, double[] x) {
		int m = A.length;
		int n = A[0].length;
		if (x.length != n) {
			throw new RuntimeException("Illegal matrix dimensions.");
		}
		double[] y = new double[m];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				y[i] += (A[i][j] * x[j]);
			}
		}
		return y;
	}

	/**
	 *
	 * @param A
	 * @param x
	 * @return
	 */
	public static double[][] scale(double[][] A, double[] x) {
		int m = A.length;
		int n = A[0].length;
		if (x.length != A.length) {
			throw new RuntimeException("Illegal matrix dimensions.");
		}
		double[][] y = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				y[i][j] = A[i][j] * x[i];
			}
		}
		return y;
	}

	/**
	 *
	 * @param A
	 * @param x
	 * @return
	 */
	public static double[][] scale(double[][] A, double x) {
		int m = A.length;
		int n = A[0].length;

		double[][] y = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				y[i][j] = A[i][j] * x;
			}
		}
		return y;
	}


	// vector-matrix multiplication (y = x^T A)
	/**
	 *
	 * @param x
	 * @param A
	 * @return
	 */
	public static double[] multiply(double[] x, double[][] A) {
		int m = A.length;
		int n = A[0].length;
		if (x.length != m) {
			throw new RuntimeException("Illegal matrix dimensions.");
		}
		double[] y = new double[n];
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				y[j] += (A[i][j] * x[i]);
			}
		}
		return y;
	}

	// test client
	/**
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("D");
		System.out.println("--------------------");
		double[][] d = {{1, 2, 3}, {4, 5, 6}, {9, 1, 3}};
		StdArrayIO.print(d);
		System.out.println();

		System.out.println("I");
		System.out.println("--------------------");
		double[][] c = Matrix.identity(5);
		StdArrayIO.print(c);
		System.out.println();

		System.out.println("A");
		System.out.println("--------------------");
		double[][] a = Matrix.random(5, 5);
		StdArrayIO.print(a);
		System.out.println();

		System.out.println("A^T");
		System.out.println("--------------------");
		double[][] b = Matrix.transpose(a);
		StdArrayIO.print(b);
		System.out.println();

		System.out.println("A + A^T");
		System.out.println("--------------------");
		double[][] e = Matrix.add(a, b);
		StdArrayIO.print(e);
		System.out.println();

		System.out.println("A * A^T");
		System.out.println("--------------------");
		double[][] f = Matrix.multiply(a, b);
		StdArrayIO.print(f);
		System.out.println();
	}
}
