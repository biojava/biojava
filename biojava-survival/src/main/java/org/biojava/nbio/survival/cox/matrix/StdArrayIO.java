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

/*************************************************************************
 *  Compilation:  javac StdArrayIO.java
 *  Execution:    java StdArrayIO < input.txt
 *
 *  A library for reading in 1D and 2D arrays of integers, doubles,
 *  and booleans from standard input and printing them out to
 *  standard output.
 *
 *  % more tinyDouble1D.txt
 *  4
 *    .000  .246  .222  -.032
 *
 *  % more tinyDouble2D.txt
 *  4 3
 *    .000  .270  .000
 *    .246  .224 -.036
 *    .222  .176  .0893
 *   -.032  .739  .270
 *
 *  % more tinyBoolean2D.txt
 *  4 3
 *    1 1 0
 *    0 0 0
 *    0 1 1
 *    1 1 1
 *
 *  % cat tinyDouble1D.txt tinyDouble2D.txt tinyBoolean2D.txt | java StdArrayIO
 *  4
 *    0.00000   0.24600   0.22200  -0.03200
 *
 *  4 3
 *    0.00000   0.27000   0.00000
 *    0.24600   0.22400  -0.03600
 *    0.22200   0.17600   0.08930
 *    0.03200   0.73900   0.27000
 *
 *  4 3
 *  1 1 0
 *  0 0 0
 *  0 1 1
 *  1 1 1
 *
 *************************************************************************/


/**
 *  <i>Standard array IO</i>. This class provides methods for reading
 *  in 1D and 2D arrays from standard input and printing out to
 *  standard output.
 *  <p>
 *  For additional documentation, see
 *  <a href="http://introcs.cs.princeton.edu/22libary">Section 2.2</a> of
 *  <i>Introduction to Programming in Java: An Interdisciplinary Approach</i>
 *  by Robert Sedgewick and Kevin Wayne.
 */
public class StdArrayIO {



	/**
	 * Print an array of doubles to standard output.
	 * @param a
	 */
	public static void print(double[] a) {
		int N = a.length;
		System.out.println(N);
		for (int i = 0; i < N; i++) {
		 //   System.out.printf("%9.5f ", a[i]);
			System.out.print(a[i] + " ");
		}
		System.out.println();
	}




	/**
	 * Print the M-by-N array of doubles to standard output.
	 * @param a
	 */
	public static void print(double[][] a) {
		int M = a.length;
		int N = a[0].length;
		System.out.println(M + "x" + N);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
	//            System.out.printf("%9.5f ", a[i][j]);
				System.out.print(a[i][j] + " ");
			}
			System.out.println();
		}
	}




	/**
	 * Print an array of ints to standard output.
	 * @param a
	 */
	public static void print(int[] a) {
		int N = a.length;
		System.out.println(N);
		for (int i = 0; i < N; i++) {
			System.out.printf("%9d ", a[i]);
		}
		System.out.println();
	}




	/**
	 * Print the M-by-N array of ints to standard output.
	 * @param a
	 */
	public static void print(int[][] a) {
		int M = a.length;
		int N = a[0].length;
		System.out.println(M + " " + N);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				System.out.printf("%9d ", a[i][j]);
			}
			System.out.println();
		}
	}




	/**
	 * Print an array of booleans to standard output.
	 * @param a
	 */
	public static void print(boolean[] a) {
		int N = a.length;
		System.out.println(N);
		for (int i = 0; i < N; i++) {
			if (a[i]) System.out.print("1 ");
			else      System.out.print("0 ");
		}
		System.out.println();
	}



	/**
	 * Print the  M-by-N array of booleans to standard output.
	 * @param a
	 */
	public static void print(boolean[][] a) {
		int M = a.length;
		int N = a[0].length;
		System.out.println(M + " " + N);
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				if (a[i][j]) System.out.print("1 ");
				else         System.out.print("0 ");
			}
			System.out.println();
		}
	}


	/**
	 * Test client.
	 * @param args
	 */
	public static void main(String[] args) {


	}

}
