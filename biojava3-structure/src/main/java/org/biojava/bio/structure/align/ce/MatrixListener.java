package org.biojava.bio.structure.align.ce;



public interface MatrixListener {

	public double[][] matrixInOptimizer(double[][] max);

	public boolean[][] initializeBreakFlag(boolean[][] brkFlag);
}
