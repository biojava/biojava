package org.biojava.nbio.structure.align.quaternary;

/**
 * The parameter bean for the {@link QsAlign} algorithm.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class QsAlignParameters {

	private double dCutoff = 10.0;
	private double maxRmsd = 7.0;
	private double minOrientationMetric = Math.PI / 8; // 45 degree

	/**
	 * The maximum allowed distance between the centroids of two equivalent
	 * Subunits, in A.
	 * 
	 * @return dCutoff
	 */
	public double getdCutoff() {
		return dCutoff;
	}

	/**
	 * The maximum allowed distance between the centroids of two equivalent
	 * Subunits, in A.
	 * 
	 * @param dCutoff
	 */
	public void setdCutoff(double dCutoff) {
		this.dCutoff = dCutoff;
	}

	
	public double getMaxRmsd() {
		return maxRmsd;
	}

	public void setMaxRmsd(double maxRmsd) {
		this.maxRmsd = maxRmsd;
	}

	public double getMinOrientationMetric() {
		return minOrientationMetric;
	}

	public void setMinOrientationMetric(double minOrientationMetric) {
		this.minOrientationMetric = minOrientationMetric;
	}

}
