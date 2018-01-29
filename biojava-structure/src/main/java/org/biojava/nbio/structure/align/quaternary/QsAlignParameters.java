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
	private double maxRmsd = 10.0;
	private double maxOrientationAngle = Math.PI / 6; // 30 degree

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

	/**
	 * The maximum allowed RMSD of the alignment, in A.
	 * 
	 * @return maxRmsd
	 */
	public double getMaxRmsd() {
		return maxRmsd;
	}

	/**
	 * The maximum allowed RMSD of the alignment, in A.
	 * 
	 * @param maxRmsd
	 */
	public void setMaxRmsd(double maxRmsd) {
		this.maxRmsd = maxRmsd;
	}

	/**
	 * The maximum orientation angle between two equivalent Subunits, in
	 * radians. Range [0, Pi].
	 * 
	 * @return the maximum orientation angle
	 */
	public double getMaxOrientationAngle() {
		return maxOrientationAngle;
	}

	/**
	 * The maximum orientation angle between two equivalent Subunits, in
	 * radians. Range [0, Pi].
	 * 
	 * @param maxOrientationAngle
	 *            maximum orientation angle
	 */
	public void setMaxOrientationAngle(double maxOrientationAngle) {
		this.maxOrientationAngle = maxOrientationAngle;
	}

}
