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
