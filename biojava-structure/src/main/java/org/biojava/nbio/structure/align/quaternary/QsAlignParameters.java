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

	/**
	 * The maximum allowed distance between two aligned residues of equivalent
	 * Subunits.
	 * 
	 * @return
	 */
	public double getdCutoff() {
		return dCutoff;
	}

	public void setdCutoff(double dCutoff) {
		this.dCutoff = dCutoff;
	}

}
