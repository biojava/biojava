package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;

/**
 * This class stores all the relevant information of an internal symmetry result
 * obtained with CeSymm. The purpose is to carry all the information packed
 * during the calculations and return a single object.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmResult {

	private MultipleAlignment multipleAlignment;
	private AFPChain selfAlignment;
	private Atom[] atoms;

	private CESymmParameters params = new CESymmParameters();
	private SymmetryAxes axes;
	private QuatSymmetryResults PG;

	private int symmOrder = 1;
	private boolean refined;
	private SymmetryType type;

	/**
	 * Return true if the symmetry result is significant, false
	 * otherwise.
	 */
	public boolean isSignificant() {

		// In any case if the order is 1 it is asymmetric
		if (symmOrder < 1)
			return false;

		// If the old version was turned ON
		if (params.getRefineMethod() == RefineMethod.NOT_REFINED) {
			if (selfAlignment.getTMScore() < params.getScoreThreshold())
				return false;
			else
				return true;
		}

		// For the new version
		if (refined) {
			// Condition is over min size and score threshold
			if (multipleAlignment.getScore(MultipleAlignmentScorer.AVGTM_SCORE) > params
					.getScoreThreshold()
					&& multipleAlignment.getCoreLength() > params
							.getMinCoreLength())
				return true;
			else
				return false;
		} else
			return false;
	}

	public MultipleAlignment getMultipleAlignment() {
		return multipleAlignment;
	}

	public void setMultipleAlignment(MultipleAlignment multipleAlignment) {
		this.multipleAlignment = multipleAlignment;
	}

	public AFPChain getSelfAlignment() {
		return selfAlignment;
	}

	public void setSelfAlignment(AFPChain selfAlignment) {
		this.selfAlignment = selfAlignment;
	}

	public CESymmParameters getParams() {
		return params;
	}

	public void setParams(CESymmParameters params) {
		this.params = params.clone();
	}

	public SymmetryAxes getAxes() {
		return axes;
	}

	public void setAxes(SymmetryAxes axes) {
		this.axes = axes;
	}

	public int getSymmOrder() {
		return symmOrder;
	}

	public void setSymmOrder(int symmOrder) {
		this.symmOrder = symmOrder;
	}

	public boolean isRefined() {
		return refined;
	}

	public void setRefined(boolean refined) {
		this.refined = refined;
	}

	public QuatSymmetryResults getPG() {
		return PG;
	}

	public void setPG(QuatSymmetryResults pG) {
		PG = pG;
	}

	public SymmetryType getType() {
		return type;
	}

	public void setType(SymmetryType type) {
		this.type = type;
	}

	public Atom[] getAtoms() {
		return atoms;
	}

	public void setAtoms(Atom[] atoms) {
		this.atoms = atoms;
	}

}
