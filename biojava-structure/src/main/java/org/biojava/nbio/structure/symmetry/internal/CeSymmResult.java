package org.biojava.nbio.structure.symmetry.internal;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.ResidueNumber;
import org.biojava.nbio.structure.ResidueRange;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureIdentifier;
import org.biojava.nbio.structure.SubstructureIdentifier;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetryResults;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;

/**
 * This Class stores all the relevant information of an internal symmetry result
 * obtained with CeSymm. The purpose is to carry all the information packed
 * during the calculations and return a single object.
 * 
 * @author Aleix Lafita
 *
 */
public class CeSymmResult {

	private MultipleAlignment multipleAlignment;
	private AFPChain selfAlignment;

	private StructureIdentifier structureId;
	private Atom[] atoms;

	private CESymmParameters params;
	private SymmetryAxes axes;
	private QuatSymmetryResults symmGroup;

	private int symmOrder;
	private int symmLevels;
	private boolean refined;
	private SymmetryType type;

	/**
	 * Return true if the symmetry result is significant, false otherwise.
	 * 
	 * @return true if significant, false otherwise
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
			// Length is over min length
			if (symmLevels > 1)
				return true; // significance already checked
			else if (multipleAlignment.getCoreLength() < params
					.getMinCoreLength())
				return false;
			else if (multipleAlignment
					.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < params
					.getScoreThreshold())
				return false;
			else
				return true;
		} else
			return false;
	}

	/**
	 * Return the symmetric repeats as structure identifiers.
	 * 
	 * @return List of StructureIdentifiers
	 * @throws StructureException
	 */
	public List<StructureIdentifier> getRepeatsID()
			throws StructureException {

		List<StructureIdentifier> protodomains = new ArrayList<StructureIdentifier>(
				symmOrder);

		if (!isSignificant())
			return protodomains;

		String pdbId = structureId.toCanonical().getPdbId();
		Block align = multipleAlignment.getBlocks().get(0);

		for (int su = 0; su < symmOrder; su++) {
			// Get the start and end residues of the repeat
			ResidueNumber res1 = atoms[align.getStartResidue(su)].getGroup()
					.getResidueNumber();
			ResidueNumber res2 = atoms[align.getFinalResidue(su)].getGroup()
					.getResidueNumber();
			ResidueRange range = new ResidueRange(res1.getChainId(), res1, res2);

			StructureIdentifier id = new SubstructureIdentifier(pdbId,
					Arrays.asList(range));

			protodomains.add(id);
		}
		return protodomains;
	}
	
	@Override
	public String toString() {
		return "CeSymmResult [structureId=" + structureId + ", symmGroup="
				+ symmGroup + ", symmOrder=" + symmOrder + ", symmLevels="
				+ symmLevels + ", refined=" + refined + ", type=" + type + "]";
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

	/**
	 * Return the symmetry order determined by the order detector if the
	 * symmetry is significant. Return 1 otherwise.
	 * 
	 * @return the order of symmetry if the result is significant
	 */
	public int getSymmOrder() {
		if (isSignificant())
			return symmOrder;
		else
			return 1;
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

	public QuatSymmetryResults getSymmGroup() {
		// Lazily calculate the symmetry group
		if (symmGroup == null) {
			try {
				symmGroup = SymmetryTools.getQuaternarySymmetry(this);
			} catch (StructureException e) {
				e.printStackTrace();
			}
		}
		return symmGroup;
	}

	public void setSymmGroup(QuatSymmetryResults symmGroup) {
		this.symmGroup = symmGroup;
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

	public int getSymmLevels() {
		return symmLevels;
	}

	public void setSymmLevels(int symmLevels) {
		this.symmLevels = symmLevels;
	}

	public StructureIdentifier getStructureId() {
		return structureId;
	}

	public void setStructureId(StructureIdentifier structureId) {
		this.structureId = structureId;
	}

	/**
	 * Returns the TM-Score of the symmetry alignment, independently of the
	 * state of the result (MultipleAlignment or AFPChain).
	 * 
	 * @return TM-score
	 */
	public double getTMScore() {
		if (multipleAlignment != null)
			return multipleAlignment
					.getScore(MultipleAlignmentScorer.AVGTM_SCORE);
		else
			return selfAlignment.getTMScore();
	}

	/**
	 * Returns the RMSD of the symmetry alignment, independently of the state of
	 * the result (MultipleAlignment or AFPChain).
	 * 
	 * @return RMSD
	 */
	public double getRMSD() {
		if (multipleAlignment != null)
			return multipleAlignment.getScore(MultipleAlignmentScorer.RMSD);
		else
			return selfAlignment.getTotalRmsdOpt();
	}

}
