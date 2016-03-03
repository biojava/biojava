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
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;
import org.biojava.nbio.structure.symmetry.utils.SymmetryTools;

/**
 * This Class stores all the relevant information of an internal symmetry result
 * obtained with CeSymm. The purpose is to carry all the information packed
 * during the calculations and return a single object.
 * 
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class CeSymmResult {

	private MultipleAlignment multipleAlignment;
	private AFPChain selfAlignment;

	private StructureIdentifier structureId;
	private Atom[] atoms;

	private CESymmParameters params;
	private SymmetryAxes axes;
	private String symmGroup;

	private int symmOrder;
	private int symmLevels;
	private boolean refined;
	private SymmetryType type;

	/**
	 * Conditions checked are: score above the threshold, order higher than 1
	 * and refinement succeeded.
	 * 
	 * @return true if significant, false otherwise
	 */
	public boolean isSignificant() {

		// In any case if the order is 1 it is asymmetric
		if (symmOrder < 2)
			return false;

		// If the TM-Score is below the threshold
		if (selfAlignment.getTMScore() < params.getScoreThreshold())
			return false;

		// If the refinement was attempted
		if (params.getRefineMethod() != RefineMethod.NOT_REFINED) {
			// Check for minimum length
			if (multipleAlignment.getCoreLength() < params.getMinCoreLength())
				return false;
			// Allow 90% of original TM-Score theshold
			if (multipleAlignment.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < params
					.getScoreThreshold() * 0.9)
				return false;
			return true;
		}

		return true;
	}

	/**
	 * Return the symmetric repeats as structure identifiers, if the result is
	 * symmetric and it was refined, return null otherwise.
	 * 
	 * @return List of StructureIdentifiers or null if not defined
	 * @throws StructureException
	 */
	public List<StructureIdentifier> getRepeatsID() throws StructureException {

		if (!isRefined())
			return null;

		List<StructureIdentifier> repeats = new ArrayList<StructureIdentifier>(
				symmOrder);

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

			repeats.add(id);
		}
		return repeats;
	}

	@Override
	public String toString() {
		return structureId + ", symmGroup=" + getSymmGroup() + ", symmOrder="
				+ symmOrder + ", symmLevels=" + symmLevels + ", refined="
				+ refined + ", type=" + type + " | " + params;
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

	public String getSymmGroup() {
		// Lazily calculate the symmetry group if significant
		if (symmGroup == null) {
			if (isSignificant()) {
				if (isRefined()) {
					try {
						symmGroup = SymmetryTools.getQuaternarySymmetry(this)
								.getSymmetry();
					} catch (StructureException e) {
						symmGroup = "C1";
					}
					if (symmGroup.equals("C1"))
						symmGroup = "R"; // could not find group
				} else { 
					// in case significant but not refined
					if (type.equals(SymmetryType.CLOSED))
						symmGroup = "C" + symmOrder;
					else 
						symmGroup = "R";
				}
			} else // case asymmetric
				symmGroup = "C1";
		}
		return symmGroup;
	}

	public void setSymmGroup(String symmGroup) {
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

}
