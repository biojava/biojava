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
	private SymmetryAxes axes = new SymmetryAxes();
	private String symmGroup;

	private int numRepeats;
	private boolean refined;

	/**
	 * Conditions checked are: score above the threshold, number of repeats
	 * higher than 1 and refinement succeeded.
	 *
	 * @return true if significant, false otherwise
	 */
	public boolean isSignificant() {

		// In any case if the order is 1 it is asymmetric
		if (numRepeats < 2)
			return false;

		// If the TM-Score before refinement is below the threshold
		if (selfAlignment.getTMScore() < params.getUnrefinedScoreThreshold())
			return false;

		// If the refinement was attempted
		if (params.getRefineMethod() != RefineMethod.NOT_REFINED) {
			// Check for minimum length
			if (multipleAlignment == null || multipleAlignment.getCoreLength() < params.getMinCoreLength())
				return false;
			// Allow 90% of original TM-Score theshold
			if (multipleAlignment.getScore(MultipleAlignmentScorer.AVGTM_SCORE) < params
					.getRefinedScoreThreshold())
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
				numRepeats);

		String pdbId = structureId.toCanonical().getPdbId();
		Block align = multipleAlignment.getBlocks().get(0);

		for (int su = 0; su < numRepeats; su++) {
			// Get the start and end residues of the repeat
			ResidueNumber res1 = atoms[align.getStartResidue(su)].getGroup()
					.getResidueNumber();
			ResidueNumber res2 = atoms[align.getFinalResidue(su)].getGroup()
					.getResidueNumber();
			ResidueRange range = new ResidueRange(res1.getChainName(), res1, res2);

			StructureIdentifier id = new SubstructureIdentifier(pdbId,
					Arrays.asList(range));

			repeats.add(id);
		}
		return repeats;
	}

	@Override
	public String toString() {
		return structureId + ", symmGroup=" + getSymmGroup() + ", numRepeats="
				+ numRepeats + ", symmLevels=" + axes.getNumLevels()
				+ ", refined=" + refined + " | " + params;
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
	public int getNumRepeats() {
		if (isSignificant())
			return numRepeats;
		else
			return 1;
	}

	public void setNumRepeats(int symmOrder) {
		this.numRepeats = symmOrder;
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
					if (axes.getElementaryAxis(0).getSymmType()
							.equals(SymmetryType.CLOSED))
						symmGroup = "C" + numRepeats;
					else
						symmGroup = "R";
				}
			} else
				// case asymmetric
				symmGroup = "C1";
		}
		return symmGroup;
	}

	public void setSymmGroup(String symmGroup) {
		this.symmGroup = symmGroup;
	}

	public Atom[] getAtoms() {
		return atoms;
	}

	public void setAtoms(Atom[] atoms) {
		this.atoms = atoms;
	}

	public int getSymmLevels() {
		return axes.getNumLevels();
	}

	public StructureIdentifier getStructureId() {
		return structureId;
	}

	public void setStructureId(StructureIdentifier structureId) {
		this.structureId = structureId;
	}

	/**
	 * Return a String describing the reasons for the CE-Symm final decision in
	 * this particular result.
	 * 
	 * @return String decision reason
	 */
	public String getReason() {
		// Cases:
		// 1. Asymmetric because insignificant self-alignment (1itb.A_1-100)
		double tm = selfAlignment.getTMScore();
		if (tm < params.getUnrefinedScoreThreshold()) {
			return String.format("Insignificant self-alignment (TM=%.2f)", tm);
		}
		// 2. Asymmetric because order detector returned 1
		if (numRepeats == 1) {
			return String.format(
					"Order detector found asymmetric alignment (TM=%.2f)", tm);
		}

		// Check that the user requested refinement
		if (params.getRefineMethod() != RefineMethod.NOT_REFINED) {
			// 3. Asymmetric because refinement failed
			if (!refined) {
				return "Refinement failed";
			}
			tm = multipleAlignment
					.getScore(MultipleAlignmentScorer.AVGTM_SCORE);
			// 4. Asymmetric because refinement & optimization were not
			// significant
			if (!isSignificant()) {
				return String.format(
						"Refinement was not significant (TM=%.2f)", tm);
			}
		} else {
			// 4. Not refined, but result was not significant
			if (!isSignificant()) {
				return String
						.format("Result was not significant (TM=%.2f)", tm);
			}
		}

		String hierarchical = "";
		if (axes.getNumLevels() > 1) {
			hierarchical = String.format("; Contains %d levels of symmetry",
					axes.getNumLevels());
		}
		// 5. Symmetric.
		// a. Open. Give # repeats (1n0r.A)
		if (axes.getElementaryAxis(0).getSymmType() == SymmetryType.OPEN) {
			return String.format("Contains %d open repeats (TM=%.2f)%s",
					getNumRepeats(), tm, hierarchical);
		}
		// b. Closed, non-hierarchical (1itb.A)
		// c. Closed, heirarchical (4gcr)
		return String.format("Significant (TM=%.2f)%s", tm, hierarchical);
	}

}
