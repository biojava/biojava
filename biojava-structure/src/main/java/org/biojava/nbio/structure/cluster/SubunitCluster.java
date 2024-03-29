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
package org.biojava.nbio.structure.cluster;

import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.multiple.Block;
import org.biojava.nbio.structure.align.multiple.BlockImpl;
import org.biojava.nbio.structure.align.multiple.BlockSet;
import org.biojava.nbio.structure.align.multiple.BlockSetImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.align.multiple.util.ReferenceSuperimposer;
import org.biojava.nbio.structure.symmetry.core.QuatSymmetrySubunits;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters;
import org.biojava.nbio.structure.symmetry.internal.CeSymm;
import org.biojava.nbio.structure.symmetry.internal.CeSymmResult;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;

/**
 * A SubunitCluster contains a set of equivalent {@link QuatSymmetrySubunits},
 * the set of equivalent residues (EQR) between {@link Subunit} and a
 * {@link Subunit} representative. It also stores the method used for
 * clustering.
 * <p>
 * This class allows the comparison and merging of SubunitClusters.
 *
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class SubunitCluster {

	private static final Logger logger = LoggerFactory.getLogger(SubunitCluster.class);

	private List<Subunit> subunits = new ArrayList<>();
	private List<List<Integer>> subunitEQR = new ArrayList<>();
	private int representative = -1;

	private SubunitClustererMethod method = SubunitClustererMethod.SEQUENCE;
	private boolean pseudoStoichiometric = false;

	/**
	 * A letter that is assigned to this cluster in stoichiometry.
	*/
	private String alpha = "";

	/**
	 * A letter that is assigned to this cluster in stoichiometry.
	 *
	 * @return alpha
	 *          String
	 */

	public String getAlpha() {
		return alpha;
	}

	/**
	 * A letter that is assigned to this cluster in stoichiometry.
	 *
	 * @param  alpha
	 *          String
	 */
	public void setAlpha(String alpha) {
		this.alpha = alpha;
	}

	/**
	 * Setter for method, to be used by SubunitClusterMerge.
	 *
	 * @param  method
	 *          SubunitClustererMethod
	 */
	public void setMethod(SubunitClustererMethod method) {
		this.method = method;
	}

	/**
	 * Setter for pseudoStoichiometric, to be used by SubunitClusterMerge.
	 *
	 * @param  pseudoStoichiometric
	 *          boolean
	 */
	public void setPseudoStoichiometric(boolean pseudoStoichiometric) {
		this.pseudoStoichiometric = pseudoStoichiometric;
	}

	/**
	 * Setter for representative, to be used by SubunitClusterMerge.
	 *
	 * @param  representative
	 *          boolean
	 */
	public void setRepresentative(int representative) {
		this.representative = representative;
	}

	/**
	 * A constructor from a single Subunit. To obtain a
	 * SubunitCluster with multiple Subunits, initialize different
	 * SubunitClusters and merge them.
	 *
	 * @param subunit
	 *            initial Subunit
	 */
	public SubunitCluster(Subunit subunit) {

		subunits.add(subunit);

		List<Integer> identity = new ArrayList<>();
		for (int i = 0; i < subunit.size(); i++)
			identity.add(i);
		subunitEQR.add(identity);

		representative = 0;
	}

	/**
	 * A copy constructor with the possibility of removing subunits.
	 * No re-clustering is done.
	 *
	 * @param other
	 *            reference SubunitCluster
	 * @param subunitsToRetain
	 *            which subunits to copy to this cluster
	 */
	public SubunitCluster(SubunitCluster other, List<Integer> subunitsToRetain) {
		method = other.method;
		pseudoStoichiometric = other.pseudoStoichiometric;
		for (int i = 0; i < other.subunits.size(); i++) {
			if(subunitsToRetain.contains(i)) {
				subunits.add(other.subunits.get(i));
				subunitEQR.add(other.subunitEQR.get(i));
			}
		}
		representative = 0;
		for (int i=1; i<subunits.size(); i++) {
			if (subunits.get(i).size() > subunits.get(representative).size()) {
				representative = i;
			}
		}
		setAlpha(other.getAlpha());
	}

	/**
	 * Subunits contained in the SubunitCluster.
	 *
	 * @return an unmodifiable view of the original List
	 */
	public List<Subunit> getSubunits() {
		return subunits;
	}

	/**
	 * Analyze the internal symmetry of the SubunitCluster and divide its
	 * {@link Subunit} into the internal repeats (domains) if they are
	 * internally symmetric.
	 *
	 * @param clusterParams {@link SubunitClustererParameters} with fields used as follows:
	 * structureCoverageThreshold
	 *            the minimum coverage of all repeats in the Subunit
	 * rmsdThreshold
	 *            the maximum allowed RMSD between the repeats
	 * minimumSequenceLength
	 *            the minimum length of the repeating units
	 * @return true if the cluster was internally symmetric, false otherwise
	 * @throws StructureException
	 */
	public boolean divideInternally(SubunitClustererParameters clusterParams)
			throws StructureException {

		CESymmParameters cesym_params = new CESymmParameters();
		cesym_params.setMinCoreLength(clusterParams.getMinimumSequenceLength());
		cesym_params.setGaps(false); // We want no gaps between the repeats

		// Analyze the internal symmetry of the representative subunit
		CeSymmResult result = CeSymm.analyze(subunits.get(representative)
				.getRepresentativeAtoms(), cesym_params);

		if (!result.isSignificant())
			return false;

		double rmsd = result.getMultipleAlignment().getScore(
				MultipleAlignmentScorer.RMSD);
		if (rmsd > clusterParams.getRMSDThreshold())
			return false;

		double coverage = result.getMultipleAlignment().getCoverages().get(0)
				* result.getNumRepeats();
		if (coverage < clusterParams.getStructureCoverageThreshold())
			return false;

		logger.info("SubunitCluster is internally symmetric with {} repeats, "
				+ "{} RMSD and {} coverage", result.getNumRepeats(), rmsd,
				coverage);

		// Divide if symmety was significant with RMSD and coverage sufficient
		List<List<Integer>> alignedRes = result.getMultipleAlignment()
				.getBlock(0).getAlignRes();

		List<List<Integer>> columns = new ArrayList<>();
		for (int s = 0; s < alignedRes.size(); s++)
			columns.add(new ArrayList<>(alignedRes.get(s).size()));

		// Extract the aligned columns of each repeat in the Subunit
		for (int col = 0; col < alignedRes.get(0).size(); col++) {

			// Check that all aligned residues are part of the Cluster
			boolean missing = false;
			for (int s = 0; s < alignedRes.size(); s++) {
				if (!subunitEQR.get(representative).contains(
						alignedRes.get(s).get(col))) {
					missing = true;
					break;
				}
			}

			// Skip the column if any residue was not part of the cluster
			if (missing)
				continue;

			for (int s = 0; s < alignedRes.size(); s++) {
				columns.get(s).add(
						subunitEQR.get(representative).indexOf(
								alignedRes.get(s).get(col)));
			}
		}

		// Divide the Subunits in their repeats
		List<Subunit> newSubunits = new ArrayList<Subunit>(subunits.size()
				* columns.size());
		List<List<Integer>> newSubunitEQR = new ArrayList<List<Integer>>(
				subunits.size() * columns.size());

		for (int s = 0; s < subunits.size(); s++) {
			for (int r = 0; r < columns.size(); r++) {

				// Calculate start and end residues of the new Subunit
				int start = subunitEQR.get(s).get(columns.get(r).get(0));
				int end = subunitEQR.get(s).get(
						columns.get(r).get(columns.get(r).size() - 1));

				Atom[] reprAtoms = Arrays.copyOfRange(subunits.get(s)
						.getRepresentativeAtoms(), start, end + 1);

				newSubunits.add(new Subunit(reprAtoms, subunits.get(s)
						.getName(), subunits.get(s).getIdentifier(), subunits
						.get(s).getStructure()));

				// Recalculate equivalent residues
				List<Integer> eqr = new ArrayList<Integer>();
				for (int p = 0; p < columns.get(r).size(); p++) {
					eqr.add(subunitEQR.get(s).get(columns.get(r).get(p))
							- start);
				}
				newSubunitEQR.add(eqr);
			}
		}

		subunits = newSubunits;
		subunitEQR = newSubunitEQR;

		// Update representative
		for (int s = 0; s < subunits.size(); s++) {
			if (subunits.get(s).size() > subunits.get(representative).size())
				representative = s;
		}

		method = SubunitClustererMethod.STRUCTURE;
		pseudoStoichiometric = true;
		return true;
	}

	/**
	 * @return the number of Subunits in the cluster
	 */
	public int size() {
		return subunits.size();
	}

	/**
	 * @return the number of aligned residues between Subunits of the cluster
	 */
	public int length() {
		return subunitEQR.get(representative).size();
	}

	/**
	 * @return the representative number
	 */
	public int getRepresentative() {
		return representative;
	}

	/**
	 * @return the subunitEQR list
	 */
	public List<List<Integer>> getSubunitEQR() {
		return subunitEQR;
	}

	/**
	 * @return the {@link SubunitClustererMethod} used for clustering the
	 *         Subunits
	 */
	public SubunitClustererMethod getClustererMethod() {
		return method;
	}

	/**
	 * @return A List of size {@link #size()} of Atom arrays of length
	 *         {@link #length()} with the aligned Atoms for each Subunit in the
	 *         cluster
	 */
	public List<Atom[]> getAlignedAtomsSubunits() {

		List<Atom[]> alignedAtoms = new ArrayList<>();

		// Loop through all subunits and add the aligned positions
		for (int s = 0; s < subunits.size(); s++)
			alignedAtoms.add(getAlignedAtomsSubunit(s));

		return alignedAtoms;
	}

	/**
	 * @param index
	 *            Subunit index in the Cluster
	 * @return An Atom array of length {@link #length()} with the aligned Atoms
	 *         from the selected Subunit in the Cluster
	 */
	public Atom[] getAlignedAtomsSubunit(int index) {

		Atom[] aligned = new Atom[subunitEQR.get(index).size()];

		// Add only the aligned positions of the Subunit in the Cluster
		for (int p = 0; p < subunitEQR.get(index).size(); p++) {
			aligned[p] = subunits.get(index).getRepresentativeAtoms()[subunitEQR
					.get(index).get(p)];
		}

		return aligned;
	}

	/**
	 * The multiple alignment is calculated from the equivalent residues in the
	 * SubunitCluster. The alignment is recalculated every time the method is
	 * called (no caching).
	 *
	 * @return MultipleAlignment representation of the aligned residues in this
	 *         Subunit Cluster
	 * @throws StructureException
	 */
	public MultipleAlignment getMultipleAlignment() throws StructureException {

		// Create a multiple alignment with the atom arrays of the Subunits
		MultipleAlignment msa = new MultipleAlignmentImpl();
		msa.setEnsemble(new MultipleAlignmentEnsembleImpl());
		msa.getEnsemble().setAtomArrays(
				subunits.stream().map(s -> s.getRepresentativeAtoms())
						.collect(Collectors.toList()));

		// Fill in the alignment information
		BlockSet bs = new BlockSetImpl(msa);
		Block b = new BlockImpl(bs);
		b.setAlignRes(subunitEQR);

		// Fill in the transformation matrices
		new ReferenceSuperimposer(representative).superimpose(msa);

		// Calculate some scores
		MultipleAlignmentScorer.calculateScores(msa);

		return msa;

	}

	@Override
	public String toString() {
		return "SubunitCluster [Size=" + size() + ", Length=" + length()
				+ ", Representative=" + representative + ", Method=" + method
				+ "]";
	}

	/**
	 * @return true if this cluster is considered pseudo-stoichiometric (i.e.,
	 * 		   was either clustered by structure, or by sequence with low scores),
	 *         false otherwise.
	 */
	public boolean isPseudoStoichiometric() {
		return pseudoStoichiometric;
	}

}
