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

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.StructureAlignmentFactory;
import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;
import org.biojava.nbio.structure.align.model.AFPChain;
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

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
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

	private static final Logger logger = LoggerFactory
			.getLogger(SubunitCluster.class);

	private List<Subunit> subunits = new ArrayList<Subunit>();
	private List<List<Integer>> subunitEQR = new ArrayList<List<Integer>>();
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
	 * A constructor from a single Subunit. To obtain a
	 * SubunitCluster with multiple Subunits, initialize different
	 * SubunitClusters and merge them.
	 * 
	 * @param subunit
	 *            initial Subunit
	 */
	public SubunitCluster(Subunit subunit) {

		subunits.add(subunit);

		List<Integer> identity = new ArrayList<Integer>();
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
		return Collections.unmodifiableList(subunits);
	}

	/**
	 * Tells whether the other SubunitCluster contains exactly the same Subunit.
	 * This is checked by String equality of their residue one-letter sequences.
	 * 
	 * @param other
	 *            SubunitCluster
	 * @return true if the SubunitClusters are identical, false otherwise
	 */
	public boolean isIdenticalTo(SubunitCluster other) {
		String thisSequence = this.subunits.get(this.representative)
				.getProteinSequenceString();
		String otherSequence = other.subunits.get(other.representative)
				.getProteinSequenceString();
		return thisSequence.equals(otherSequence);
	}

	/**
	 * Merges the other SubunitCluster into this one if it contains exactly the
	 * same Subunit. This is checked by {@link #isIdenticalTo(SubunitCluster)}.
	 * 
	 * @param other
	 *            SubunitCluster
	 * @return true if the SubunitClusters were merged, false otherwise
	 */
	public boolean mergeIdentical(SubunitCluster other) {

		if (!isIdenticalTo(other))
			return false;

		logger.info("SubunitClusters are identical");

		this.subunits.addAll(other.subunits);
		this.subunitEQR.addAll(other.subunitEQR);

		return true;
	}

	/**
	 * Merges the other SubunitCluster into this one if their representatives
	 * sequences are similar (according to the criteria in params).
	 * <p>
	 * The sequence alignment is performed using linear {@link SimpleGapPenalty} and
	 * BLOSUM62 as scoring matrix.
	 * 
	 * @param other
	 *            SubunitCluster
	 * @param params
	 *            SubunitClustererParameters, with information whether to use local
	 *            or global alignment, sequence identity and coverage thresholds.
	 *            Threshold values lower than 0.7 are not recommended.
	 *            Use {@link #mergeStructure} for lower values.
	 * @return true if the SubunitClusters were merged, false otherwise
	 * @throws CompoundNotFoundException
	 */

	public boolean mergeSequence(SubunitCluster other, SubunitClustererParameters params) throws CompoundNotFoundException {
		PairwiseSequenceAlignerType alignerType = PairwiseSequenceAlignerType.LOCAL;
		if (params.isUseGlobalMetrics()) {
			alignerType = PairwiseSequenceAlignerType.GLOBAL;
		}
		return mergeSequence(other, params,alignerType
				, new SimpleGapPenalty(),
				SubstitutionMatrixHelper.getBlosum62());
	}

	/**
	 * Merges the other SubunitCluster into this one if their representatives
	 * sequences are similar (according to the criteria in params).
	 * <p>
	 * The sequence alignment is performed using linear {@link SimpleGapPenalty} and
	 * BLOSUM62 as scoring matrix.
	 * 
	 * @param other
	 *            SubunitCluster
	 * @param params
	 *            {@link SubunitClustererParameters}, with information whether to use local
	 *            or global alignment, sequence identity and coverage thresholds.
	 *            Threshold values lower than 0.7 are not recommended.
	 *            Use {@link #mergeStructure} for lower values.
	 * @param alignerType
	 *            parameter for the sequence alignment algorithm
	 * @param gapPenalty
	 *            parameter for the sequence alignment algorithm
	 * @param subsMatrix
	 *            parameter for the sequence alignment algorithm
	 * @return true if the SubunitClusters were merged, false otherwise
	 * @throws CompoundNotFoundException
	 */

	public boolean mergeSequence(SubunitCluster other, SubunitClustererParameters params,
								 PairwiseSequenceAlignerType alignerType,
								 GapPenalty gapPenalty,
								 SubstitutionMatrix<AminoAcidCompound> subsMatrix)
			throws CompoundNotFoundException {

		// Extract the protein sequences as BioJava alignment objects
		ProteinSequence thisSequence = this.subunits.get(this.representative)
				.getProteinSequence();
		ProteinSequence otherSequence = other.subunits
				.get(other.representative).getProteinSequence();

		// Perform the alignment with provided parameters
		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> aligner = Alignments
				.getPairwiseAligner(thisSequence, otherSequence, alignerType,
						gapPenalty, subsMatrix);

		double sequenceIdentity;
		if(params.isUseGlobalMetrics()) {
			sequenceIdentity = aligner.getPair().getPercentageOfIdentity(true);
		} else {
			sequenceIdentity = aligner.getPair().getPercentageOfIdentity(false);
		}

		if (sequenceIdentity < params.getSequenceIdentityThreshold())
			return false;

		double sequenceCoverage = 0;
		if(params.isUseSequenceCoverage()) {
			// Calculate real coverage (subtract gaps in both sequences)
			double gaps1 = aligner.getPair().getAlignedSequence(1)
					.getNumGapPositions();
			double gaps2 = aligner.getPair().getAlignedSequence(2)
					.getNumGapPositions();
			double lengthAlignment = aligner.getPair().getLength();
			double lengthThis = aligner.getQuery().getLength();
			double lengthOther = aligner.getTarget().getLength();
			sequenceCoverage = (lengthAlignment - gaps1 - gaps2)
					/ Math.max(lengthThis, lengthOther);

			if (sequenceCoverage < params.getSequenceCoverageThreshold())
				return false;
		}

		logger.info(String.format("SubunitClusters are similar in sequence "
						+ "with %.2f sequence identity and %.2f coverage", sequenceIdentity,
				sequenceCoverage));

		// If coverage and sequence identity sufficient, merge other and this
		List<Integer> thisAligned = new ArrayList<Integer>();
		List<Integer> otherAligned = new ArrayList<Integer>();

		// Extract the aligned residues of both Subunit
		for (int p = 1; p < aligner.getPair().getLength() + 1; p++) {

			// Skip gaps in any of the two sequences
			if (aligner.getPair().getAlignedSequence(1).isGap(p))
				continue;
			if (aligner.getPair().getAlignedSequence(2).isGap(p))
				continue;

			int thisIndex = aligner.getPair().getIndexInQueryAt(p) - 1;
			int otherIndex = aligner.getPair().getIndexInTargetAt(p) - 1;

			// Only consider residues that are part of the SubunitCluster
			if (this.subunitEQR.get(this.representative).contains(thisIndex)
					&& other.subunitEQR.get(other.representative).contains(
							otherIndex)) {
				thisAligned.add(thisIndex);
				otherAligned.add(otherIndex);
			}
		}

		// Do a List intersection to find out which EQR columns to remove
		List<Integer> thisRemove = new ArrayList<Integer>();
		List<Integer> otherRemove = new ArrayList<Integer>();

		for (int t = 0; t < this.subunitEQR.get(this.representative).size(); t++) {
			// If the index is aligned do nothing, otherwise mark as removing
			if (!thisAligned.contains(this.subunitEQR.get(this.representative)
					.get(t)))
				thisRemove.add(t);
		}

		for (int t = 0; t < other.subunitEQR.get(other.representative).size(); t++) {
			// If the index is aligned do nothing, otherwise mark as removing
			if (!otherAligned.contains(other.subunitEQR.get(
					other.representative).get(t)))
				otherRemove.add(t);
		}
		// Now remove unaligned columns, from end to start
		Collections.sort(thisRemove);
		Collections.reverse(thisRemove);
		Collections.sort(otherRemove);
		Collections.reverse(otherRemove);

		for (int t = 0; t < thisRemove.size(); t++) {
			for (List<Integer> eqr : this.subunitEQR) {
				int column = thisRemove.get(t);
				eqr.remove(column);
			}
		}

		for (int t = 0; t < otherRemove.size(); t++) {
			for (List<Integer> eqr : other.subunitEQR) {
				int column = otherRemove.get(t);
				eqr.remove(column);
			}
		}

		// The representative is the longest sequence
		if (this.subunits.get(this.representative).size() < other.subunits.get(
				other.representative).size())
			this.representative = other.representative + subunits.size();

		this.subunits.addAll(other.subunits);
		this.subunitEQR.addAll(other.subunitEQR);

		this.method = SubunitClustererMethod.SEQUENCE;

		pseudoStoichiometric = !params.isHighConfidenceScores(sequenceIdentity,sequenceCoverage);

		return true;
	}

	/**
	 * Merges the other SubunitCluster into this one if their representative
	 * Atoms are structurally similar (according to the criteria in params).
	 * <p>
	 *
	 * @param other
	 *            SubunitCluster
	 * @param params
	 *            {@link SubunitClustererParameters}, with information on what alignment
	 *            algorithm to use, RMSD/TMScore and structure coverage thresholds.
	 * @return true if the SubunitClusters were merged, false otherwise
	 * @throws StructureException
	 */

	public boolean mergeStructure(SubunitCluster other, SubunitClustererParameters params) throws StructureException {

		StructureAlignment aligner = StructureAlignmentFactory.getAlgorithm(params.getSuperpositionAlgorithm());
		ConfigStrucAligParams aligner_params = aligner.getParameters();

		Method setOptimizeAlignment = null;
		try {
			setOptimizeAlignment = aligner_params.getClass().getMethod("setOptimizeAlignment", boolean.class);
		} catch (NoSuchMethodException e) {
			//alignment algorithm does not have an optimization switch, moving on
		}
		if (setOptimizeAlignment != null) {
			try {
				setOptimizeAlignment.invoke(aligner_params, params.isOptimizeAlignment());
			} catch (IllegalAccessException|InvocationTargetException e) {
				logger.warn("Could not set alignment optimisation switch");
			}
		}

		AFPChain afp = aligner.align(this.subunits.get(this.representative)
				.getRepresentativeAtoms(),
				other.subunits.get(other.representative)
						.getRepresentativeAtoms());

		// Convert AFPChain to MultipleAlignment for convenience
		MultipleAlignment msa = new MultipleAlignmentEnsembleImpl(
				afp,
				this.subunits.get(this.representative).getRepresentativeAtoms(),
				other.subunits.get(other.representative)
						.getRepresentativeAtoms(), false)
				.getMultipleAlignment(0);

		double structureCoverage = Math.min(msa.getCoverages().get(0), msa
				.getCoverages().get(1));

		if(params.isUseStructureCoverage() && structureCoverage < params.getStructureCoverageThreshold()) {
			return false;
		}

		double rmsd = afp.getTotalRmsdOpt();
		if (params.isUseRMSD() && rmsd > params.getRMSDThreshold()) {
			return false;
		}

		double tmScore = afp.getTMScore();
		if (params.isUseTMScore() && tmScore < params.getTMThreshold()) {
			return false;
		}

		logger.info(String.format("SubunitClusters are structurally similar with "
				+ "%.2f RMSD %.2f coverage", rmsd, structureCoverage));

		// Merge clusters
		List<List<Integer>> alignedRes = msa.getBlock(0).getAlignRes();
		List<Integer> thisAligned = new ArrayList<Integer>();
		List<Integer> otherAligned = new ArrayList<Integer>();

		// Extract the aligned residues of both Subunit
		for (int p = 0; p < msa.length(); p++) {

			// Skip gaps in any of the two sequences
			if (alignedRes.get(0).get(p) == null)
				continue;
			if (alignedRes.get(1).get(p) == null)
				continue;

			int thisIndex = alignedRes.get(0).get(p);
			int otherIndex = alignedRes.get(1).get(p);

			// Only consider residues that are part of the SubunitCluster
			if (this.subunitEQR.get(this.representative).contains(thisIndex)
					&& other.subunitEQR.get(other.representative).contains(
							otherIndex)) {
				thisAligned.add(thisIndex);
				otherAligned.add(otherIndex);
			}
		}

		// Do a List intersection to find out which EQR columns to remove
		List<Integer> thisRemove = new ArrayList<Integer>();
		List<Integer> otherRemove = new ArrayList<Integer>();

		for (int t = 0; t < this.subunitEQR.get(this.representative).size(); t++) {
			// If the index is aligned do nothing, otherwise mark as removing
			if (!thisAligned.contains(this.subunitEQR.get(this.representative)
					.get(t)))
				thisRemove.add(t);
		}

		for (int t = 0; t < other.subunitEQR.get(other.representative).size(); t++) {
			// If the index is aligned do nothing, otherwise mark as removing
			if (!otherAligned.contains(other.subunitEQR.get(
					other.representative).get(t)))
				otherRemove.add(t);
		}

		// Now remove unaligned columns, from end to start
		Collections.sort(thisRemove);
		Collections.reverse(thisRemove);
		Collections.sort(otherRemove);
		Collections.reverse(otherRemove);

		for (int t = 0; t < thisRemove.size(); t++) {
			for (List<Integer> eqr : this.subunitEQR) {
				int column = thisRemove.get(t);
				eqr.remove(column);
			}
		}

		for (int t = 0; t < otherRemove.size(); t++) {
			for (List<Integer> eqr : other.subunitEQR) {
				int column = otherRemove.get(t);
				eqr.remove(column);
			}
		}

		// The representative is the longest sequence
		if (this.subunits.get(this.representative).size() < other.subunits.get(
				other.representative).size())
			this.representative = other.representative + subunits.size();

		this.subunits.addAll(other.subunits);
		this.subunitEQR.addAll(other.subunitEQR);

		this.method = SubunitClustererMethod.STRUCTURE;
		pseudoStoichiometric = true;

		return true;
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

		List<List<Integer>> columns = new ArrayList<List<Integer>>();
		for (int s = 0; s < alignedRes.size(); s++)
			columns.add(new ArrayList<Integer>(alignedRes.get(s).size()));

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

		List<Atom[]> alignedAtoms = Collections.emptyList();

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
