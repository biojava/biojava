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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.cluster.Subunit;
import org.biojava.nbio.structure.cluster.SubunitCluster;

/**
 * Result of a Quaternary Structure Alignment {@link QsAlign}. The QsAlignResult
 * holds the original inputs of the algorithm and the results and scores of the
 * alignment.
 * 
 * @author Aleix Lafita
 * @since 5.0.0
 *
 */
public class QsAlignResult {

	private List<SubunitCluster> clusters;

	private List<Subunit> subunits1;
	private List<Subunit> subunits2;

	private Map<Integer, Integer> subunitMap;
	private MultipleAlignment alignment;

	private QsRelation relation;

	/**
	 * The Constructor of the result takes the same inputs as the
	 * {@link QsAlign} algorithm.
	 * 
	 * @param subunits1
	 * @param subunits2
	 */
	public QsAlignResult(List<Subunit> subunits1, List<Subunit> subunits2) {

		this.subunits1 = subunits1;
		this.subunits2 = subunits2;

		subunitMap = Collections.emptyMap();
		relation = QsRelation.DIFFERENT;

	}

	/**
	 * Original Subunits of the first group.
	 * 
	 * @return an unmodifiable view of the original List
	 */
	public List<Subunit> getSubunits1() {
		return Collections.unmodifiableList(subunits1);
	}

	/**
	 * Original Subunits of the second group.
	 * 
	 * @return an unmodifiable view of the original List
	 */
	public List<Subunit> getSubunits2() {
		return Collections.unmodifiableList(subunits2);
	}

	/**
	 * Map of Subunit equivalencies from the first to the second group.
	 * 
	 * @return an unmodifiable view of the original Map
	 */
	public Map<Integer, Integer> getSubunitMap() {

		if (subunitMap == null)
			return Collections.emptyMap();

		return Collections.unmodifiableMap(subunitMap);
	}

	/**
	 * Map of Subunit equivalencies from the first to the second group.
	 * 
	 * @param subunitMap
	 */
	public void setSubunitMap(Map<Integer, Integer> subunitMap) {

		// Check consistency of the map
		if (Collections.max(subunitMap.keySet()) > subunits1.size()
				| Collections.max(subunitMap.values()) > subunits2.size())
			throw new IndexOutOfBoundsException(
					"Subunit Map index higher than Subunit List size.");

		// Update the relation enum
		if (subunitMap.size() == 0) {
			relation = QsRelation.DIFFERENT;
		} else if (subunitMap.keySet().size() == subunits1.size()) {
			if (subunitMap.values().size() == subunits2.size()) {
				relation = QsRelation.EQUIVALENT;
			} else {
				relation = QsRelation.PARTIAL_COMPLETE;
			}
		} else {
			if (subunitMap.values().size() == subunits2.size()) {
				relation = QsRelation.PARTIAL_COMPLETE;
			} else {
				relation = QsRelation.PARTIAL_INCOMPLETE;
			}
		}

		this.subunitMap = subunitMap;
	}

	/**
	 * The length of the alignment is the number of Subunit equivalencies it
	 * contains. This is equivalent to the size of the Subunit Map.
	 * 
	 * @return length of the alignment
	 */
	public int length() {
		if (subunitMap == null)
			return 0;

		return subunitMap.size();
	}

	/**
	 * The transformation 4D matrix that needs to be applied to the second group
	 * of Subunits to superimpose them onto the first group of Subunits, given
	 * the equivalent residues in the SubunitCluster and the Subunit
	 * equivalencies.
	 * <p>
	 * This is equivalent to
	 * multipleAlignment.getBlockSet(0).getTransformations().get(1).
	 * 
	 * @return Matrix4d
	 */
	public Matrix4d getTransform() {

		if (alignment == null)
			return null;

		return alignment.getBlockSet(0).getTransformations().get(1);
	}

	/**
	 * The RMSD between the equivalent residues of the equivalent Subunits after
	 * superposition of the Subunit groups. This is equivalent to
	 * multipleAlignment.getScore(MultipleAlignmentScorer.RMSD).
	 * 
	 * @return rmsd
	 */
	public double getRmsd() {

		if (alignment == null)
			return -1.0;
		if (alignment.getScore(MultipleAlignmentScorer.RMSD) == null)
			return MultipleAlignmentScorer.getRMSD(alignment);

		return alignment.getScore(MultipleAlignmentScorer.RMSD);
	}

	/**
	 * The quaternary structure relation {@link QsRelation} between the two
	 * groups of Subunits.
	 * 
	 * @return relation
	 */
	public QsRelation getRelation() {
		return relation;
	}

	/**
	 * The quaternary structure relation {@link QsRelation} between the two
	 * groups of Subunits.
	 * 
	 * @param relation
	 */
	public void setRelation(QsRelation relation) {
		this.relation = relation;
	}

	/**
	 * The alignment that specifies the residue equivalencies of the equivalent
	 * Subunits.
	 * 
	 * @return alignment as a MultipleAlignment object
	 */
	public MultipleAlignment getAlignment() {
		return alignment;
	}

	/**
	 * The alignment that specifies the residue equivalencies of the equivalent
	 * Subunits.
	 * 
	 * @param alignment
	 *            a MultipleAlignment object
	 */
	public void setAlignment(MultipleAlignment alignment) {
		this.alignment = alignment;
	}

	/**
	 * Return the aligned subunits of the first Subunit group, in the alignment
	 * order.
	 * 
	 * @return a List of Subunits in the alignment order
	 */
	public List<Subunit> getAlignedSubunits1() {

		List<Subunit> aligned = new ArrayList<Subunit>(subunitMap.size());

		for (Integer key : subunitMap.keySet())
			aligned.add(subunits1.get(key));

		return aligned;
	}

	/**
	 * Return the aligned subunits of the second Subunit group, in the alignment
	 * order.
	 * 
	 * @return a List of Subunits in the alignment order
	 */
	public List<Subunit> getAlignedSubunits2() {

		List<Subunit> aligned = new ArrayList<Subunit>(subunitMap.size());

		for (Integer key : subunitMap.keySet())
			aligned.add(subunits2.get(subunitMap.get(key)));

		return aligned;
	}

	public void setClusters(List<SubunitCluster> clusters) {
		this.clusters = clusters;
	}

	public Atom[] getAlignedAtomsForSubunits1(int index) {

		// Obtain the indices of the clustered subunits
		for (SubunitCluster cluster : clusters) {
			if (cluster.getSubunits().contains(subunits1.get(index))) {
				return cluster.getAlignedAtomsSubunit(cluster.getSubunits()
						.indexOf(subunits1.get(index)));
			}
		}
		return null;
	}

	public Atom[] getAlignedAtomsForSubunits2(int index) {

		// Obtain the indices of the clustered subunits
		for (SubunitCluster cluster : clusters) {
			if (cluster.getSubunits().contains(subunits2.get(index))) {
				return cluster.getAlignedAtomsSubunit(cluster.getSubunits()
						.indexOf(subunits2.get(index)));
			}
		}
		return null;
	}

	@Override
	public String toString() {
		return "QsAlignResult [relation="
				+ relation
				+ ", rmsd="
				+ getRmsd()
				+ ", length="
				+ length()
				+ ", Aligned 1: "
				+ getAlignedSubunits1().stream().map(s -> s.getName())
						.collect(Collectors.toList())
				+ ", Aligned 2: "
				+ getAlignedSubunits2().stream().map(s -> s.getName())
						.collect(Collectors.toList()) + "]";
	}

}
