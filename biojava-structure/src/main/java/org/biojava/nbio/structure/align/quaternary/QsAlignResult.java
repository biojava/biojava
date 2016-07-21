package org.biojava.nbio.structure.align.quaternary;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentScorer;
import org.biojava.nbio.structure.cluster.Subunit;

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
			return -1.0;

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

	@Override
	public String toString() {
		return "QsAlignResult [relation=" + relation + ", rmsd=" + getRmsd()
				+ ", length=" + length() + "]";
	}

}
