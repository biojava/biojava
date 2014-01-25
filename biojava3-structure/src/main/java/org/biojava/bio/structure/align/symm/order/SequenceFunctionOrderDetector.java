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
 * @since 3.0.8
 */
package org.biojava.bio.structure.align.symm.order;

import java.util.Map;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.align.model.AFPChain;
import org.biojava.bio.structure.align.util.AlignmentTools;

/**
 * Tries to detect symmetry in an alignment.
 *
 * <p>Conceptually, an alignment is a function f:A->B between two sets of
 * integers. The function may have simple topology (meaning that if two
 * elements of A are close, then their images in B will also be close), or
 * may have more complex topology (such as a circular permutation). This
 * function checks <i>alignment</i> against a reference function
 * <i>identity</i>, which should have simple topology. It then tries to
 * determine the symmetry order of <i>alignment</i> relative to
 * <i>identity</i>, up to a maximum order of <i>maxSymmetry</i>.
 *
 *
 * <p><strong>Details</strong><br/>
 * Considers the offset (in number of residues) which a residue moves
 * after undergoing <i>n</i> alternating transforms by alignment and
 * identity. If <i>n</i> corresponds to the intrinsic order of the alignment,
 * this will be small. This algorithm tries increasing values of <i>n</i>
 * and looks for abrupt decreases in the root mean squared offset.
 * If none are found at <i>n</i><=maxSymmetry, the alignment is reported as
 * non-symmetric.
 *
 * @param alignment The alignment to test for symmetry
 * @param identity An alignment with simple topology which approximates
 *  the sequential relationship between the two proteins. Should map in the
 *  reverse direction from alignment.
 * @param maxSymmetry Maximum symmetry to consider. High values increase
 *  the calculation time and can lead to overfitting.
 * @param minimumMetricChange Percent decrease in root mean squared offsets
 *  in order to declare symmetry. 0.4f seems to work well for CeSymm.
 * @return The order of symmetry of alignment, or 1 if no order <=
 *  maxSymmetry is found.
 *
 * @author sbliven
 */
public class SequenceFunctionOrderDetector implements OrderDetector {

	private int maxSymmetry = 8;
	private float minimumMetricChange = 0.4f;
	
	public SequenceFunctionOrderDetector() {
		super();
	}

	public SequenceFunctionOrderDetector(int maxSymmetry, float minimumMetricChange) {
		super();
		this.maxSymmetry = maxSymmetry;
		this.minimumMetricChange = minimumMetricChange;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) throws OrderDetectionFailedException {
		try {
			Map<Integer,Integer> alignment = AlignmentTools.alignmentAsMap(afpChain);

			return AlignmentTools.getSymmetryOrder(alignment,
					new AlignmentTools.IdentityMap<Integer>(), maxSymmetry, minimumMetricChange);
//			return AlignmentTools.getSymmetryOrder(afpChain, maxSymmetry, minimumMetricChange);
		} catch (StructureException e) {
			throw new OrderDetectionFailedException(e);
		}
	}

}
