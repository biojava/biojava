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

import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * Calls Spencer's method for determining order.
 * This method uses the sequence alignment information
 * to guess the order of symmetry.
 * 
 * @author dmyersturnbull
 * @since 4.2.0
 * 
 */
public class SequenceFunctionOrderDetector implements OrderDetector {

	private int maxSymmetry = 8;
	private float minimumMetricChange = 0.4f;
	
	public SequenceFunctionOrderDetector() {}

	public SequenceFunctionOrderDetector(int maxSymmetry, float minimumMetricChange) {
		this.maxSymmetry = maxSymmetry;
		this.minimumMetricChange = minimumMetricChange;
	}

	@Override
	public int calculateOrder(AFPChain afpChain, Atom[] ca) 
			throws RefinerFailedException {
		try {
			Map<Integer,Integer> alignment = 
					AlignmentTools.alignmentAsMap(afpChain);

			return AlignmentTools.getSymmetryOrder(alignment,
					new AlignmentTools.IdentityMap<Integer>(), 
					maxSymmetry, minimumMetricChange);
			
		} catch (StructureException e) {
			throw new RefinerFailedException(e);
		}
	}

	@Override
	public String toString() {
		return "SequenceFunctionOrderDetector [maxSymmetry=" + maxSymmetry
				+ ", minimumMetricChange=" + minimumMetricChange + "]";
	}
}
