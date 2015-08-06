package org.biojava.nbio.structure.symmetry.internal;

import java.util.Map;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AlignmentTools;

/**
 * Calls Spencer's method for determining order.
 * 
 * @author dmyersturnbull
 * @since 4.2.0
 * 
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
	
	//Order Detector Fails: "d1osya_", "d1xq4a_", "d1yioa2", "d1dcea2"

	@Override
	public String toString() {
		return "SequenceFunctionOrderDetector [maxSymmetry=" + maxSymmetry
				+ ", minimumMetricChange=" + minimumMetricChange + "]";
	}
}
