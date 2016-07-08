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

import java.util.List;
import java.util.Random;

import org.biojava.nbio.structure.align.ce.CeParameters;

/**
 * Provides parameters to {@link CeSymm}.
 *
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.1.1
 *
 */
public class CESymmParameters extends CeParameters {

	private int maxSymmOrder;
	private int userOrder;
	private SymmetryType symmType;
	private OrderDetectorMethod orderDetectorMethod;
	private RefineMethod refineMethod;
	private boolean optimization;
	private int rndSeed;
	private int symmLevels;
	private double unrefinedScoreThreshold;
	private double refinedScoreThreshold;
	private int sseThreshold;
	private int minCoreLength;
	private double distanceCutoff;
	private boolean gaps;
	private int optimizationSteps;

	public static enum OrderDetectorMethod {
		SEQUENCE_FUNCTION, GRAPH_COMPONENT, ANGLE, USER_INPUT;
		public static final OrderDetectorMethod DEFAULT = SEQUENCE_FUNCTION;
	}

	public static enum RefineMethod {
		NOT_REFINED, SEQUENCE_FUNCTION, GRAPH_COMPONENT;
		public static final RefineMethod DEFAULT = SEQUENCE_FUNCTION;
	}

	public static final double DEFAULT_SYMMETRY_THRESHOLD = 0.4;

	/**
	 * The internal symmetry detection can be divided into two types: CLOSE:
	 * includes the circular and dihedral symmetries, and OPEN: includes the
	 * helical and protein repeats symmetries.
	 * <p>
	 * All internal symmetry cases share one property: all the repeats have the
	 * same 3D transformation.
	 * <p>
	 * AUTO option automatically identifies the type. The criterion for
	 * classification is that the CLOSE symmetry generates CeSymm alignments
	 * with circular permutations (2 blocks in AFPChain), whereas the OPEN
	 * symmetry generates alignments without a CP (only one block in AFPChain).
	 */
	public enum SymmetryType {
		CLOSED, OPEN, AUTO;
		public static final SymmetryType DEFAULT = AUTO;
	}

	public CESymmParameters() {
		reset();
	}

	@Override
	public CESymmParameters clone() {
		return new CESymmParameters(this);
	}
	
	public CESymmParameters(CESymmParameters o) {
		this.maxSymmOrder = o.maxSymmOrder;
		this.symmType = o.symmType;
		this.orderDetectorMethod = o.orderDetectorMethod;
		this.userOrder = o.userOrder;
		this.refineMethod = o.refineMethod;
		this.optimization = o.optimization;
		this.rndSeed = o.rndSeed;
		this.symmLevels = o.symmLevels;
		this.unrefinedScoreThreshold = o.unrefinedScoreThreshold;
		this.refinedScoreThreshold = o.refinedScoreThreshold;
		this.sseThreshold = o.sseThreshold;
		this.minCoreLength = o.minCoreLength;
		this.distanceCutoff = o.distanceCutoff;
		this.gaps = o.gaps;
		this.optimizationSteps = o.optimizationSteps;

		this.winSize = o.winSize;
		this.rmsdThr = o.rmsdThr;
		this.rmsdThrJoin = o.rmsdThrJoin;
		this.scoringStrategy = o.scoringStrategy;
		this.maxGapSize = o.maxGapSize;
		this.showAFPRanges = o.showAFPRanges;
		this.maxOptRMSD = o.maxOptRMSD;
		this.gapOpen = o.gapOpen;
		this.gapExtension = o.gapExtension;
		this.distanceIncrement = o.distanceIncrement;
		this.oRmsdThr = o.oRmsdThr;
		this.maxNrIterationsForOptimization = o.maxNrIterationsForOptimization;
		this.seqWeight = o.seqWeight;
	}

	@Override
	public void reset() {
		super.reset();
		maxSymmOrder = 8;
		symmType = SymmetryType.DEFAULT;
		orderDetectorMethod = OrderDetectorMethod.DEFAULT;
		userOrder = 0;
		refineMethod = RefineMethod.DEFAULT;
		optimization = true;
		rndSeed = new Random().nextInt(10000);
		symmLevels = 0;
		unrefinedScoreThreshold = DEFAULT_SYMMETRY_THRESHOLD;
		refinedScoreThreshold = DEFAULT_SYMMETRY_THRESHOLD * 0.9;
		sseThreshold = 0;
		minCoreLength = 15;
		distanceCutoff = 7.0;
		gaps = true;
		optimizationSteps = 0;
	}

	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = super.getUserConfigHelp();

		// maxSymmOrder help explanation
		params.add("Sets the maximum order of symmetry of the protein.");

		// userOrder help explanation
		params.add("Order of symmetry determined by the user. "
				+ "Use it with the USER_INPUT order option. Imposes an order"
				+ " of symmetry to the alignment. If 0 the order is set "
				+ "automatically.");

		StringBuilder symmTypes = new StringBuilder("Type of Symmetry: ");
		SymmetryType[] vals = SymmetryType.values();
		if (vals.length == 1) {
			symmTypes.append(vals[0].name());
		} else if (vals.length > 1) {
			for (int i = 0; i < vals.length - 1; i++) {
				symmTypes.append(vals[i].name());
				symmTypes.append(", ");
			}
			symmTypes.append("or ");
			symmTypes.append(vals[vals.length - 1].name());
		}
		params.add(symmTypes.toString());

		StringBuilder orderTypes = new StringBuilder("Order Detection Method: ");
		OrderDetectorMethod[] vals2 = OrderDetectorMethod.values();
		if (vals2.length == 1) {
			orderTypes.append(vals2[0].name());
		} else if (vals2.length > 1) {
			for (int i = 0; i < vals2.length - 1; i++) {
				orderTypes.append(vals2[i].name());
				orderTypes.append(", ");
			}
			orderTypes.append("or ");
			orderTypes.append(vals[vals.length - 1].name());
		}
		params.add(orderTypes.toString());

		StringBuilder refineTypes = new StringBuilder("Refinement Method: ");
		RefineMethod[] values = RefineMethod.values();
		if (values.length == 1) {
			refineTypes.append(values[0].name());
		} else if (values.length > 1) {
			for (int i = 0; i < values.length - 1; i++) {
				refineTypes.append(values[i].name());
				refineTypes.append(", ");
			}
			refineTypes.append("or ");
			refineTypes.append(values[values.length - 1].name());
		}
		params.add(refineTypes.toString());

		// optimization help explanation
		params.add("Optimize the refined alignment if true.");
		// seed help explanation
		params.add("Random seed for the Monte Carlo optimization, "
				+ "for reproducibility of results.");
		// symmetry levels
		params.add("Specify the maximum number of symmetry levels to explore "
				+ "recursively. If equal to 1, only C and H symmetries can be "
				+ "found. If equal to 2, D and two-level hierarchical C and H "
				+ "can be found, etc. If equal to 0, the number of recursive "
				+ "iterations is unbounded (until thresholds reached).");
		// unrefined score threshold
		params.add("Unrefined score threshold: TM-score values for the optimal"
				+ " self-alignment, before refinement, below the "
				+ "threshold will be considered asymmetric.");
		// refined score threshold
		params.add("Refined score threshold: TM-score values for the refined "
				+ "multiple alignment of repeats below the "
				+ "threshold will be considered asymmetric.");
		// SSE threshold
		params.add("SSE threshold: The minimum number of secondary structure "
				+ "elements (strands or helices) in each symmetrical repeat. "
				+ "If the repeats do not have enough SSE, the structure will "
				+ "be considered asymmetric. 0 means no restriction.");
		// min core repeat length
		params.add("Minimum core length: the minimum number of non-gapped "
				+ "residues in every symmetric repeat.");
		// distance cutoff
		params.add("Distance Cutoff: the maximum allowed distance (in A) "
				+ "between two aligned residues.");

		// gaps
		params.add("Internal Gaps: allow up to 50% of repeats to have gaps in "
				+ "the multiple alignment if true, "
				+ "otherwise all repeats must be aligned at each position.");

		// optimization steps
		params.add("Optimization Steps: maximum number of optimization steps:"
				+ " 0 means calculated automatically with the alignment length.");

		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = super.getUserConfigParameters();
		params.add("MaxSymmOrder");
		params.add("UserOrder");
		params.add("SymmType");
		params.add("OrderDetectorMethod");
		params.add("RefineMethod");
		params.add("Optimization");
		params.add("RndSeed");
		params.add("SymmLevels");
		params.add("UnrefinedScoreThreshold");
		params.add("RefinedScoreThreshold");
		params.add("SSEThreshold");
		params.add("MinCoreLength");
		params.add("DistanceCutoff");
		params.add("Gaps");
		params.add("OptimizationSteps");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames() {
		List<String> params = super.getUserConfigParameterNames();
		params.add("Maximum Order of Symmetry");
		params.add("User Input Order");
		params.add("Type of Symmetry");
		params.add("Order Detection Method");
		params.add("Refinement Method");
		params.add("Optimization");
		params.add("Random Seed");
		params.add("Symmetry Levels");
		params.add("Unrefined Score Threshold");
		params.add("Refined Score Threshold");
		params.add("SSE Threshold");
		params.add("Minimum Core Length");
		params.add("Distance Cutoff");
		params.add("Internal Gaps");
		params.add("Optimization Steps");
		return params;
	}

	@Override
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = super.getUserConfigTypes();
		params.add(Integer.class);
		params.add(Integer.class);
		params.add(SymmetryType.class);
		params.add(OrderDetectorMethod.class);
		params.add(RefineMethod.class);
		params.add(Boolean.class);
		params.add(Integer.class);
		params.add(Integer.class);
		params.add(Double.class);
		params.add(Double.class);
		params.add(Integer.class);
		params.add(Integer.class);
		params.add(Double.class);
		params.add(Boolean.class);
		params.add(Integer.class);
		return params;
	}

	public RefineMethod getRefineMethod() {
		return refineMethod;
	}

	public void setRefineMethod(RefineMethod refineMethod) {
		this.refineMethod = refineMethod;
	}

	@Deprecated
	public void setRefineResult(boolean doRefine) {
		if (!doRefine) {
			refineMethod = RefineMethod.NOT_REFINED;
		} else {
			refineMethod = RefineMethod.DEFAULT;
		}
	}

	public OrderDetectorMethod getOrderDetectorMethod() {
		return orderDetectorMethod;
	}

	public void setOrderDetectorMethod(OrderDetectorMethod orderDetectorMethod) {
		this.orderDetectorMethod = orderDetectorMethod;
	}

	public void setUserOrder(Integer userOrder) {
		this.userOrder = userOrder;
	}

	public int getUserOrder() {
		return userOrder;
	}

	public void setMaxSymmOrder(Integer maxSymmOrder) {
		this.maxSymmOrder = maxSymmOrder;
	}

	public int getMaxSymmOrder() {
		return maxSymmOrder;
	}

	public SymmetryType getSymmType() {
		return symmType;
	}

	public void setSymmType(SymmetryType type) {
		this.symmType = type;
	}

	public boolean getOptimization() {
		return optimization;
	}

	public void setOptimization(Boolean optimization) {
		this.optimization = optimization;
	}

	public int getRndSeed() {
		return rndSeed;
	}

	public void setRndSeed(Integer seed) {
		this.rndSeed = seed;
	}

	public int getSymmLevels() {
		return symmLevels;
	}

	public void setSymmLevels(Integer symmLevels) {
		this.symmLevels = symmLevels;
	}

	public double getUnrefinedScoreThreshold() {
		return unrefinedScoreThreshold;
	}

	public void setUnrefinedScoreThreshold(Double unrefinedScoreThreshold) {
		this.unrefinedScoreThreshold = unrefinedScoreThreshold;
	}
	
	public double getRefinedScoreThreshold() {
		return refinedScoreThreshold;
	}

	public void setRefinedScoreThreshold(Double refinedScoreThreshold) {
		this.refinedScoreThreshold = refinedScoreThreshold;
	}

	public int getSSEThreshold() {
		return sseThreshold;
	}

	public void setSSEThreshold(Integer sseThreshold) {
		this.sseThreshold = sseThreshold;
	}

	public int getMinCoreLength() {
		return minCoreLength;
	}

	public void setMinCoreLength(Integer minCoreLength) {
		this.minCoreLength = minCoreLength;
	}

	public double getDistanceCutoff() {
		return distanceCutoff;
	}

	public void setDistanceCutoff(Double distanceCutoff) {
		this.distanceCutoff = distanceCutoff;
	}

	public boolean isGaps() {
		return gaps;
	}

	public void setGaps(Boolean gaps) {
		this.gaps = gaps;
	}

	public int getOptimizationSteps() {
		return optimizationSteps;
	}

	public void setOptimizationSteps(Integer optimizationSteps) {
		this.optimizationSteps = optimizationSteps;
	}

	@Override
	public String toString() {
		return "CESymmParameters [maxSymmOrder=" + maxSymmOrder
				+ ", userOrder=" + userOrder + ", symmType=" + symmType
				+ ", orderDetectorMethod=" + orderDetectorMethod
				+ ", refineMethod=" + refineMethod + ", optimization="
				+ optimization + ", rndSeed=" + rndSeed + ", symmLevels="
				+ symmLevels + ", unrefinedScoreThreshold="
				+ unrefinedScoreThreshold + ", refinedScoreThreshold="
				+ refinedScoreThreshold + ", sseThreshold=" + sseThreshold
				+ ", minCoreLength=" + minCoreLength + ", distanceCutoff="
				+ distanceCutoff + ", gaps=" + gaps + ", optimizationSteps="
				+ optimizationSteps + "]";
	}

}
