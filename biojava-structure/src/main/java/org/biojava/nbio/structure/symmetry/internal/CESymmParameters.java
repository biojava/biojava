package org.biojava.nbio.structure.symmetry.internal;

import java.util.List;
import java.util.Random;

import org.biojava.nbio.structure.align.ce.CeParameters;

/**
 * Provides parameters to {@link CeSymm}.
 * 
 * @author Spencer Bliven
 * @author Aleix Lafita
 * @since 4.2.0
 * 
 */
public class CESymmParameters extends CeParameters {

	private int maxSymmOrder;
	private SymmetryType symmetryType;
	private OrderDetectorMethod orderDetectorMethod;
	private int userOrder;
	private RefineMethod refineMethod;
	private boolean optimization;
	private int seed; //random seed
	private boolean multipleAxes;
	private double symmetryThreshold;
	private int minSubunitLength;
	private double distanceCutoff;
	private boolean gaps;

	public static enum OrderDetectorMethod {
		SEQUENCE_FUNCTION,
		USER_INPUT;
		public static final OrderDetectorMethod DEFAULT = SEQUENCE_FUNCTION;
	}

	public static enum RefineMethod {
		NOT_REFINED,
		SINGLE,
		MULTIPLE;
		public static final RefineMethod DEFAULT = SINGLE;
	}

	public static final double DEFAULT_SYMMETRY_THRESHOLD = 0.4;

	/**
	 * The internal symmetry detection can be divided into two types: 
	 * CLOSE: includes the circular and dihedral symmetries, and
	 * OPEN: includes the helical and protein repeats symmetries.
	 * <p>
	 * All internal symmetry cases share one property: 
	 * all the subunits have the same 3D transformation.
	 * <p>
	 * AUTO option automatically identifies the type. 
	 * The criterion for classification is that the CLOSE symmetry 
	 * generates CeSymm alignments with circular permutations 
	 * (2 blocks in AFPChain), whereas the OPEN symmetry
	 * generates alignments without a CP (only one block in AFPChain).
	 */
	public enum SymmetryType {
		CLOSE,
		OPEN,
		AUTO;
		public static final SymmetryType DEFAULT = AUTO;
	}

	public CESymmParameters() {
		reset();
	}
	
	@Override
	public CESymmParameters clone(){
		
		CESymmParameters p = new CESymmParameters();
		
		p.maxSymmOrder = maxSymmOrder;
		p.symmetryType = symmetryType;
		p.orderDetectorMethod = orderDetectorMethod;
		p.userOrder = userOrder;
		p.refineMethod = refineMethod;
		p.optimization = optimization;
		p.seed = seed;
		p.multipleAxes = multipleAxes;
		p.symmetryThreshold = symmetryThreshold;
		p.minSubunitLength = minSubunitLength;
		p.distanceCutoff = distanceCutoff;
		p.gaps = gaps;
		
		p.winSize = winSize;
		p.rmsdThr = rmsdThr;
		p.rmsdThrJoin = rmsdThrJoin;
		p.scoringStrategy = scoringStrategy;
		p.maxGapSize = maxGapSize;
		p.showAFPRanges = showAFPRanges;
		p.maxOptRMSD = maxOptRMSD;
		p.gapOpen = gapOpen;
		p.gapExtension = gapExtension;
		p.distanceIncrement = distanceIncrement;
		p.oRmsdThr = oRmsdThr;
		p.maxNrIterationsForOptimization = maxNrIterationsForOptimization;
		p.seqWeight = seqWeight;
		
		return p;
	}

	@Override
	public void reset(){
		super.reset();
		maxSymmOrder = 8;
		symmetryType = SymmetryType.DEFAULT;
		orderDetectorMethod = OrderDetectorMethod.DEFAULT;
		userOrder = 0;
		refineMethod = RefineMethod.DEFAULT;
		optimization = true;
		seed = new Random().nextInt(10000);
		multipleAxes = true;
		symmetryThreshold = DEFAULT_SYMMETRY_THRESHOLD;
		minSubunitLength = 15;
		distanceCutoff = 7.0;
		gaps = true;
	}

	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = super.getUserConfigHelp();

		//maxSymmOrder help explanation
		params.add("Sets the maximum order of symmetry of the protein.");

		StringBuilder symmTypes = new StringBuilder("Type of Symmetry: ");
		SymmetryType[] vals = SymmetryType.values();
		if(vals.length == 1) {
			symmTypes.append(vals[0].name());
		} else if(vals.length > 1 ) {
			for(int i=0;i<vals.length-1;i++) {
				symmTypes.append(vals[i].name());
				symmTypes.append(", ");
			}
			symmTypes.append("or ");
			symmTypes.append(vals[vals.length-1].name());
		}
		params.add(symmTypes.toString());

		StringBuilder orderTypes = new StringBuilder("Order Detection Method: ");
		OrderDetectorMethod[] vals2 = OrderDetectorMethod.values();
		if(vals2.length == 1) {
			orderTypes.append(vals2[0].name());
		} else if(vals2.length > 1 ) {
			for(int i=0;i<vals2.length-1;i++) {
				orderTypes.append(vals2[i].name());
				orderTypes.append(", ");
			}
			orderTypes.append("or ");
			orderTypes.append(vals[vals.length-1].name());
		}
		params.add(orderTypes.toString());
		
		//userOrder help explanation
		params.add("Order of symmetry determined by the user. "
				+ "Use it with the USER_INPUT order option. Imposes an order"
				+ " of symmetry to the alignment. If 0 the order is set "
				+ "automatically.");

		StringBuilder refineTypes = new StringBuilder("Refinement Method: ");
		RefineMethod[] values = RefineMethod.values();
		if(values.length == 1) {
			refineTypes.append(values[0].name());
		} else if(values.length > 1 ) {
			for(int i=0;i<values.length-1;i++) {
				refineTypes.append(values[i].name());
				refineTypes.append(", ");
			}
			refineTypes.append("or ");
			refineTypes.append(values[values.length-1].name());
		}
		params.add(refineTypes.toString());

		//optimization help explanation
		params.add("Optimize the refined alignment if true.");
		//seed help explanation
		params.add("Random seed for the Monte Carlo optimization, "
				+ "for reproducibility of results.");
		//multiple axes explanation
		params.add("Perform multiple CeSymm iterations to find all "
				+ "symmetry axes if true.");
		//threshold
		params.add("Symmetry threshold: TM-score values below the "
				+ "threshold will be considered asymmetric.");
		//min subunit length
		params.add("Subunit length: the minimum number of non-gapped "
				+ "residues in every symmetric subunit.");
		//distance cutoff
		params.add("Distance Cutoff: the maximum allowed distance (in A) "
				+ "between two aligned residues.");
		
		//gaps
		params.add("Gaps: allow gaps in the alignment if true.");

		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = super.getUserConfigParameters();
		params.add("MaxSymmOrder");
		params.add("SymmetryType");
		params.add("OrderDetectorMethod");
		params.add("UserOrder");
		params.add("RefineMethod");
		params.add("Optimization");
		params.add("Seed");
		params.add("MultipleAxes");
		params.add("SymmetryThreshold");
		params.add("MinSubunitLength");
		params.add("DistanceCutoff");
		params.add("Gaps");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames(){
		List<String> params = super.getUserConfigParameterNames();
		params.add("Maximum Order of Symmetry");
		params.add("Type of Symmetry");
		params.add("Order Detection Method");
		params.add("User Input Order");
		params.add("Refinement Method");
		params.add("Optimization");
		params.add("Random Seed");
		params.add("Multiple Axes");
		params.add("Symmetry Threshold");
		params.add("Minimum Subunit Length");
		params.add("Distance Cutoff");
		params.add("Gaps");
		return params;
	}

	@Override
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = super.getUserConfigTypes();
		params.add(Integer.class);
		params.add(SymmetryType.class);
		params.add(OrderDetectorMethod.class);
		params.add(Integer.class);
		params.add(RefineMethod.class);
		params.add(Boolean.class);
		params.add(Integer.class);
		params.add(Boolean.class);
		params.add(Double.class);
		params.add(Integer.class);
		params.add(Double.class);
		params.add(Boolean.class);
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
		if (!doRefine){
			refineMethod = RefineMethod.NOT_REFINED;
		}
		else{
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

	public SymmetryType getSymmetryType() {
		return symmetryType;
	}

	public void setSymmetryType(SymmetryType type) {
		this.symmetryType = type;
	}

	public boolean getOptimization() {
		return optimization;
	}

	public void setOptimization(Boolean optimization) {
		this.optimization = optimization;
	}

	public int getSeed() {
		return seed;
	}

	public void setSeed(Integer seed) {
		this.seed = seed;
	}

	public boolean isMultipleAxes() {
		return multipleAxes;
	}

	public void setMultipleAxes(Boolean multipleAxes) {
		this.multipleAxes = multipleAxes;
	}

	public double getSymmetryThreshold() {
		return symmetryThreshold;
	}

	public void setSymmetryThreshold(Double symmetryThreshold) {
		this.symmetryThreshold = symmetryThreshold;
	}

	public int getMinSubunitLength() {
		return minSubunitLength;
	}

	public void setMinSubunitLength(Integer minSubunitLength) {
		this.minSubunitLength = minSubunitLength;
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

}
