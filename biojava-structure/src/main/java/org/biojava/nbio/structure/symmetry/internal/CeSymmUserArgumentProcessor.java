package org.biojava.nbio.structure.symmetry.internal;

import org.biojava.nbio.structure.align.StructureAlignment;
import org.biojava.nbio.structure.align.ce.CeUserArgumentProcessor;
import org.biojava.nbio.structure.align.ce.StartupParameters;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.OrderDetectorMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.RefineMethod;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;

public class CeSymmUserArgumentProcessor extends CeUserArgumentProcessor{
	
	protected static class CeSymmStartupParams 
	extends CeUserArgumentProcessor.CeStartupParams {
		
		protected int maxSymmOrder;
		private SymmetryType symmetryType;
		protected OrderDetectorMethod orderDetectorMethod;
		protected RefineMethod refineMethod;
		private boolean optimization;
		private int seed;
		private boolean multipleAxes;
		private double symmetryThreshold;

		public CeSymmStartupParams() {
			super();
			maxSymmOrder = 8;
			symmetryType = SymmetryType.DEFAULT;
			orderDetectorMethod = OrderDetectorMethod.DEFAULT;
			refineMethod = RefineMethod.DEFAULT;
			optimization = true;
			seed = 0;
			multipleAxes = true;
			symmetryThreshold = CESymmParameters.DEFAULT_SYMMETRY_THRESHOLD;
		}

		public OrderDetectorMethod getOrderDetectorMethod() {
			return orderDetectorMethod;
		}
		
		public SymmetryType getSymmetryType() {
			return symmetryType;
		}

		public void setSymmetryType(SymmetryType symmetryType) {
			this.symmetryType = symmetryType;
		}

		public boolean getOptimization() {
			return optimization;
		}

		public void setOptimization(boolean optimization) {
			this.optimization = optimization;
		}

		public int getSeed() {
			return seed;
		}

		public void setSeed(int seed) {
			this.seed = seed;
		}

		public void setOrderDetectorMethod(OrderDetectorMethod method) {
			this.orderDetectorMethod = method;
		}
		
		public RefineMethod getRefineMethod() {
			return refineMethod;
		}
		
		public void setRefineMethod(RefineMethod refineMethod) {
			this.refineMethod = refineMethod;
		}

		public int getMaxSymmOrder() {
			return maxSymmOrder;
		}
		
		public void setMaxSymmOrder(int maxSymmOrder) {
			this.maxSymmOrder = maxSymmOrder;
		}

		public boolean isMultipleAxes() {
			return multipleAxes;
		}

		public void setMultipleAxes(boolean multipleAxes) {
			this.multipleAxes = multipleAxes;
		}

		public double getSymmetryThreshold() {
			return symmetryThreshold;
		}

		public void setSymmetryThreshold(double symmetryThreshold) {
			this.symmetryThreshold = symmetryThreshold;
		}

		@Override
		public String toString() {
			return "CeSymmStartupParams [maxSymmOrder=" + maxSymmOrder
					+ ", symmetryType=" + symmetryType
					+ ", orderDetectorMethod=" + orderDetectorMethod
					+ ", refineMethod=" + refineMethod + ", optimization="
					+ optimization + ", seed=" + seed + ", multipleAxes="
					+ multipleAxes + ", symmetryThreshold=" + symmetryThreshold
					+ ", maxGapSize=" + maxGapSize + ", winSize=" + winSize
					+ ", scoringStrategy=" + scoringStrategy + ", maxOptRMSD="
					+ maxOptRMSD + ", gapOpen=" + gapOpen + ", gapExtension="
					+ gapExtension + ", showAFPRanges=" + showAFPRanges + "]";
		}

	}
	
	@Override
	protected StartupParameters getStartupParametersInstance() {
		return new CeSymmStartupParams();
	}
	
	@Override
	public StructureAlignment getAlgorithm() {
		return new CeSymm();
	}
	

	@Override
	public Object getParameters() {
		
		StructureAlignment alignment = getAlgorithm();
		
		CESymmParameters aligParams = 
				(CESymmParameters) alignment.getParameters();
		CeSymmStartupParams startParams = (CeSymmStartupParams) params;
		
		if ( aligParams == null)
			aligParams = new CESymmParameters();
		
		// Copy relevant parameters from the startup parameters
		aligParams.setMaxGapSize(startParams.getMaxGapSize());
		aligParams.setWinSize(startParams.getWinSize());
		aligParams.setScoringStrategy(startParams.getScoringStrategy());
		aligParams.setMaxOptRMSD(startParams.getMaxOptRMSD());
		aligParams.setGapOpen(startParams.getGapOpen());
		aligParams.setGapExtension(startParams.getGapExtension());
		aligParams.setShowAFPRanges(startParams.isShowAFPRanges());
		aligParams.setMaxSymmOrder(startParams.getMaxSymmOrder());
		aligParams.setOrderDetectorMethod(startParams.getOrderDetectorMethod());
		aligParams.setRefineMethod(startParams.getRefineMethod());
		aligParams.setSymmetryType(startParams.getSymmetryType());
		aligParams.setOptimization(startParams.getOptimization());
		aligParams.setSeed(startParams.getSeed());
		aligParams.setMultipleAxes(startParams.isMultipleAxes());
		aligParams.setSymmetryThreshold(startParams.getSymmetryThreshold());
		
		return aligParams;
	}

	@Override
	public String getDbSearchLegend(){
		/*String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2"
				+ "\tsim1\tsim2\t " ;
		return legend;*/
		
		return "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2"
				+ "\tcov1\tcov2\t%ID\tDescription\t ";
		
	}
}
