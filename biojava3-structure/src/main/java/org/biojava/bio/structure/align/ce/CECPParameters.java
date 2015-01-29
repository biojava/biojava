package org.biojava.bio.structure.align.ce;

import java.util.List;

/**
 * Provides parameters to {@link CeCPMain}
 * 
 * @author Spencer Bliven
 *
 */
public class CECPParameters extends CeParameters {

	public static final int DEFAULT_MIN_CP_LENGTH = 5; //The minimum block length for CPs. Blocks shorter than this will be ignored.

	public static enum DuplicationHint {
		SHORTER("Shorter of the two"),
		LEFT("Left"),
		RIGHT("Right");
		
		
		private String name;
		private DuplicationHint(String name) {
			this.name = name;
		}
		@Override
		public String toString() {
			return name;
		}
	}
	
	protected DuplicationHint duplicationHint;
	protected Integer minCPLength;
	
	public CECPParameters() {
		super();
		// super calls reset();
	}

	@Override
	public String toString() {
		return "CECPParameters [scoringStrategy=" + scoringStrategy 
		+ ", maxGapSize=" + maxGapSize 
		+ ", rmsdThr=" + rmsdThr 
		+ ", rmsdThrJoin="+ rmsdThrJoin 
		+ ", winSize=" + winSize 
		+ ", showAFPRanges=" + showAFPRanges 
		+ ", maxOptRMSD=" + maxOptRMSD
		+ ", seqWeight=" + seqWeight
		+ ", duplicationHint=" + duplicationHint
		+ ", minCPLength=" + minCPLength
		+ "]";
	}


	@Override
	public void reset(){
		super.reset();
		duplicationHint = DuplicationHint.SHORTER;
		minCPLength = DEFAULT_MIN_CP_LENGTH;
		setMaxGapSize(0);
	}


	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = super.getUserConfigHelp();
		params.add("Direction to duplicate: SHORTER, LEFT, or RIGHT");
		params.add("Minimum length of a CP block to consider");
		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = super.getUserConfigParameters();
		params.add("DuplicationHint");
		params.add("MinCPLength");
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames(){
		List<String> params = super.getUserConfigParameterNames();
		
		params.add("Which to duplicate");
		params.add("Min CP Length");
		return params;
	}

	@Override
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = super.getUserConfigTypes();
		params.add(DuplicationHint.class);
		params.add(Integer.class);
		return params;
	}
	
	public DuplicationHint getDuplicationHint() {
		return duplicationHint;
	}

	public void setDuplicationHint(DuplicationHint duplicationHint) {
		this.duplicationHint = duplicationHint;
	}


	public Integer getMinCPLength() {
		return minCPLength;
	}


	public void setMinCPLength(Integer minCPLength) {
		this.minCPLength = minCPLength;
	}
}
