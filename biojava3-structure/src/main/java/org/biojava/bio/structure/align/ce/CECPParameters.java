package org.biojava.bio.structure.align.ce;

import java.util.List;

public class CECPParameters extends CeParameters {

	public static final int DEFAULT_MIN_CP_LENGTH = 5; //The minimum block length for CPs. Blocks shorter than this will be ignored.

	public static enum DuplicationHint {
		DUPLICATE_SHORTER,
		DUPLICATE_LEFT,
		DUPLICATE_RIGHT,
	}
	
	protected DuplicationHint duplicationHint;
	protected int minCPLength;
	
	public CECPParameters() {
		super();
		duplicationHint = DuplicationHint.DUPLICATE_SHORTER;
		minCPLength = DEFAULT_MIN_CP_LENGTH;
		setMaxGapSize(0);
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
		duplicationHint = DuplicationHint.DUPLICATE_SHORTER;
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


	public int getMinCPLength() {
		return minCPLength;
	}


	public void setMinCPLength(int minCPLength) {
		this.minCPLength = minCPLength;
	}
}
