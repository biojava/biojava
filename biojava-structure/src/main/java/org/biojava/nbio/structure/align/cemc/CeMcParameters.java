package org.biojava.nbio.structure.align.cemc;

import org.biojava.nbio.structure.align.ce.CeParameters;

/** 
 * Contains the parameters that can be sent to the CEMC algorithm.
 * It extends the {@link CeParameters} because it uses the CE algorithm internally.
 * 
 * @author Aleix Lafita
 *
 */
public class CeMcParameters extends CeParameters {

	public CeMcParameters(){
		super();
	}

	@Override
	public String toString() {
		return "CeMcParameters [winSize=" + winSize + ", rmsdThr=" + rmsdThr
				+ ", rmsdThrJoin=" + rmsdThrJoin + ", maxOptRMSD=" + maxOptRMSD
				+ ", scoringStrategy=" + scoringStrategy + ", maxGapSize="
				+ maxGapSize + ", showAFPRanges=" + showAFPRanges
				+ ", sideChainScoringType=" + sideChainScoringType
				+ ", gapOpen=" + gapOpen + ", gapExtension=" + gapExtension
				+ ", distanceIncrement=" + distanceIncrement + ", oRmsdThr="
				+ oRmsdThr + ", maxNrIterationsForOptimization="
				+ maxNrIterationsForOptimization + ", substitutionMatrix="
				+ substitutionMatrix + ", seqWeight=" + seqWeight + "]";
	}
}
