package org.biojava.bio.structure.align.ce;

import org.biojava.bio.structure.align.StructureAlignment;

public class CeSideChainUserArgumentProcessor extends CeUserArgumentProcessor {

	@Override
	public StructureAlignment getAlgorithm() {
		return new CeSideChainMain();
	}

	@Override
	public Object getParameters() {
	   
		CeParameters params = new CeParameters();
	
		params.setScoringStrategy(CeParameters.ScoringStrategy.SIDE_CHAIN_SCORING);
		//params.setMaxGapSize(0);
		return params;
	}


	@Override
	public String getDbSearchLegend(){
		String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		return legend;
	}

}
