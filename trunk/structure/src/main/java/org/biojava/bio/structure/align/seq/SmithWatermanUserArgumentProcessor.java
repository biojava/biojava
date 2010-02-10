package org.biojava.bio.structure.align.seq;


import org.biojava.bio.structure.align.StructureAlignment;
import org.biojava.bio.structure.align.ce.AbstractUserArgumentProcessor;


public class SmithWatermanUserArgumentProcessor extends AbstractUserArgumentProcessor{

	

	
	public StructureAlignment getAlgorithm() {
		return new SmithWaterman3Daligner();
	}

	

	@Override
	public Object getParameters() {
		return null;
	}
	
	public String getDbSearchLegend(){
		String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		return legend;
	}

	
}
