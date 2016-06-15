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
package org.biojava.nbio.structure.align.ce;

import org.biojava.nbio.structure.align.StructureAlignment;

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
