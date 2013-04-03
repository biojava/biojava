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
 * Created on Nov 2, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;


import org.biojava.bio.structure.align.StructureAlignment;

/** process the arguments from command line
 * 
 * @author Andreas Prlic
 *
 */
public class CeUserArgumentProcessor extends AbstractUserArgumentProcessor {

	public StructureAlignment getAlgorithm() {
		return new CeMain();
	}


	@Override
	public Object getParameters() {
		
		StructureAlignment alignment = getAlgorithm();
		
		CeParameters p = (CeParameters) alignment.getParameters();
		
		if ( p == null)
			p = new CeParameters();
		
		p.setMaxOptRMSD(params.getMaxOptRMSD());
		p.setMaxGapSize(params.getMaxGapSize());
		p.setShowAFPRanges(params.isShowAFPRanges());
		return p;
	}


	public String getDbSearchLegend(){
		//String legend = "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tsim1\tsim2\t " ;
		//return legend;
		
		return "# name1\tname2\tscore\tz-score\trmsd\tlen1\tlen2\tcov1\tcov2\t%ID\tDescription\t " ;
		
		
	}
	

}
