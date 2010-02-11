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
 * Created on Sep 15, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.ce;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import org.biojava.bio.structure.StructureTools;



/** Contains the parameters that can be sent to CE
 * 
 * @author Andreas Prlic
 *
 */
public class CeParameters implements ConfigStrucAligParams  {
	int winSize;
	double rmsdThr;
	double rmsdThrJoin;
	String[] alignmentAtoms;
	private int maxGapSize;

	boolean showAFPRanges;
	boolean checkCircular;

	public CeParameters(){
		reset();
	}



	@Override
	public String toString() {
		return "CeParameters [alignmentAtoms="
		+ Arrays.toString(alignmentAtoms) + ", maxGapSize="
		+ maxGapSize + ", rmsdThr=" + rmsdThr + ", rmsdThrJoin="
		+ rmsdThrJoin + ", winSize=" + winSize + "]";
	}



	public void reset(){
		winSize = 8;
		rmsdThr = 3.0;
		rmsdThrJoin = 4.0;
		alignmentAtoms = new String[]{StructureTools.caAtomName};
		maxGapSize = 30;
		showAFPRanges = false;
		checkCircular = false;
	}

	/** The window size to look at
	 * 
	 * @return
	 */
	public Integer getWinSize() {
		return winSize;
	}
	public void setWinSize(Integer winSize) {
		this.winSize = winSize;
	}

	/** RMSD Threshold
	 * 
	 * @return
	 */
	public Double getRmsdThr() {
		return rmsdThr;
	}
	public void setRmsdThr(Double rmsdThr) {
		this.rmsdThr = rmsdThr;
	}

	/** RMSD threshold for joining of AFPs
	 * 
	 * @return
	 */
	public Double getRmsdThrJoin() {
		return rmsdThrJoin;
	}
	public void setRmsdThrJoin(Double rmsdThrJoin) {
		this.rmsdThrJoin = rmsdThrJoin;
	}

	public String[] getAlignmentAtoms() {
		return alignmentAtoms;
	}

	public void setAlignmentAtoms(String[] alignmentAtoms) {
		this.alignmentAtoms = alignmentAtoms;
	}


	public void setMaxGapSize(Integer maxGapSize){
		this.maxGapSize = maxGapSize;
	}

	/** the Max gap size parameter G . default is 30, which was
	 * described to obtained empirically in the CE paper.
	 * the larger the max gap size, the longer the compute time,
	 * but in same cases drastically improved results. Set to -1 to get the 
	 * §
	 * @return
	 */
	public Integer getMaxGapSize() {
		return maxGapSize;
	}


	public List<String> getUserConfigHelp() {
		List<String> params =new ArrayList<String>();
		String helpMaxGap = "This parameter configures the maximum gap size G, that is applied during the AFP extension. The larger the value, the longer the calculation time can become, Default value is 30. Set to 0 for no limit. " ;
		String helpRmsdThr = "This configures the RMSD threshold applied during the trace of the fragment matrix.";
		String helpWinSize = "This configures the fragment size m of Aligned Fragment Pairs (AFPs).";
		params.add(helpMaxGap);
		params.add(helpRmsdThr);
		params.add(helpWinSize);
		return params;
	}

	public List<String> getUserConfigParameters() {
		List<String> params = new ArrayList<String>();
		params.add("MaxGapSize");
		params.add("RmsdThr");
		params.add("WinSize");
		return params;
	}

	public List<String> getUserConfigParameterNames(){
		List<String> params = new ArrayList<String>();
		params.add("max. gap size G (during AFP extension).");
		params.add("RMSD threshold during trace of the fragment matrix.");
		params.add("fragment size m");
		return params;
	}

	public List<Class> getUserConfigTypes() {
		List<Class> params = new ArrayList<Class>();
		params.add(Integer.class);
		params.add(Double.class);
		params.add(Integer.class);
		return params;
	}



	/**
	 * @return whether information about AFPs should be printed
	 */
	public boolean isShowAFPRanges()
	{
		return showAFPRanges;
	}
	public void setShowAFPRanges(boolean showAFPRanges)
	{
		this.showAFPRanges = showAFPRanges;
	}


	/**
	 * @return whether the protein should be checked for circular permutations.
	 */
	public boolean isCheckCircular() {
		return checkCircular;
	}
	/**
	 * @param checkCircular whether the protein should be checked for circular permutations
	 */
	public void setCheckCircular(boolean checkCircular) {
		this.checkCircular = checkCircular;
	}

}
