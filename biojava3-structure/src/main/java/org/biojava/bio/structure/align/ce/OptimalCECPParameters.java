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

import java.util.List;

/** Contains the parameters that can be sent to CE
 * 
 * @author Andreas Prlic
 *
 */
public class OptimalCECPParameters extends CeParameters {
	/**
	 * If true, ignores {@link #cpPoint} and tries all possible cp points.
	 */
	protected Boolean tryAllCPs;
	/**
	 * The CP point, specified as a residue index
	 * 
	 * <p>TODO make this a ResidueNumber 
	 */
	protected Integer cpPoint;

	protected OptimalCECPParameters params;

	@Override
	public String toString() {
		return "OptimalCECPParameters [scoringStrategy=" + scoringStrategy 
		+ ", maxGapSize=" + maxGapSize 
		+ ", rmsdThr=" + rmsdThr 
		+ ", rmsdThrJoin="+ rmsdThrJoin 
		+ ", winSize=" + winSize 
		+ ", showAFPRanges=" + showAFPRanges 
		+ ", maxOptRMSD=" + maxOptRMSD
		+ ", seqWeight=" + seqWeight
		+ ", tryAllCPs" + tryAllCPs
		+ ", cpPoint" + cpPoint
		+ "]";
	}


	@Override
	public void reset(){
		super.reset();
		tryAllCPs = true;
		cpPoint = 0;
	}


	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = super.getUserConfigHelp();
		params.add("Should we try all CP sites? Otherwise, only try the site specified by CPPoint.");
		params.add("Index of the CP site. Ignored unless TryAllCPs=false");
		return params;
	}

	public List<String> getUserConfigParameters() {
		List<String> params = super.getUserConfigParameters();
		params.add("TryAllCPs");
		params.add("CPPoint");
		return params;
	}

	public List<String> getUserConfigParameterNames(){
		List<String> params = super.getUserConfigParameterNames();
		
		params.add("Try all CPs");
		params.add("CP Point");
		return params;
	}

	@SuppressWarnings("unchecked")
	public List<Class> getUserConfigTypes() {
		List<Class> params = super.getUserConfigTypes();
		params.add(Boolean.class);
		params.add(Integer.class);
		return params;
	}


	/**
	 * @return Whether we should try all CP sites
	 */
	public Boolean isTryAllCPs() {
		return tryAllCPs;
	}


	/**
	 * @param tryAllCPs Set whether we should try all CP sites
	 */
	public void setTryAllCPs(Boolean tryAllCPs) {
		this.tryAllCPs = tryAllCPs;
	}


	/**
	 * @return the cpPoint
	 */
	public Integer getCPPoint() {
		return cpPoint;
	}

	/**
	 * @param cpPoint the cpPoint to set
	 */
	public void setCPPoint(Integer cpPoint) {
		this.cpPoint = cpPoint;
	}



}
