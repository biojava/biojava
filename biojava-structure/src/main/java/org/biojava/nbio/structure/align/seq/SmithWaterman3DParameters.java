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
 * Created on Jun 7, 2010
 * Author: ap3
 *
 */

package org.biojava.nbio.structure.align.seq;

import org.biojava.nbio.structure.align.ce.ConfigStrucAligParams;

import java.util.ArrayList;
import java.util.List;

public class SmithWaterman3DParameters implements ConfigStrucAligParams {

	private short gapOpen;    // gap opening penalty for sequence alignment
	private short gapExtend;  // gap extension penalty for sequence alignment
	private double maxRmsd;   // maximum RMSD of superposition allowed
	private int minLen;       // minimum alignment length allowed

	public SmithWaterman3DParameters() {
		reset();
	}

	@Override
	public List<String> getUserConfigHelp() {
		List<String> params = new ArrayList<String>();
		params.add("The Gap open penalty");
		params.add("The Gap extension penalty");
		params.add("The maximum RMSD of superposition allowed");
		params.add("The minimum alignment length allowed");
		
		// TODO Auto-generated method stub
		return params;
	}

	@Override
	public List<String> getUserConfigParameterNames() {
		List<String> params = new ArrayList<String>();
		params.add("Gap Open");
		params.add("Gap Extension");
		params.add("Maximum RMSD");
		params.add("Minimum Alignment Length");
		
		return params;
	}

	@Override
	public List<String> getUserConfigParameters() {
		List<String> params = new ArrayList<String>();
		params.add("GapOpen");
		params.add("GapExtend");
		params.add("MaxRmsd");
		params.add("MinLen");
		
		return params;
	}

	@Override
	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes() {
		List<Class> params = new ArrayList<Class>();
		params.add(Short.class);
		params.add(Short.class);
		params.add(Double.class);
		params.add(Integer.class);

		return params;
	}

	@Override
	public void reset() {
		gapOpen = (short) 8;
		gapExtend = (short) 1;
		maxRmsd = 99;
		minLen = 30;

	}

	public Short getGapExtend() {
		return gapExtend;
	}

	public void setGapExtend(Short gapExtend) {
		this.gapExtend = gapExtend;
	}

	public Short getGapOpen() {
		return gapOpen;
	}

	public void setGapOpen(Short gapOpen) {
		this.gapOpen = gapOpen;
	}

	public Double getMaxRmsd() {
		return maxRmsd;
	}

	public void setMaxRmsd(Double maxRmsd) {
		this.maxRmsd = maxRmsd;
	}

	public Integer getMinLen() {
		return minLen;
	}

	public void setMinLen(Integer minLen) {
		this.minLen = minLen;
	}

}
