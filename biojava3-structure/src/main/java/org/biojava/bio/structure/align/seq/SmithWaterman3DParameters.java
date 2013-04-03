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

package org.biojava.bio.structure.align.seq;

import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.structure.align.ce.ConfigStrucAligParams;

public class SmithWaterman3DParameters implements ConfigStrucAligParams
{


	short gapOpen ;      
	short gapExtend ;    

	public SmithWaterman3DParameters(){
		reset();
	}



	public List<String> getUserConfigHelp()
	{
		List<String> params =new ArrayList<String>();
		params.add("The Gap open penalty");
		params.add("The Gap Extension penalty");



		// TODO Auto-generated method stub
		return params;
	}

	public List<String> getUserConfigParameterNames()
	{
		List<String> params =new ArrayList<String>();
		params.add("Gap Open");
		params.add("Gap Extension");


		return params;
	}

	public List<String> getUserConfigParameters()
	{
		List<String> params =new ArrayList<String>();
		params.add("GapOpen");
		params.add("GapExtend");      

		return params;
	}

	@SuppressWarnings("rawtypes")
	public List<Class> getUserConfigTypes()
	{
		List<Class> params = new ArrayList<Class>();
		params.add(Short.class);
		params.add(Short.class);

		return params;
	}

	public void reset()
	{
		gapOpen = (short) 8;
		gapExtend = (short) 1;      

	}


	public Short getGapExtend()
	{
		return gapExtend;
	}



	public void setGapExtend(Short gapExtend)
	{
		this.gapExtend = gapExtend;
	}



	public Short getGapOpen() {
		return gapOpen;
	}



	public void setGapOpen(Short gapOpen) {
		this.gapOpen = gapOpen;
	}




}
