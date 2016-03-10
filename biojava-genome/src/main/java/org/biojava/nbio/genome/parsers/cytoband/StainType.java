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
 * created at 20 Feb 2014
 * Author: ap3
 */

package org.biojava.nbio.genome.parsers.cytoband;

import java.io.Serializable;



public enum StainType implements Serializable{

	acen("acen"),
	gneg("gneg"),
	gpos100("gpos100"),
	gpos25("gpos25"),
	gpos50("gpos50"),
	gpos75("gpos75"),
	gvar("gvar"),
	stalk("stalk");

	public final String type;

	StainType(String type){
		this.type = type;
	}

	public static StainType getStainTypeFromString(String type){
		 for(StainType st : StainType.values()){
			 if (st.type.equals(type)){
				 return st;
			 }
		 }
		 return null;
	}
}
