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
 * Created on July 27, 2010
 *
 */

package org.biojava.bio.structure;

import java.io.Serializable;

/**
 * ElementType is an enumeration of the types of elements found in the periodic table.
 * Each element type is further classified into Metal, Metalloid, and Non-Metal.
 * 
 * Element types based on definition at http://www.ptable.com/
 * 
 * @author Peter Rose
 * @version %I% %G%
 * @since 3.0
 */

public enum ElementType implements Serializable {
	
	METALLOID(false),
	OTHER_NONMETAL(false),
	HALOGEN(false),
	NOBLE_GAS(false),
	ALKALI_METAL(true),
	ALKALINE_EARTH_METAL(true),
	LANTHANOID(true),
	ACTINOID(true),
	TRANSITION_METAL(true),
	POST_TRANSITION_METAL(true),
	UNKNOWN(false);
	
	private boolean metal;
	
	private ElementType(boolean metal) {
		this.metal = metal;
	}
	
	/**
     * Returns <CODE>true</CODE> if ElementType is a metal.
     * @return <CODE>true</CODE> if ElementType is a metal.
     */
	public boolean isMetal() {
		return metal;
	}
	
	/**
     * Returns <CODE>true</CODE> if ElementType is a metalloid.
     * @return <CODE>true</CODE> if ElementType is a metalloid.
     */
	public boolean isMetalloid() {
		return this.equals(METALLOID);
	}
	
	/**
     * Returns <CODE>true</CODE> if ElementType is a non-metal.
     * @return <CODE>true</CODE> if ElementType is a non-metal.
     */
	public boolean isNonMetal() {
		return this.equals(OTHER_NONMETAL) || this.equals(HALOGEN) || this.equals(NOBLE_GAS);
	}
}
