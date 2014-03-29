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
 * created at 27 Mar 2014
 * Author: ap3 
 */

package org.biojava.bio.structure.xtal.io;

import javax.xml.bind.annotation.XmlElement;

import org.biojava.bio.structure.xtal.SpaceGroup;

public class SpaceGroupMapElements {

	@XmlElement 
	public Integer key;
	
	@XmlElement(name="SpaceGroup", namespace="http://www.biojava.org")
	public SpaceGroup value;
	
	private SpaceGroupMapElements(){
		
	}
	
	public SpaceGroupMapElements(Integer key, SpaceGroup value){
		this.key = key;
		this.value = value;
	}
}
