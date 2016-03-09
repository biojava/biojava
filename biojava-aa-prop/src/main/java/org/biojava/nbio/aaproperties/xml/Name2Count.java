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
package org.biojava.nbio.aaproperties.xml;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;

@XmlAccessorType(XmlAccessType.FIELD)
public class Name2Count{
	@XmlAttribute(name = "name", required = true)
	private String name;
	@XmlAttribute(name = "count", required = true)
	private int count;

	public Name2Count(){}

	public Name2Count(String n, int c){
		if(c <= 0){
			throw new Error("Count must be > 0.");
		}
		this.name = n;
		this.count = c;
	}

	public String getName(){
		return this.name;
	}

	public int getCount(){
		return this.count;
	}
}
