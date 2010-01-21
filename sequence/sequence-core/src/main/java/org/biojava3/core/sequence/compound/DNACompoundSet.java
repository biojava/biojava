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
 * Created on 01-21-2010
 *
 * @author Richard Holland
 * @auther Scooter Willis
 *
 */
package org.biojava3.core.sequence.compound;

import java.util.HashMap;
import java.util.Map;

import org.biojava3.core.sequence.template.CompoundSet;


public class DNACompoundSet implements CompoundSet<NucleotideCompound>{

	private final Map<String,NucleotideCompound> dnaCompoundCache = new HashMap<String,NucleotideCompound>();
	
	public DNACompoundSet() {
		dnaCompoundCache.put("A", new NucleotideCompound("A",this,"T"));
		dnaCompoundCache.put("T", new NucleotideCompound("T",this,"A"));
		dnaCompoundCache.put("C", new NucleotideCompound("C",this,"G"));
		dnaCompoundCache.put("G", new NucleotideCompound("G",this,"C"));
		dnaCompoundCache.put("a", new NucleotideCompound("a",this,"t"));
		dnaCompoundCache.put("t", new NucleotideCompound("t",this,"a"));
		dnaCompoundCache.put("c", new NucleotideCompound("c",this,"g"));
		dnaCompoundCache.put("g", new NucleotideCompound("g",this,"c"));
	}
	
	public String getStringForCompound(NucleotideCompound compound) {
		return compound.toString();
	}

	public NucleotideCompound getCompoundForString(String string) {
		if (string.length()==0) {
			return null;
		}
		if (string.length()>this.getMaxSingleCompoundStringLength()) {
			throw new IllegalArgumentException("String supplied is too long.");
		}
		return this.dnaCompoundCache.get(string);
	}

	public int getMaxSingleCompoundStringLength() {
		return 1;
	}

 
        private final static DNACompoundSet dnaCompoundSet = new DNACompoundSet();

        static public DNACompoundSet getDNACompoundSet(){
            return dnaCompoundSet;
        }

}
