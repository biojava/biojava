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
package org.biojava.nbio.structure.symmetry.core;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Bean for a single sequence. These are intended to be unique sequences (100% id)
 * as an imput to clustering.
 *
 */
public class UniqueSequenceList implements Cloneable {
	private String sequenceString = "";
	private String seqResSequence = "";
    private List<Integer> alignment1 = null;
    private List<Integer> alignment2 = null;
    private Atom[] caAtoms = null;
    private String chainId = null;
    private int modelNumber = -1;
    private int structureId = -1;
    
    public UniqueSequenceList(Atom[] cAlphaAtoms, String chainId, int modelNumber, int structureId, String seqResSequence) {
    	this.caAtoms = cAlphaAtoms;
    	this.chainId = chainId;
    	this.modelNumber = modelNumber;
    	this.structureId = structureId;
    	this.seqResSequence = seqResSequence;
    	this.sequenceString =  getSequenceString(cAlphaAtoms);
    	this.alignment1 = new ArrayList<Integer>(cAlphaAtoms.length);
    	this.alignment2 = new ArrayList<Integer>(cAlphaAtoms.length);
    	for (int i = 0; i < cAlphaAtoms.length; i++) {
    		this.alignment1.add(i);
    		this.alignment2.add(i);
    	}
    }
    
    /**
     * Return true is the sequence and residues numbers of the passed in array of
     * atoms matches those of this unique sequence list
     * 
     * @param caAlphaAtoms
     * @return
     */
    public boolean isMatch(Atom[] caAlphaAtoms) {
    	return sequenceString.equals(getSequenceString(caAlphaAtoms));
    }
    
    public String getChainId() {
    	return chainId;
    }
    
    public int getModelNumber() {
    	return modelNumber;
    }
     
    public int getStructureId() {
    	return structureId;
    }
    
    public Atom[] getCalphaAtoms() {
    	return caAtoms;
    }
	
	public String getSeqResSequence() {
		return seqResSequence;
	}
	
	/**
	 * @param sequenceString the sequenceString to set
	 */
	public void setSequenceString(String sequenceString) {
		this.sequenceString = sequenceString;
	}
	/**
	 * @return the alignment1
	 */
	public List<Integer> getAlignment1() {
		return alignment1;
	}
	/**
	 * @param alignment1 the alignment1 to set
	 */
	public void setAlignment1(List<Integer> alignment1) {
		this.alignment1 = alignment1;
	}
	/**
	 * @return the alignment2
	 */
	public List<Integer> getAlignment2() {
		return alignment2;
	}
	/**
	 * @param alignment2 the alignment2 to set
	 */
	public void setAlignment2(List<Integer> alignment2) {
		this.alignment2 = alignment2;
	}
	
	@Override
	public Object clone() {
		UniqueSequenceList copy = null;
		try {
			copy = (UniqueSequenceList) super.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// deep copy lists and arrays
		copy.alignment1 = new ArrayList<Integer>(this.alignment1);
		copy.alignment2 = new ArrayList<Integer>(this.alignment2);
		copy.caAtoms = Arrays.copyOf(this.caAtoms, this.caAtoms.length); // note, that atoms in this array will be identical (this is intended)
		return copy;
	}
	
	public static String getSequenceString(Atom[] caAlphaAtoms) {
		StringBuilder builder = new StringBuilder();

		for (Atom a:  caAlphaAtoms) {
			Group g = a.getGroup();
			// TODO is the check for UNK required? UNK should have been filtered already in ChainClusterer?
			if (! g.getPDBName().equals("UNK")) {
				builder.append(g.getResidueNumber());
				builder.append(g.getPDBName());
			}
		}
		
//		System.out.println("getSequenceString: " + builder.toString());
		return builder.toString();
	}
     
	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("length: ");
		builder.append(caAtoms.length);
		builder.append(" seq: ");
		builder.append(sequenceString);
		builder.append("\n");
		builder.append(alignment1);
		builder.append("\n");
		builder.append(alignment2);
		return builder.toString();
	}
	
}
