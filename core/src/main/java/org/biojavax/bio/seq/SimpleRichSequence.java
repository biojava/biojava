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

package org.biojavax.bio.seq;

import org.biojava.bio.BioException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojavax.Namespace;

/**
 * A simple implementation of RichSequence.
 * @author Richard Holland
 * @author Bubba Puryear
 * @since 1.5
 */
public class SimpleRichSequence extends ThinRichSequence {
    
    private SymbolList symList;
        
    /**
     * Creates a new instance of SimpleRichSequence. Note the use of Double for
     * seqversion, which indicates that it is nullable.
     * @param ns the namespace for this sequence.
     * @param name the name of the sequence.
     * @param accession the accession of the sequence.
     * @param version the version of the sequence.
     * @param symList the symbols for the sequence.
     * @param seqversion the version of the symbols for the sequence.
     */
    public SimpleRichSequence(Namespace ns, String name, String accession, int version, SymbolList symList, Double seqversion) {
        super(ns,name,accession,version,symList.getAlphabet(),seqversion);
        this.symList = symList;
    }

    // Hibernate requirement - not for public use.
    protected SimpleRichSequence() {} 
                    
    // Hibernate requirement - not for public use.
    protected void setAlphabetName(String alphaname) throws IllegalSymbolException, BioException {
        super.setAlphabetName(alphaname);
        this.checkMakeSequence();
    }
        
    // Hibernate requirement - not for public use.
    private String seqstring;
    
    // Hibernate requirement - not for public use.
    protected void setStringSequence(String seq) throws IllegalSymbolException, BioException {
        this.seqstring = seq;
        // convert the string into a symbollist
        this.checkMakeSequence();
    }
    
    // Hibernate requirement - not for public use.
    protected String getStringSequence() { 
        // convert the symbollist into a string
        return (this.symList==SymbolList.EMPTY_LIST?null:this.seqString()); 
    }
    
    // when both alphabet and sequence have been received, make a symbollist out of them
    private void checkMakeSequence() throws IllegalSymbolException, BioException {
        if (this.getAlphabet()!=null && this.seqstring!=null) {
            // Make the symbol list and assign it.
            this.symList = new SimpleSymbolList(this.getAlphabet().getTokenization("token"), seqstring);
        }
    }    
    
    /**
     * {@inheritDoc}
     */
    public int length() { return this.symList.length(); }
        
    // Hibernate requirement - not for public use.
    protected void setSequenceLength(int length) {} // ignore this calculated field
    
    // Hibernate requirement - not for public use.
    protected int getSequenceLength() { return this.length(); }
    
    /**
     * {@inheritDoc}
     */
    public SymbolList getInternalSymbolList() {
        return this.symList;
    }
}
