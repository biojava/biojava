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
 * Created on DATE
 *
 */

package org.biojava3.core.sequence;

/**
 *
 * @author Scooter Willis
 */

import java.util.LinkedHashMap;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.loader.SequenceStringProxyLoader;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.SequenceProxyLoader;

public class DNASequence extends AbstractSequence<NucleotideCompound> {

    private DNASequence parentDNASequence = null;
    private Integer begin = null;
    private Integer end = null;
    private LinkedHashMap<String, GeneSequence> geneSequenceHashMap = new LinkedHashMap<String, GeneSequence>();

    public DNASequence() {
        throw new UnsupportedOperationException("Null constructor not supported");
    }

    public DNASequence(String seqString) {
        super(seqString, DNACompoundSet.getDNACompoundSet());
    }

    public DNASequence(SequenceProxyLoader<NucleotideCompound> proxyLoader) {
        super(proxyLoader, DNACompoundSet.getDNACompoundSet());
    }

    public void setParentDNASequence(DNASequence parentDNASequence) {
        this.parentDNASequence = parentDNASequence;
    }

    public DNASequence getParentDNASequence(){
        return parentDNASequence;
    }

    public GeneSequence removeGeneSequence(String accession) {
        return geneSequenceHashMap.remove(accession);
    }

    public GeneSequence addGene(AccessionID accession, int begin, int end) {
        GeneSequence geneSequence = new GeneSequence(this, begin, end);
        geneSequence.setAccession(accession);
        geneSequenceHashMap.put(accession.toString(), geneSequence);
        return geneSequence;
    }

    public GeneSequence getGene(String accession) {
        return geneSequenceHashMap.get(accession);
    }

    public DNASequence getReverseComplement() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public DNASequence getReverse() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public DNASequence getComplement() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public RNASequence getRNASequence() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * @return the begin
     */
    public int getBegin() {
        return begin;
    }

    /**
     * @param begin the begin to set
     */
    public void setBegin(Integer begin) {
        this.begin = begin;
    }

    /**
     * @return the end
     */
    public Integer getEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(Integer end) {
        this.end = end;
    }

    public static void main(String[] args) {
        DNASequence dnaSequence = new DNASequence("ATCG");
        System.out.println(dnaSequence.toString());

        SequenceStringProxyLoader<NucleotideCompound> sequenceStringProxyLoader = new SequenceStringProxyLoader("GCTA", DNACompoundSet.getDNACompoundSet());
        DNASequence dnaSequenceFromProxy = new DNASequence(sequenceStringProxyLoader);
        System.out.println(dnaSequenceFromProxy.toString());

    }
}
