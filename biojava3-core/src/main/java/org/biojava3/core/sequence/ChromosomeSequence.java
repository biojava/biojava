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

import java.util.LinkedHashMap;
import org.biojava3.core.sequence.compound.DNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;

/**
 *
 * @author Scooter Willis
 */
public class ChromosomeSequence extends DNASequence {

    private int chromosomeNumber;
    private LinkedHashMap<String, GeneSequence> geneSequenceHashMap = new LinkedHashMap<String, GeneSequence>();


        public ChromosomeSequence() {
//        throw new UnsupportedOperationException("Null constructor not supported");
    }

    public ChromosomeSequence(String seqString) {
        super(seqString, DNACompoundSet.getDNACompoundSet());
    }

    public ChromosomeSequence(ProxySequenceReader<NucleotideCompound> proxyLoader) {
        super(proxyLoader, DNACompoundSet.getDNACompoundSet());
    }

    public ChromosomeSequence(String seqString, CompoundSet<NucleotideCompound> compoundSet) {
        super(seqString, compoundSet);
    }

    public ChromosomeSequence(ProxySequenceReader<NucleotideCompound> proxyLoader, CompoundSet<NucleotideCompound> compoundSet) {
        super(proxyLoader, compoundSet);
    }

    /**
     * @return the chromosomeNumber
     */
    public int getChromosomeNumber() {
        return chromosomeNumber;
    }

    /**
     * @param chromosomeNumber the chromosomeNumber to set
     */
    public void setChromosomeNumber(int chromosomeNumber) {
        this.chromosomeNumber = chromosomeNumber;
    }

    public LinkedHashMap<String, GeneSequence> getGeneSequences() {
        return geneSequenceHashMap;
    }

    public GeneSequence removeGeneSequence(String accession) {
        return geneSequenceHashMap.remove(accession);
    }

    public GeneSequence addGene(AccessionID accession, int begin, int end, Strand strand) {
        GeneSequence geneSequence = new GeneSequence(this, begin, end, strand);
        geneSequence.setAccession(accession);
        geneSequenceHashMap.put(accession.toString(), geneSequence);
        return geneSequence;
    }

    public GeneSequence getGene(String accession) {
        return geneSequenceHashMap.get(accession);
    }
}
