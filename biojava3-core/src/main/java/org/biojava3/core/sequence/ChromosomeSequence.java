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
import org.biojava3.core.sequence.template.SequenceReader;

/**
 * A ChromosomeSequence is a DNASequence but keeps track of geneSequences
 * @author Scooter Willis
 */
public class ChromosomeSequence extends DNASequence {

    private int chromosomeNumber;
    private LinkedHashMap<String, GeneSequence> geneSequenceHashMap = new LinkedHashMap<String, GeneSequence>();

    /**
     * Empty constructor used by tools that need a proper Bean that allows the actual
     * sequence data to be set after construction. Not recommended
     */
    public ChromosomeSequence() {
//        throw new UnsupportedOperationException("Null constructor not supported");
    }

    /**
     * String is king and assume DNA
     * @param seqString
     */
    public ChromosomeSequence(String seqString) {
        super(seqString, DNACompoundSet.getDNACompoundSet());
    }

    /**
     * Fairly important constructor given the size of a ChromsomeSequence where the
     * ProxySequenceReader could load from disk via RandomAccessFile so that the sequence
     * doesn't need to be kept in memory. Could also be a NCBI proxy to load sequence
     * data as needed from remote web server.
     * @param proxyLoader
     */
    public ChromosomeSequence(SequenceReader<NucleotideCompound> proxyLoader) {
        super(proxyLoader, DNACompoundSet.getDNACompoundSet());
    }

    /**
     * Allows the creation of a ChromosomeSequence using String for the sequence with a custom CompoundSet
     * @param seqString
     * @param compoundSet
     */
    public ChromosomeSequence(String seqString, CompoundSet<NucleotideCompound> compoundSet) {
        super(seqString, compoundSet);
    }

    /**
     * Allows the creation of a ChromosomeSequence using a ProxyResequenceReader for the sequence with a custom CompoundSet
     * @param proxyLoader
     * @param compoundSet
     */
    public ChromosomeSequence(SequenceReader<NucleotideCompound> proxyLoader, CompoundSet<NucleotideCompound> compoundSet) {
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

    /**
     * Get the list of genes that have been added to the ChromosomeSequence where accession.toString is the key.
     * The list retains the order the genes are added
     * @return
     */

    public LinkedHashMap<String, GeneSequence> getGeneSequences() {
        return geneSequenceHashMap;
    }

    /**
     *
     * @param accession
     * @return
     */
    public GeneSequence removeGeneSequence(String accession) {
        return geneSequenceHashMap.remove(accession);
    }

    /**
     * Add a gene to the chromosome sequence using bioIndexing starts at 1 instead of 0. The
     * GeneSequence that is returned will have a reference to parent chromosome sequence
     * which actually contains the sequence data. Strand is important for positive and negative
     * direction where negative strand means we need reverse complement. If negative strand then
     * bioBegin will be greater than bioEnd
     *
     *
     * @param accession
     * @param begin
     * @param end
     * @param strand
     * @return
     */
    public GeneSequence addGene(AccessionID accession, int bioBegin, int bioEnd, Strand strand) {
        GeneSequence geneSequence = new GeneSequence(this, bioBegin, bioEnd, strand);
        geneSequence.setAccession(accession);
        geneSequenceHashMap.put(accession.toString(), geneSequence);
        return geneSequence;
    }

    /**
     * Get the gene based on accession. Will return null if not found
     * @param accession
     * @return
     */
    public GeneSequence getGene(String accession) {
        return geneSequenceHashMap.get(accession);
    }
}
