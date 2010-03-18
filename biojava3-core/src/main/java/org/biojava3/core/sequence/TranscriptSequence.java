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

import java.util.ArrayList;
import java.util.Collections;
import java.util.logging.Logger;


import org.biojava3.core.sequence.transcription.TranscriptionEngine;

/**
 *
 * @author Scooter Willis
 */
public class TranscriptSequence extends DNASequence {

    private static final Logger log = Logger.getLogger(TranscriptSequence.class.getName());

    public enum Sense {

        POSITIVE, NEGATIVE, UNDEFINED
    }
    private StartCodonSequence startCodonSequence = null;
    private StopCodonSequence stopCodonSequence = null;
    boolean intronAdded = false; // need to deal with the problem that typically introns are not added when validating the list and adding in introns as the regions not included in exons
    private final ArrayList<IntronSequence> intronSequenceList = new ArrayList<IntronSequence>();
    private final ArrayList<ExonSequence> exonSequenceList = new ArrayList<ExonSequence>();
    private Sense sense = Sense.UNDEFINED;

    /**
     *
     * @param parentDNASequence
     * @param begin
     * @param end inclusive of end
     */
    public TranscriptSequence(DNASequence parentDNASequence, int begin, int end, Sense sense) {
        setParentDNASequence(parentDNASequence);
        setBegin(begin);
        setEnd(end);
        this.sense = sense;
    }

    public void validate() {
        ExonComparator exonComparator = new ExonComparator();
        //sort based on start position and sense;
        Collections.sort(exonSequenceList, exonComparator);
        if (intronAdded) {
            log.severe( this.getAccession() + " has introns added which will not be handled properly trying to fill in introns gaps from validate method");
        }

        
    //    log.severe("Add in support for building introns based on added exons");

    }

    public Sense getSense() {
        return sense;
    }

    /**
     *
     * @param accession
     * @return
     */
    public IntronSequence removeIntron(String accession) {
        for (IntronSequence intronSequence : intronSequenceList) {
            if (intronSequence.getAccession().getID().equals(accession)) {
                intronSequenceList.remove(intronSequence);
                return intronSequence;
            }
        }
        return null;
    }

    /**
     *
     * @param accession
     * @param begin
     * @param end
     * @return
     */
    public IntronSequence addIntron(AccessionID accession, int begin, int end) {
        intronAdded = true;
        IntronSequence intronSequence = new IntronSequence(this, begin, end, 0, sense); // working off the assumption that intron frame is always 0 or doesn't matter and same sense as parent
        intronSequence.setAccession(accession);
        intronSequenceList.add(intronSequence);
        return intronSequence;
    }

    /**
     *
     * @param accession
     * @return
     */
    public ExonSequence removeExon(String accession) {
        for (ExonSequence exonSequence : exonSequenceList) {
            if (exonSequence.getAccession().getID().equals(accession)) {
                exonSequenceList.remove(exonSequence);
                validate();
                return exonSequence;
            }
        }
        return null;
    }

    /**
     *
     * @param accession
     * @param begin
     * @param end
     * @return
     */
    public ExonSequence addExon(AccessionID accession, int begin, int end) {
        ExonSequence exonSequence = new ExonSequence(this, begin, end, sense); //sense should be the same as parent
        exonSequence.setAccession(accession);
        exonSequenceList.add(exonSequence);
        validate();
        return exonSequence;
    }

    public RNASequence getRNACodingSequence() {
        StringBuilder sb = new StringBuilder();
        System.err.println("Need to sort exon list");
        for (ExonSequence exonSequence : exonSequenceList) {
            sb.append(exonSequence.toString());
        }
        return new RNASequence(sb.toString());
    }

    public ProteinSequence getProteinSequence() {
        return getProteinSequence(TranscriptionEngine.getDefault());
    }

    public ProteinSequence getProteinSequence(TranscriptionEngine engine) {
        RNASequence rnaCodingSequence = getRNASequence();
        return rnaCodingSequence.getProteinSequence(engine);
    }

    /**
     * @return the startCodonSequence
     */
    public StartCodonSequence getStartCodonSequence() {
        return startCodonSequence;
    }

    /**
     * @param startCodonSequence the startCodonSequence to set
     */
    public void addStartCodonSequence(AccessionID accession, int begin, int end) {
        this.startCodonSequence = new StartCodonSequence(this, begin, end);
        startCodonSequence.setAccession(accession);
    }

    /**
     * @return the stopCodonSequence
     */
    public StopCodonSequence getStopCodonSequence() {
        return stopCodonSequence;
    }

    /**
     * @param stopCodonSequence the stopCodonSequence to set
     */
    public void addStopCodonSequence(AccessionID accession, int begin, int end) {
        this.stopCodonSequence = new StopCodonSequence(this, begin, end);
        stopCodonSequence.setAccession(accession);
    }
}
