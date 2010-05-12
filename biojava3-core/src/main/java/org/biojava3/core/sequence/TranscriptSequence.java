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
import java.util.LinkedHashMap;
import java.util.logging.Logger;

import org.biojava3.core.sequence.transcription.TranscriptionEngine;

/**
 *
 * @author Scooter Willis
 */
public class TranscriptSequence extends DNASequence {


    private Strand strand = Strand.UNDEFINED;
    private static final Logger log = Logger.getLogger(TranscriptSequence.class.getName());
    private final ArrayList<CDSSequence> cdsSequenceList = new ArrayList<CDSSequence>();
    private final LinkedHashMap<String, CDSSequence> cdsSequenceHashMap = new LinkedHashMap<String, CDSSequence>();
    private StartCodonSequence startCodonSequence = null;
    private StopCodonSequence stopCodonSequence = null;

    /**
     *
     * @param parentDNASequence
     * @param begin
     * @param end inclusive of end
     */
    public TranscriptSequence(DNASequence parentDNASequence, int begin, int end, Strand strand) {
        setParentSequence(parentDNASequence);
        setBioBegin(begin);
        setBioEnd(end);
        setStrand(strand);
    }

    /**
     * @return the strand
     */
    public Strand getStrand() {
        return strand;
    }

    /**
     * @param strand the strand to set
     */
    public void setStrand(Strand strand) {
        this.strand = strand;
    }

    /**
     *
     * @param accession
     * @return
     */
    public CDSSequence removeCDS(String accession) {
        for (CDSSequence cdsSequence : cdsSequenceList) {
            if (cdsSequence.getAccession().getID().equals(accession)) {
                cdsSequenceList.remove(cdsSequence);
                cdsSequenceHashMap.remove(accession);
                return cdsSequence;
            }
        }
        return null;
    }

    /**
     *
     * @param accession
     * @param begin
     * @param end
     * @param phase 0,1,2
     * @return
     */
    public CDSSequence addCDS(AccessionID accession, int begin, int end, int phase) throws Exception {
        if (cdsSequenceHashMap.containsKey(accession.getID())) {
            throw new Exception("Duplicate accesion id " + accession.getID());
        }
        CDSSequence cdsSequence = new CDSSequence(this, begin, end, phase); //sense should be the same as parent
        cdsSequence.setAccession(accession);
        cdsSequenceList.add(cdsSequence);
        Collections.sort(cdsSequenceList, new CDSComparator());
        cdsSequenceHashMap.put(accession.getID(), cdsSequence);
        return cdsSequence;
    }

    public DNASequence getDNACodingSequence() {
        StringBuilder sb = new StringBuilder();
        for (CDSSequence cdsSequence : cdsSequenceList) {
            sb.append(cdsSequence.getCodingSequence());
        }
        return new DNASequence(sb.toString().toUpperCase());
    }

    public ProteinSequence getProteinSequence() {
        return getProteinSequence(TranscriptionEngine.getDefault());
    }

    public ProteinSequence getProteinSequence(TranscriptionEngine engine) {
        DNASequence dnaCodingSequence = getDNACodingSequence();
        RNASequence rnaCodingSequence = dnaCodingSequence.getRNASequence(engine);
        ProteinSequence proteinSequence = rnaCodingSequence.getProteinSequence(engine);
        proteinSequence.setAccession(new AccessionID(this.getAccession().getID()));
        return proteinSequence;
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
