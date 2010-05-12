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



/**
 *
 * @author Scooter Willis
 */
public class GeneSequence extends DNASequence {

    private final LinkedHashMap<String, TranscriptSequence> transcriptSequenceHashMap = new LinkedHashMap<String, TranscriptSequence>();
    private static final Logger log = Logger.getLogger(GeneSequence.class.getName());
    private final LinkedHashMap<String, IntronSequence> intronSequenceHashMap = new LinkedHashMap<String, IntronSequence>();
    private final LinkedHashMap<String, ExonSequence> exonSequenceHashMap = new LinkedHashMap<String, ExonSequence>();
    private final ArrayList<IntronSequence> intronSequenceList = new ArrayList<IntronSequence>();
    private final ArrayList<ExonSequence> exonSequenceList = new ArrayList<ExonSequence>();
    boolean intronAdded = false; // need to deal with the problem that typically introns are not added when validating the list and adding in introns as the regions not included in exons

    /**
     *
     * @param parentDNASequence
     * @param begin
     * @param end inclusive of end
     */
    public GeneSequence(DNASequence parentDNASequence, int begin, int end) {
        setParentSequence(parentDNASequence);
        setBioBegin(begin);
        setBioEnd(end);
    }

    public void validate() {
        ExonComparator exonComparator = new ExonComparator();
        //sort based on start position and sense;
        Collections.sort(exonSequenceList, exonComparator);
        if (intronAdded) {
            log.severe(this.getAccession() + " has introns added which will not be handled properly trying to fill in introns gaps from validate method");
        }


        //    log.severe("Add in support for building introns based on added exons");

    }

    public TranscriptSequence getTranscript(String accession) {
        return transcriptSequenceHashMap.get(accession);
    }

    public LinkedHashMap<String, TranscriptSequence> getTranscripts() {
        return transcriptSequenceHashMap;
    }

    public TranscriptSequence removeTranscript(String accession) {


        return transcriptSequenceHashMap.remove(accession);
    }

    public TranscriptSequence addTranscript(AccessionID accession, int begin, int end, Strand strand) throws Exception {
        if (transcriptSequenceHashMap.containsKey(accession.getID())) {
            throw new Exception("Duplicate accesion id " + accession.getID());
        }
        TranscriptSequence transcriptSequence = new TranscriptSequence(this, begin, end, strand);
        transcriptSequence.setAccession(accession);
        transcriptSequenceHashMap.put(accession.getID(), transcriptSequence);
        return transcriptSequence;
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
                intronSequenceHashMap.remove(accession);
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
    public IntronSequence addIntron(AccessionID accession, int begin, int end) throws Exception {
        if (intronSequenceHashMap.containsKey(accession.getID())) {
            throw new Exception("Duplicate accesion id " + accession.getID());
        }
        intronAdded = true;
        IntronSequence intronSequence = new IntronSequence(this, begin, end); // working off the assumption that intron frame is always 0 or doesn't matter and same sense as parent
        intronSequence.setAccession(accession);
        intronSequenceList.add(intronSequence);
        intronSequenceHashMap.put(accession.getID(), intronSequence);
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
                exonSequenceHashMap.remove(accession);
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
    public ExonSequence addExon(AccessionID accession, int begin, int end) throws Exception {
        if (exonSequenceHashMap.containsKey(accession.getID())) {
            throw new Exception("Duplicate accesion id " + accession.getID());
        }

        ExonSequence exonSequence = new ExonSequence(this, begin, end); //sense should be the same as parent
        exonSequence.setAccession(accession);
        exonSequenceList.add(exonSequence);
        exonSequenceHashMap.put(accession.getID(), exonSequence);
        validate();
        return exonSequence;
    }
}
