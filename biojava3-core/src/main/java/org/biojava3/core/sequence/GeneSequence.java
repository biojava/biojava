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
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;

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
    private Strand strand = Strand.UNDEFINED;
    private ChromosomeSequence chromosomeSequence;

    /**
     * A class that keeps track of the details of a GeneSequence which is difficult to properly model. Two important concepts that is difficult
     * to make everything flexible but still work. You can have GFF features that only describe Exons or Exons/Introns or CDS regions and one
     * or more Transcriptions. You can have exon sequences but that does not imply transcription to the actual protein.
     *
     * The GeneSequence will keep track of Exons and Introns but to get a Protein sequence you need to start with a
     * TranscriptSequence and then add CDS sequences.
     *
     * This is also a key class in the biojava-3-genome module for reading and writing GFF3 files
     *
     * @param parentDNASequence
     * @param begin
     * @param end inclusive of end
     * @param strand force a gene to have strand and transcription sequence will inherit
     */
    public GeneSequence(ChromosomeSequence parentSequence, int begin, int end, Strand strand) {
        chromosomeSequence = parentSequence;
        setParentSequence(parentSequence);
        setBioBegin(begin);
        setBioEnd(end);
        setStrand(strand);
    }

    /**
     * The parent ChromosomeSequence which contains the actual DNA sequence data
     * @return Chromosome sequence
     */
    public ChromosomeSequence getParentChromosomeSequence() {
        return chromosomeSequence;
    }

    @Override
    public int getLength() {
        return Math.abs(this.getBioEnd() - this.getBioBegin()) + 1;
    }

    /**
     * Once everything has been added to the gene sequence where you might have added exon sequences only then you
     * can infer the intron sequences and add them. You may also have the case where you only added one or more
     * TranscriptSequences and from that you can infer the exon sequences and intron sequences.
     * Currently not implement
     */
    public void addIntronsUsingExons() throws Exception {
        if (intronAdded) { //going to assume introns are correct
            return;
        }
        if (exonSequenceList.size() == 0) {
            return;
        }
        ExonComparator exonComparator = new ExonComparator();
        //sort based on start position and sense;
        Collections.sort(exonSequenceList, exonComparator);
        int shift = -1;
        if (getStrand() == Strand.NEGATIVE) {
            shift = 1;
        }
        ExonSequence firstExonSequence = exonSequenceList.get(0);
        int intronIndex = 1;
 //       if (firstExonSequence.getBioBegin().intValue() != getBioBegin().intValue()) {
 //           this.addIntron(new AccessionID(this.getAccession().getID() + "-" + "intron" + intronIndex), getBioBegin(), firstExonSequence.getBioBegin() + shift);
 //           intronIndex++;
 //       }
        for (int i = 0; i < exonSequenceList.size() - 1; i++) {
            ExonSequence exon1 = exonSequenceList.get(i);
            ExonSequence exon2 = exonSequenceList.get(i + 1);
            this.addIntron(new AccessionID(this.getAccession().getID() + "-" + "intron" + intronIndex), exon1.getBioEnd() - shift, exon2.getBioBegin() + shift);
            intronIndex++;
        }

 //       ExonSequence lastExonSequence = exonSequenceList.get(exonSequenceList.size() - 1);
 //       if (lastExonSequence.getBioEnd().intValue() != getBioEnd().intValue()) {
 //           this.addIntron(new AccessionID(this.getAccession().getID() + "-" + "intron" + intronIndex), lastExonSequence.getBioEnd() - shift, getBioEnd());
 //           intronIndex++;
 //       }

        //    log.severe("Add in support for building introns based on added exons");

    }

    /**
     * A gene should have Strand
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
     * Get the transcript sequence by accession
     * @param accession
     * @return the transcript
     */
    public TranscriptSequence getTranscript(String accession) {
        return transcriptSequenceHashMap.get(accession);
    }

    /**
     * Get the collection of transcription sequences assigned to this gene
     * @return transcripts
     */
    public LinkedHashMap<String, TranscriptSequence> getTranscripts() {
        return transcriptSequenceHashMap;
    }

    /**
     * Remove the transcript sequence from the gene
     * @param accession
     * @return transcriptsequence
     */
    public TranscriptSequence removeTranscript(String accession) {


        return transcriptSequenceHashMap.remove(accession);
    }

    /**
     * Add a transcription sequence to a gene which describes a ProteinSequence
     * @param accession
     * @param begin
     * @param end
     * @return transcript sequence
     * @throws Exception If the accession id is already used
     */
    public TranscriptSequence addTranscript(AccessionID accession, int begin, int end) throws Exception {
        if (transcriptSequenceHashMap.containsKey(accession.getID())) {
            throw new Exception("Duplicate accesion id " + accession.getID());
        }
        TranscriptSequence transcriptSequence = new TranscriptSequence(this, begin, end);
        transcriptSequence.setAccession(accession);
        transcriptSequenceHashMap.put(accession.getID(), transcriptSequence);
        return transcriptSequence;
    }

    /**
     * Remove the intron by accession
     * @param accession
     * @return intron sequence
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
     * Add an Intron Currently used to mark an IntronSequence as a feature
     * @param accession
     * @param begin
     * @param end
     * @return intron sequence
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
     * Remove the exon sequence
     * @param accession
     * @return exon sequence
     */
    public ExonSequence removeExon(String accession) {
        for (ExonSequence exonSequence : exonSequenceList) {
            if (exonSequence.getAccession().getID().equals(accession)) {
                exonSequenceList.remove(exonSequence);
                exonSequenceHashMap.remove(accession);
                // we now have a new gap which creates an intron
                intronSequenceList.clear();
                intronSequenceHashMap.clear();
                intronAdded = false;
                try{
                    addIntronsUsingExons();
                }catch(Exception e){
                    log.severe("Remove Exon validate() error " + e.getMessage());
                }
                return exonSequence;
            }
        }
        return null;
    }

    /**
     * Add an ExonSequence mainly used to mark as a feature
     * @param accession
     * @param begin
     * @param end
     * @return exon sequence
     */
    public ExonSequence addExon(AccessionID accession, int begin, int end) throws Exception {
        if (exonSequenceHashMap.containsKey(accession.getID())) {
            throw new Exception("Duplicate accesion id " + accession.getID());
        }

        ExonSequence exonSequence = new ExonSequence(this, begin, end); //sense should be the same as parent
        exonSequence.setAccession(accession);
        exonSequenceList.add(exonSequence);
        exonSequenceHashMap.put(accession.getID(), exonSequence);
        return exonSequence;
    }

    /**
     * Get the exons as an ArrayList
     * @return exons
     */
    public ArrayList<ExonSequence> getExonSequences() {
        return exonSequenceList;
    }

    /**
     * Get the introns as an ArrayList
     * @return introns 
     */
    public ArrayList<IntronSequence> getIntronSequences() {
        return intronSequenceList;
    }

    /**
     * Try to give method clarity where you want a DNASequence coding in the 5' to 3' direction
     * Returns the DNASequence representative of the 5' and 3' reading based on strand
     * @return dna sequence
     */
    public DNASequence getSequence5PrimeTo3Prime() {
        String sequence = getSequenceAsString(this.getBioBegin(), this.getBioEnd(), this.getStrand());
        if (getStrand() == Strand.NEGATIVE) {
            //need to take complement of sequence because it is negative and we are returning the gene sequence from the opposite strand
            StringBuilder b = new StringBuilder(getLength());
            CompoundSet<NucleotideCompound> compoundSet = this.getCompoundSet();
            for (int i = 0; i < sequence.length(); i++) {
                String nucleotide = sequence.charAt(i) + "";
                NucleotideCompound nucleotideCompound = compoundSet.getCompoundForString(nucleotide);
                b.append(nucleotideCompound.getComplement().getShortName());
            }
            sequence = b.toString();
        }
        DNASequence dnaSequence = new DNASequence(sequence.toUpperCase());
        dnaSequence.setAccession(new AccessionID(this.getAccession().getID()));
        return dnaSequence;
    }
}
