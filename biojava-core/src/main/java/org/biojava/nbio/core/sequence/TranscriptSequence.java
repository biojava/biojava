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
package org.biojava.nbio.core.sequence;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.transcription.TranscriptionEngine;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;

/**
 * This is the sequence if you want to go from a gene sequence to a protein sequence. Need to start with a
 * ChromosomeSequence then getting a GeneSequence and then a TranscriptSequence
 * @author Scooter Willis
 */
public class TranscriptSequence extends GenomeSequence {

	private final static Logger logger = LoggerFactory.getLogger(TranscriptSequence.class);

	private final ArrayList<CDSSequence> cdsSequenceList = new ArrayList<CDSSequence>();
	private final LinkedHashMap<String, CDSSequence> cdsSequenceHashMap = new LinkedHashMap<String, CDSSequence>();
	private StartCodonSequence startCodonSequence = null;
	private StopCodonSequence stopCodonSequence = null;
	private GeneSequence parentGeneSequence = null;

	/**
	 * Use {@code}public TranscriptSequence(GeneSequence parentDNASequence, AccessionID accessionID, int begin, int end){@code}
	 * that requires an explicit accessionID
	 * @deprecated
	 */
	public TranscriptSequence(GeneSequence parentDNASequence, int begin, int end) {
		setCompoundSet(DNACompoundSet.getDNACompoundSet());
		try {
			initSequenceStorage(parentDNASequence.getSequenceAsString());
		} catch (CompoundNotFoundException e) {
			throw new IllegalArgumentException(e);
		}
		setParentSequence(parentDNASequence);
		this.parentGeneSequence = parentDNASequence;
		setBioBegin(begin);
		setBioEnd(end);
	}

	/**
	 *
	 * @param parentDNASequence
	 * @param accessionID
	 * @param begin
	 * @param end inclusive of end
	 * @throws  IllegalArgumentException if the parentDNASequence is incompatible with DNACompoundSet
	 */
	public TranscriptSequence(GeneSequence parentDNASequence, AccessionID accessionID, int begin, int end) {
		this(parentDNASequence, begin, end);
		setAccession(accessionID);
	}

	/**
	 * @return the strand
	 */
	public Strand getStrand() {
		return parentGeneSequence.getStrand();
	}

	/**
	 * Remove a CDS or coding sequence from the transcript sequence
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
	 * Get the CDS sequences that have been added to the TranscriptSequences
	 * @return
	 */
	public LinkedHashMap<String, CDSSequence> getCDSSequences() {
		return cdsSequenceHashMap;
	}

	/**
	 * Add a Coding Sequence region with phase to the transcript sequence
	 * @param accession
	 * @param begin
	 * @param end
	 * @param phase 0,1,2
	 * @return
	 */
	public CDSSequence addCDS(AccessionID accession, int begin, int end, int phase) throws Exception {
		if (cdsSequenceHashMap.containsKey(accession.getID())) {
			throw new Exception("Duplicate accession id " + accession.getID());
		}
		CDSSequence cdsSequence = new CDSSequence(this, begin, end, phase); //sense should be the same as parent
		cdsSequence.setAccession(accession);
		cdsSequenceList.add(cdsSequence);
		Collections.sort(cdsSequenceList, new CDSComparator());
		cdsSequenceHashMap.put(accession.getID(), cdsSequence);
		return cdsSequence;
	}

	/**
	 * http://www.sequenceontology.org/gff3.shtml
	 * http://biowiki.org/~yam/bioe131/GFF.ppt
	 * @return
	 */
	/**
	 * Return a list of protein sequences based on each CDS sequence
	 * where the phase shift between two CDS sequences is assigned to the
	 * CDS sequence that starts the triplet. This can be used to map
	 * a CDS/exon region of a protein sequence back to the DNA sequence
	 * If you have a protein sequence and a predicted gene you can take the
	 * predict CDS protein sequences and align back to the protein sequence.
	 * If you have errors in mapping the predicted protein CDS regions to
	 * an the known protein sequence then you can identify possible errors
	 * in the prediction
	 *
	 * @return
	 */
	public ArrayList<ProteinSequence> getProteinCDSSequences() {
		ArrayList<ProteinSequence> proteinSequenceList = new ArrayList<ProteinSequence>();
		for (int i = 0; i < cdsSequenceList.size(); i++) {
			CDSSequence cdsSequence = cdsSequenceList.get(i);
			String codingSequence = cdsSequence.getCodingSequence();
			//          logger.debug("CDS {} {} = {}", getStrand(), cdsSequence.getPhase(), codingSequence);
			if (this.getStrand() == Strand.NEGATIVE) {
				if (cdsSequence.phase == 1) {
					codingSequence = codingSequence.substring(1, codingSequence.length());
				} else if (cdsSequence.phase == 2) {
					codingSequence = codingSequence.substring(2, codingSequence.length());
				}
				if (i < cdsSequenceList.size() - 1) {
					CDSSequence nextCDSSequence = cdsSequenceList.get(i + 1);
					if (nextCDSSequence.phase == 1) {
						String nextCodingSequence = nextCDSSequence.getCodingSequence();
						codingSequence = codingSequence + nextCodingSequence.substring(0, 1);
					} else if (nextCDSSequence.phase == 2) {
						String nextCodingSequence = nextCDSSequence.getCodingSequence();
						codingSequence = codingSequence + nextCodingSequence.substring(0, 2);
					}
				}
			} else {
				if (cdsSequence.phase == 1) {
					codingSequence = codingSequence.substring(1, codingSequence.length());
				} else if (cdsSequence.phase == 2) {
					codingSequence = codingSequence.substring(2, codingSequence.length());
				}
				if (i < cdsSequenceList.size() - 1) {
					CDSSequence nextCDSSequence = cdsSequenceList.get(i + 1);
					if (nextCDSSequence.phase == 1) {
						String nextCodingSequence = nextCDSSequence.getCodingSequence();
						codingSequence = codingSequence + nextCodingSequence.substring(0, 1);
					} else if (nextCDSSequence.phase == 2) {
						String nextCodingSequence = nextCDSSequence.getCodingSequence();
						codingSequence = codingSequence + nextCodingSequence.substring(0, 2);
					}
				}
			}


			//    logger.debug("Coding Sequence: {}", codingSequence);

			DNASequence dnaCodingSequence = null;
			try {
				dnaCodingSequence = new DNASequence(codingSequence.toUpperCase());
			} catch (CompoundNotFoundException e) {
				// if I understand this should not happen, please correct if I'm wrong - JD 2014-10-24
				logger.error("Could not create DNA coding sequence, {}. This is most likely a bug.", e.getMessage());
			}
			RNASequence rnaCodingSequence = dnaCodingSequence.getRNASequence(TranscriptionEngine.getDefault());
			ProteinSequence proteinSequence = rnaCodingSequence.getProteinSequence(TranscriptionEngine.getDefault());
			proteinSequence.setAccession(new AccessionID(cdsSequence.getAccession().getID()));
			proteinSequence.setParentDNASequence(cdsSequence, 1, cdsSequence.getLength());
			proteinSequenceList.add(proteinSequence);
		}
		return proteinSequenceList;
	}

	/**
	 * Get the stitched together CDS sequences then maps to the cDNA
	 * @return
	 */
	public DNASequence getDNACodingSequence() {
		StringBuilder sb = new StringBuilder();
		for (CDSSequence cdsSequence : cdsSequenceList) {
			sb.append(cdsSequence.getCodingSequence());
		}

		DNASequence dnaSequence = null;
		try {
			dnaSequence = new DNASequence(sb.toString().toUpperCase());
		} catch (CompoundNotFoundException e) {
			// if I understand this should not happen, please correct if I'm wrong - JD 2014-10-24
			logger.error("Could not create DNA coding sequence, {}. This is most likely a bug.", e.getMessage());
		}
		dnaSequence.setAccession(new AccessionID(this.getAccession().getID()));
		return dnaSequence;
	}

	/**
	 * Get the protein sequence
	 * @return
	 */
	public ProteinSequence getProteinSequence() {
		return getProteinSequence(TranscriptionEngine.getDefault());
	}

	/**
	 * Get the protein sequence with user defined TranscriptEngine
	 * @param engine
	 * @return
	 */
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
	 * Sets the start codon sequence at given begin /  end location. Note that calling this method multiple times
	 * will replace any existing value.
	 * @param accession
	 * @param begin
	 * @param end
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
	 * Sets the stop codon sequence at given begin /  end location. Note that calling this method multiple times
	 * will replace any existing value.
	 * @param accession
	 * @param begin
	 * @param end
	 */
	public void addStopCodonSequence(AccessionID accession, int begin, int end) {
		this.stopCodonSequence = new StopCodonSequence(this, begin, end);
		stopCodonSequence.setAccession(accession);
	}
}
