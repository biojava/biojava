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
 * created at 28 Jan 2014
 * Author: ap3 
 */

package org.biojava.nbio.genome.parsers.genename;

import java.io.Serializable;
import java.util.List;

public class GeneChromosomePosition implements Comparable<GeneChromosomePosition>, Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = -6886306238993367835L;
	String geneName;
	String genebankId;
	String chromosome;
	Character orientation;
	Integer transcriptionStart;
	Integer transcriptionEnd;
	Integer cdsStart;
	Integer cdsEnd;
	int exonCount;
	List<Integer> exonStarts;
	List<Integer> exonEnds;
	
	public String getGeneName() {
		return geneName;
	}


	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}


	public String getGenebankId() {
		return genebankId;
	}


	public void setGenebankId(String genebankId) {
		this.genebankId = genebankId;
	}


	public String getChromosome() {
		return chromosome;
	}


	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}


	public Character getOrientation() {
		return orientation;
	}


	public void setOrientation(Character orientation) {
		this.orientation = orientation;
	}


	public Integer getTranscriptionStart() {
		return transcriptionStart;
	}


	public void setTranscriptionStart(Integer transcriptionStart) {
		this.transcriptionStart = transcriptionStart;
	}


	public Integer getTranscriptionEnd() {
		return transcriptionEnd;
	}


	public void setTranscriptionEnd(Integer transcriptionEnd) {
		this.transcriptionEnd = transcriptionEnd;
	}


	public Integer getCdsStart() {
		return cdsStart;
	}


	public void setCdsStart(Integer cdsStart) {
		this.cdsStart = cdsStart;
	}


	public Integer getCdsEnd() {
		return cdsEnd;
	}


	public void setCdsEnd(Integer cdsEnd) {
		this.cdsEnd = cdsEnd;
	}


	public int getExonCount() {
		return exonCount;
	}


	public void setExonCount(int exonCount) {
		this.exonCount = exonCount;
	}


	public List<Integer> getExonStarts() {
		return exonStarts;
	}


	public void setExonStarts(List<Integer> exonStarts) {
		this.exonStarts = exonStarts;
	}


	public List<Integer> getExonEnds() {
		return exonEnds;
	}


	public void setExonEnds(List<Integer> exonEnds) {
		this.exonEnds = exonEnds;
	}


	@Override
	public int compareTo(GeneChromosomePosition o) {
		return geneName.compareTo(o.getGeneName());
	}


	@Override
	public String toString() {
		return "GeneChromosomePosition [geneName=" + geneName + ", genebankId="
				+ genebankId + ", chromosome=" + chromosome + ", orientation="
				+ orientation + ", transcriptionStart=" + transcriptionStart
				+ ", transcriptionEnd=" + transcriptionEnd + ", cdsStart="
				+ cdsStart + ", cdsEnd=" + cdsEnd + ", exonCount=" + exonCount
				+ ", exonStarts=" + exonStarts + ", exonEnds=" + exonEnds + "]";
	}
	
}
