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
 * Author: Andreas Prlic
 */

package org.biojava.nbio.genome.parsers.genename;

import java.io.Serializable;

/** 
 * A simple bean that contains gene name information as available from www.genenames.org
 *
 * @author Andreas Prlic
 *
 */
public class GeneName implements Serializable, Comparable<GeneName>{
	//[HGNC ID, Approved Symbol, Approved Name, Status, Previous Symbols, Previous Names, Synonyms,
	// Chromosome, Accession Numbers, RefSeq IDs,Uniprot]

	private static final long serialVersionUID = -7163977639324764020L;
	
	private String hgncId;
	private String approvedSymbol;
	private String approvedName;
	private String status;
	private String previousSymbols;
	private String previousNames;
	private String synonyms;
	private String chromosome;
	private String accessionNr;
	private String refseqIds;
	private String uniprot;
	private String omimId;
	private String ensemblGeneId;

	public String getHgncId() {
		return hgncId;
	}
	public void setHgncId(String hgncId) {
		this.hgncId = hgncId;
	}
	public String getApprovedSymbol() {
		return approvedSymbol;
	}
	public void setApprovedSymbol(String approvedSymbol) {
		this.approvedSymbol = approvedSymbol;
	}
	public String getApprovedName() {
		return approvedName;
	}
	public void setApprovedName(String approvedName) {
		this.approvedName = approvedName;
	}
	public String getStatus() {
		return status;
	}
	public void setStatus(String status) {
		this.status = status;
	}
	public String getPreviousSymbols() {
		return previousSymbols;
	}
	public void setPreviousSymbols(String previousSymbols) {
		this.previousSymbols = previousSymbols;
	}
	public String getPreviousNames() {
		return previousNames;
	}
	public void setPreviousNames(String previousNames) {
		this.previousNames = previousNames;
	}
	public String getSynonyms() {
		return synonyms;
	}
	public void setSynonyms(String synonyms) {
		this.synonyms = synonyms;
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public String getAccessionNr() {
		return accessionNr;
	}
	public void setAccessionNr(String accessionNr) {
		this.accessionNr = accessionNr;
	}
	public String getRefseqIds() {
		return refseqIds;
	}
	public void setRefseqIds(String refseqIds) {
		this.refseqIds = refseqIds;
	}
	public String getUniprot() {
		return uniprot;
	}
	public void setUniprot(String uniprot) {
		this.uniprot = uniprot;
	}

	public String getEnsemblGeneId() {
		return ensemblGeneId;
	}
	public void setEnsemblGeneId(String ensemblGeneId) {
		this.ensemblGeneId = ensemblGeneId;
	}
	public String getOmimId() {
		return omimId;
	}
	public void setOmimId(String omimId) {
		this.omimId = omimId;
	}
	@Override
	public int compareTo(GeneName o) {
		return hgncId.compareTo(o.getHgncId());
	}

	public boolean equals(GeneName o){
		return hgncId.equals(o.getHgncId());
	}
	@Override
	public String toString() {
		return "GeneName [hgncId=" + hgncId + ", approvedSymbol="
				+ approvedSymbol + ", approvedName=" + approvedName
				+ ", status=" + status + ", previousSymbols=" + previousSymbols
				+ ", previousNames=" + previousNames + ", synonyms=" + synonyms
				+ ", chromosome=" + chromosome + ", accessionNr=" + accessionNr
				+ ", refseqIds=" + refseqIds + ", uniprot=" + uniprot
				+ ", omimId=" + omimId + ", ensemblGeneId=" + ensemblGeneId
				+ "]";
	}

}
