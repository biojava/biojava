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
 */
package org.biojava.nbio.core.sequence.io.embl;

import java.util.LinkedList;
import java.util.List;


/**
 * this class contains the parsed data of embl file
 *
 * @author Noor Aldeen Al Mbaidin
 * @since 5.0.0
 */

public class EmblRecord {

    private EmblId emblId;
    private List<EmblReference> emblReference;
    private List<String> accessionNumber = new LinkedList<>();
    private String projectIdentifier;
    private String orGanelle;
    private String createdDate;
    private String featureHeader;
    private String featureTable;
    private String lastUpdatedDate;
    private String sequenceDescription;
    private List<String> keyword = new LinkedList<>();
    private String organismSpecies;
    private String organismClassification;
    private String databaseCrossReference;
    private String assemblyHeader;
    private String assemblyInformation;
    private String constructedSequence;
    private String sequenceHeader;
    private String sequence;

    /**
     * The ID (IDentification) line
     * The tokens represent:
     * 1. Primary accession number
     * 2. Sequence version number
     * 3. Topology: 'circular' or 'linear'
     * 4. Molecule type
     * 5. Data class
     * 6. Taxonomic division
     * 7. Sequence length
     *
     * @return EmblId
     */
    public EmblId getEmblId() {
        return emblId;
    }

    public void setEmblId(EmblId emblId) {
        this.emblId = emblId;
    }

    /**
     * The Reference (RN, RC, RP, RX, RG, RA, RT, RL) Lines
     * These lines comprise the literature citations within the database.
     * The citations provide access to the papers from which the data has been
     * abstracted.
     *
     * @return EmblReference
     */
    public List<EmblReference> getEmblReference() {
        return emblReference;
    }

    public void setEmblReference(List<EmblReference> emblReference) {
        this.emblReference = emblReference;
    }

    /**
     * The AC (Accession number) line lists the accession numbers associated with
     * the entry.
     *
     * @return List<String>
     */
    public List<String> getAccessionNumber() {
        return accessionNumber;
    }

    public void setAccessionNumber(List<String> accessionNumber) {
        this.accessionNumber = accessionNumber;
    }

    /**
     * @return String
     */
    public String getProjectIdentifier() {
        return projectIdentifier;
    }

    public void setProjectIdentifier(String projectIdentifier) {
        this.projectIdentifier = projectIdentifier;
    }

    /**
     * The OG (OrGanelle) linetype indicates the sub-cellular location of non-nuclear
     * sequences.
     *
     * @return String
     */
    public String getOrGanelle() {
        return orGanelle;
    }

    public void setOrGanelle(String orGanelle) {
        this.orGanelle = orGanelle;
    }

    /**
     * The DT  line shows when an entry first appeared in the database
     *
     * @return String
     */
    public String getCreatedDate() {
        return createdDate;
    }

    public void setCreatedDate(String createdDate) {
        this.createdDate = createdDate;
    }

    /**
     * The FH (Feature Header) lines are present only to improve readability of
     * an entry when it is printed or displayed on a terminal screen.
     *
     * @return String
     */
    public String getFeatureHeader() {
        return featureHeader;
    }

    public void setFeatureHeader(String featureHeader) {
        this.featureHeader = featureHeader;
    }

    /**
     * The FT (Feature Table) lines provide a mechanism for the annotation of the
     * sequence data. Regions or sites in the sequence which are of interest are
     * listed in the table.
     *
     * @return String
     */
    public String getFeatureTable() {
        return featureTable;
    }

    public void setFeatureTable(String featureTable) {
        this.featureTable = featureTable;
    }

    /**
     * The DT (DaTe) line shows when an entry was last updated in the database.
     *
     * @return String
     */
    public String getLastUpdatedDate() {
        return lastUpdatedDate;
    }

    public void setLastUpdatedDate(String lastUpdatedDate) {
        this.lastUpdatedDate = lastUpdatedDate;
    }

    /**
     * The DE (Description) lines contain general descriptive information about the
     * sequence stored. This may include the designations of genes for which the
     * sequence codes, the region of the genome from which it is derived, or other
     * information which helps to identify the sequence.
     *
     * @return String
     */
    public String getSequenceDescription() {
        return sequenceDescription;
    }

    public void setSequenceDescription(String sequenceDescription) {
        this.sequenceDescription = sequenceDescription;
    }

    /**
     * The KW (KeyWord) lines provide information which can be used to generate
     * cross-reference indexes of the sequence entries based on functional,
     * structural, or other categories deemed important.
     *
     * @return List<String>
     */
    public List<String> getKeyword() {
        return keyword;
    }

    public void setKeyword(List<String> keyword) {
        this.keyword = keyword;
    }

    /**
     * The OS (Organism Species) line specifies the preferred scientific name of
     * the organism which was the source of the stored sequence. In most
     * cases this is done by giving the Latin genus and species designations,
     * followed (in parentheses) by the preferred common name in English where known.
     *
     * @return String
     */
    public String getOrganismSpecies() {
        return organismSpecies;
    }

    public void setOrganismSpecies(String organismSpecies) {
        this.organismSpecies = organismSpecies;
    }

    /**
     * The OC (Organism Classification) lines contain the taxonomic classification
     * Of the source organism
     *
     * @return String
     */
    public String getOrganismClassification() {
        return organismClassification;
    }

    public void setOrganismClassification(String organismClassification) {
        this.organismClassification = organismClassification;
    }

    /**
     * The DR (Database Cross-reference) line cross-references other databases which
     * contain information related to the entry in which the DR line appears.
     *
     * @return String
     */
    public String getDatabaseCrossReference() {
        return databaseCrossReference;
    }

    public void setDatabaseCrossReference(String databaseCrossReference) {
        this.databaseCrossReference = databaseCrossReference;
    }

    /**
     * The AH (Assembly Header) line provides column headings for the assembly information.
     *
     * @return String
     */
    public String getAssemblyHeader() {
        return assemblyHeader;
    }

    public void setAssemblyHeader(String assemblyHeader) {
        this.assemblyHeader = assemblyHeader;
    }

    /**
     * The AS (Assembly Information) lines provide information on the composition of
     * a TPA or TSA sequence.
     *
     * @return String
     */
    public String getAssemblyInformation() {
        return assemblyInformation;
    }

    public void setAssemblyInformation(String assemblyInformation) {
        this.assemblyInformation = assemblyInformation;
    }

    /**
     * Con(structed) sequences in the CON data classes represent complete
     * chromosomes, genomes and other long sequences constructed from segment entries.
     *
     * @return String
     */
    public String getConstructedSequence() {
        return constructedSequence;
    }

    public void setConstructedSequence(String constructedSequence) {
        this.constructedSequence = constructedSequence;
    }

    /**
     * The SQ (SeQuence header) line marks the beginning of the sequence data and
     * Gives a summary of its content.
     *
     * @return String
     */
    public String getSequenceHeader() {
        return sequenceHeader;
    }

    public void setSequenceHeader(String sequenceHeader) {
        this.sequenceHeader = sequenceHeader;
    }

    /**
     * The Sequence Data Line
     *
     * @return String
     */
    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

}
