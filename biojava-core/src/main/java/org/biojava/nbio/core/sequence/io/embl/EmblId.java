package org.biojava.nbio.core.sequence.io.embl;

/**
 *
 * This class contains the processed data of embl file
 * Primary accession number
 * Sequence version number
 * Topology: 'circular' or 'linear'
 * Molecule type
 * Data class
 * Taxonomic division
 * Sequence length
 * @since 5.0.0
 * @author Noor Aldeen Al Mbaidin
 */
public class EmblId {


    private String primaryAccession;
    private String sequenceVersion;
    private String topology;
    private String moleculeType;
    private String dataClass;
    private String taxonomicDivision;
    private String sequenceLength;

    public EmblId() {
    }

    /**
     *
     * @return String
     */
    public String getPrimaryAccession() {
        return primaryAccession;
    }

    public void setPrimaryAccession(String primaryAccession) {
        this.primaryAccession = primaryAccession;
    }

    /**
     * return the sequence version
     * @return String
     */
    public String getSequenceVersion() {
        return sequenceVersion;
    }

    public void setSequenceVersion(String sequenceVersion) {
        this.sequenceVersion = sequenceVersion;
    }

    public String getTopology() {
        return topology;
    }

    public void setTopology(String topology) {
        this.topology = topology;
    }

    /**
     * Molecule type this represents the type of molecule as stored
     * @return String
     */
    public String getMoleculeType() {
        return moleculeType;
    }

    public void setMoleculeType(String moleculeType) {
        this.moleculeType = moleculeType;
    }

    public String getDataClass() {
        return dataClass;
    }

    public void setDataClass(String dataClass) {
        this.dataClass = dataClass;
    }

    /**
     *
     * @return String
     */
    public String getTaxonomicDivision() {
        return taxonomicDivision;
    }

    public void setTaxonomicDivision(String taxonomicDivision) {
        this.taxonomicDivision = taxonomicDivision;
    }

    /**
     * Sequence length The last item on the ID line is the length of the
     * sequence (the total number of bases in the sequence). This number includes
     * base positions reported as present but undetermined (coded as "N").
     * @return String
     */
    public String getSequenceLength() {
        return sequenceLength;
    }

    public void setSequenceLength(String sequenceLength) {
        this.sequenceLength = sequenceLength;
    }

}
