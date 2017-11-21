package org.biojava.nbio.core.sequence.io.embl;

/**
 * This class should parse the data of embl file identification
 * @author Noor Aldeen Al Mbaidin
 */
public class EmblId {


    private String PrimaryAccession;
    private String SequenceVersion;
    private String Topology;
    private String MoleculeType;
    private String DataClass;
    private String TaxonomicDivision;
    private String SequenceLength;

    public EmblId() {
    }

    public String getPrimaryAccession() {
        return PrimaryAccession;
    }

    public void setPrimaryAccession(String primaryAccession) {
        PrimaryAccession = primaryAccession;
    }

    public String getSequenceVersion() {
        return SequenceVersion;
    }

    public void setSequenceVersion(String sequenceVersion) {
        SequenceVersion = sequenceVersion;
    }

    public String getTopology() {
        return Topology;
    }

    public void setTopology(String topology) {
        Topology = topology;
    }

    public String getMoleculeType() {
        return MoleculeType;
    }

    public void setMoleculeType(String moleculeType) {
        MoleculeType = moleculeType;
    }

    public String getDataClass() {
        return DataClass;
    }

    public void setDataClass(String dataClass) {
        DataClass = dataClass;
    }

    public String getTaxonomicDivision() {
        return TaxonomicDivision;
    }

    public void setTaxonomicDivision(String taxonomicDivision) {
        TaxonomicDivision = taxonomicDivision;
    }

    public String getSequenceLength() {
        return SequenceLength;
    }

    public void setSequenceLength(String sequenceLength) {
        SequenceLength = sequenceLength;
    }

}
