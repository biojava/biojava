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

/**
 * This class contains the processed data of embl file
 * Primary accession number
 * Sequence version number
 * Topology: 'circular' or 'linear'
 * Molecule type
 * Data class
 * Taxonomic division
 * Sequence length
 *
 * @author Noor Aldeen Al Mbaidin
 * @since 5.0.0
 */
public class EmblId {


    private final String primaryAccession;
    private final String sequenceVersion;
    private final String topology;
    private final String moleculeType;
    private final String dataClass;
    private final String taxonomicDivision;
    private final String sequenceLength;

    public EmblId(String primaryAccession, String sequenceVersion, String topology,
                  String moleculeType, String dataClass, String taxonomicDivision,
                  String sequenceLength) {
        this.primaryAccession = primaryAccession;
        this.sequenceVersion = sequenceVersion;
        this.topology = topology;
        this.moleculeType = moleculeType;
        this.dataClass = dataClass;
        this.taxonomicDivision = taxonomicDivision;
        this.sequenceLength = sequenceLength;
    }

    /**
     * @return String
     */
    public String getPrimaryAccession() {
        return primaryAccession;
    }

    /**
     * return the sequence version
     *
     * @return String
     */
    public String getSequenceVersion() {
        return sequenceVersion;
    }

    public String getTopology() {
        return topology;
    }

    /**
     * Molecule type this represents the type of molecule as stored
     *
     * @return String
     */
    public String getMoleculeType() {
        return moleculeType;
    }

    public String getDataClass() {
        return dataClass;
    }

    /**
     * @return String
     */
    public String getTaxonomicDivision() {
        return taxonomicDivision;
    }

    /**
     * Sequence length The last item on the ID line is the length of the
     * sequence (the total number of bases in the sequence). This number includes
     * base positions reported as present but undetermined (coded as "N").
     *
     * @return String
     */
    public String getSequenceLength() {
        return sequenceLength;
    }

}
