/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojavax.bio.seq.io;

/**
 * This class is used by <code>FastaFormat</code> to determine which fields are in the 
 * fasta header. By default they all are except for the sequence name. This is for
 * compliance with fasta files that come from Genbank where the name is derived
 * from the accession number so need not be repeated.
 * The class can be used to customise
 * what appears. Eg if you only want the accession set everything else false.
 * Note that if fields in the <code>RichSequence</code> being parsed by the 
 * <code>FastaFormat</code> object then they may not be in the header even if
 * they are specified in this class. 
 * @author Mark Schreiber
 * @since 1.6
 */
public class FastaHeader {
    private boolean showIdentifier = true;
    private boolean showNamespace = true;
    private boolean showAccession = true;
    private boolean showVersion = true;
    private boolean showName = false;
    private boolean showDescription = true;

    public boolean isShowIdentifier() {
        return showIdentifier;
    }

    public void setShowIdentifier(boolean showIdentifier) {
        this.showIdentifier = showIdentifier;
    }

    public boolean isShowNamespace() {
        return showNamespace;
    }

    public void setShowNamespace(boolean showNamespace) {
        this.showNamespace = showNamespace;
    }

    public boolean isShowAccession() {
        return showAccession;
    }

    public void setShowAccession(boolean showAccession) {
        this.showAccession = showAccession;
    }

    public boolean isShowVersion() {
        return showVersion;
    }

    /**
     * Determines if the version number of a sequence should be displayed. If
     * there is no accession this may not make much sense.
     * @param showVersion
     */
    public void setShowVersion(boolean showVersion) {
        this.showVersion = showVersion;
    }

    public boolean isShowName() {
        return showName;
    }

    public void setShowName(boolean showName) {
        this.showName = showName;
    }

    public boolean isShowDescription() {
        return showDescription;
    }

    public void setShowDescription(boolean showDescription) {
        this.showDescription = showDescription;
    }
}
