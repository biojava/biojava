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

/**
 *
 * @author Scooter Willis
 */
public class AccessionID {
// GenBank                           gi|gi-number|gb|accession|locus
// EMBL Data Library                 gi|gi-number|emb|accession|locus
// DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
// NBRF PIR                          pir||entry
// Protein Research Foundation       prf||name
// SWISS-PROT UNIPROT                sp|accession|name
// Brookhaven Protein Data Bank (1)  pdb|entry|chain
// Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
// Patents                           pat|country|number
// GenInfo Backbone Id               bbs|number
// General database identifier       gnl|database|identifier
// NCBI Reference Sequence           ref|accession|locus
// Local Sequence identifier         lcl|identifier
    public enum Source {

        GENBANK,EMBL,DDBJ,NBRF,PRF,PDB1,PDB2,PDBe,PATENTS,GENINFO,GENERAL,NCBI, UNIPROT, PFAM, LOCAL, UNKNOWN
    }
    private String id = null;
    private Source source = Source.LOCAL;

    public AccessionID(){
        id = "";
        
    }

    public AccessionID(String id) {
        this.id = id;
        this.source = Source.LOCAL;
    }

    public AccessionID(String id, Source source) {
        this.id = id;
        this.source = source;
    }

    /**
     * @return the id
     */
    public String getID() {
        return id;
    }

    /**
     * @return the source
     */
    public Source getSource() {
        return source;
    }

    @Override
    public String toString() {
        return id;
    }
}
