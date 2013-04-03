/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

/**
 *<pre>
 * GenBank                           gi|gi-number|gb|accession|locus
 * ENA Data Library                  gi|gi-number|emb|accession|locus
 * DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
 * NBRF PIR                          pir||entry
 * Protein Research Foundation       prf||name
 * SWISS-PROT UNIPROT                sp|accession|name
 * Brookhaven Protein Data Bank (1)  pdb|entry|chain
 * Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
 * Patents                           pat|country|number
 * GenInfo Backbone Id               bbs|number
 * General database identifier       gnl|database|identifier
 * NCBI Reference Sequence           ref|accession|locus
 * Local Sequence identifier         lcl|identifier
 * </pr>
 * @author Scooter Willis <willishf at gmail dot com>
 */

public enum DataSource {

    GENBANK, ENA, DDBJ, NBRF, PRF, PDB1, PDB2, PDBe, PATENTS, GENINFO, GENERAL, NCBI, UNIPROT, PFAM, LOCAL, UNKNOWN
}
