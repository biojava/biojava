/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import org.biojava3.core.exceptions.HeaderParseException;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.io.template.HeaderParserInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.AbstractSequence.AnnotationType;

/**
 * GenBank                           gi|gi-number|gb|accession|locus
 * EMBL Data Library                 gi|gi-number|emb|accession|locus
 * DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
 * NBRF PIR                          pir||entry
 * Protein Research Foundation       prf||name
 * SWISS-PROT                        sp|accession|name
 * Brookhaven Protein Data Bank (1)  pdb|entry|chain
 * Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
 * RDB EBI                           PDB:1ECY_A mol:protein length:142  ECOTIN
 * Patents                           pat|country|number
 * GenInfo Backbone Id               bbs|number
 * General database identifier       gnl|database|identifier
 * NCBI Reference Sequence           ref|accession|locus
 * Local Sequence identifier         lcl|identifier
 * 
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GenericFastaHeaderParser implements HeaderParserInterface {

    public void parseHeader(String header, AbstractSequence sequence) {
        //uniptrot 
        // tr|Q0TET7|Q0TET7_ECOL5 Putative uncharacterized protein OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC) GN=ECP_2553 PE=4 SV=1
        sequence.setOriginalHeader(header);
        String[] data = header.split("|");

        if (data.length == 1) {
            throw new HeaderParseException("Header length parse error " + header);
        }
        if (data[0].equalsIgnoreCase("sp") || data[0].equalsIgnoreCase("tr")) {
            if (data[0].equalsIgnoreCase("sp")) {
                sequence.setAnnotationType(AnnotationType.CURATED);
            } else {
                sequence.setAnnotationType(AnnotationType.PREDICTED);
            }

            sequence.setAccession(new AccessionID(data[1], AccessionID.Source.UNIPROT));
            if (data.length > 1) {
                sequence.setDescription(data[2]);
            }

        } else if (data[0].equalsIgnoreCase("gi")) {
            AccessionID.Source giSource = AccessionID.Source.UNKNOWN;
            if (data.length >= 3) {
                if (data[2].equalsIgnoreCase("gb")) {
                    giSource = AccessionID.Source.GENBANK;
                } else if (data[2].equalsIgnoreCase("emb")) {
                    giSource = AccessionID.Source.EMBL;
                } else if (data[2].equalsIgnoreCase("dbj")) {
                    giSource = AccessionID.Source.DDBJ;
                }
                sequence.setAccession(new AccessionID(data[3]));
            } else {
                sequence.setAccession(new AccessionID(header, giSource));
            }
        } else if (data[0].equalsIgnoreCase("pir")) {
            sequence.setAccession(new AccessionID(data[2], AccessionID.Source.NBRF));
        } else if (data[0].equalsIgnoreCase("prf")) {
            sequence.setAccession(new AccessionID(data[2], AccessionID.Source.PRF));
        }else if (data[0].equalsIgnoreCase("pdb")) {
            sequence.setAccession(new AccessionID(data[1] + ":" + data[2], AccessionID.Source.PDB1));
        }else if (data[0].startsWith("PDB")){
            String[] pdbe = data[0].split(" ");
            String[] pdbaccession = pdbe[0].split(":");
            sequence.setAccession(new AccessionID(pdbaccession[1],AccessionID.Source.PDBe));
        }else if (data[0].indexOf(":") != -1){
            sequence.setAccession(new AccessionID(data[0],AccessionID.Source.PDB2));
        }else{
            sequence.setAccession(new AccessionID(header));
        }
    }
}
