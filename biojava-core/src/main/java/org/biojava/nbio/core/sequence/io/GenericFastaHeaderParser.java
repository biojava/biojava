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
 * Created on 01-21-2010
 */
package org.biojava.nbio.core.sequence.io;

import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DataSource;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.io.template.SequenceHeaderParserInterface;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.AbstractSequence.AnnotationType;
import org.biojava.nbio.core.sequence.template.Compound;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

/**
 * The default fasta header parser where some headers are well defined based on the source
 * database which allows us to set the source of the protein sequence and the identifier
 * that can be used in future implementations to load features from external sources
 *
 * If the user has a custom header with local data then they can create their own implementation
 * of a FastaHeaderParserInterface
 *<pre>
 * GenBank                           gi|gi-number|gb|accession|locus
 * ENA Data Library                  gi|gi-number|emb|accession|locus
 * DDBJ, DNA Database of Japan       gi|gi-number|dbj|accession|locus
 * NBRF PIR                          pir||entry
 * Protein Research Foundation       prf||name
 * SWISS-PROT                        sp|accession|name
 * Brookhaven Protein Data Bank (1)  pdb|entry|chain
 * Brookhaven Protein Data Bank (2)  entry:chain|PDBID|CHAIN|SEQUENCE
 * PDB EBI                           PDB:1ECY_A mol:protein length:142  ECOTIN
 * Patents                           pat|country|number
 * GenInfo Backbone Id               bbs|number
 * General database identifier       gnl|database|identifier
 * NCBI Reference Sequence           ref|accession|locus
 * Local Sequence identifier         lcl|identifier
 *</pre>
 * @author Scooter Willis 
 */
public class GenericFastaHeaderParser<S extends AbstractSequence<C>, C extends Compound> implements SequenceHeaderParserInterface<S,C> {

	private final static Logger logger = LoggerFactory.getLogger(GenericFastaHeaderParser.class);

	/**
	 * Parse out the components where some have a | and others do not
	 * @param header
	 * @return
	 */
	private String[] getHeaderValues(String header) {
		String[] data = new String[0];
		List<String> values = new ArrayList<>();
		StringBuffer sb = new StringBuffer();
		//commented out 1/11/2012 to resolve an issue where headers do contain a length= at the end that are not recognized
		//if(header.indexOf("length=") != -1){
		//    data = new String[1];
		//    int index = header.indexOf("length=");
		//    data[0] = header.substring(0, index).trim();
	//        logger.debug("accession=" + data[0]);
		//    return data;
		//} else
		 if (!header.startsWith("PDB:")) {
			for (int i = 0; i < header.length(); i++) {
				if (header.charAt(i) == '|') {
					values.add(sb.toString());
					sb.setLength(0);//faster than  = new StringBuffer();
				} else if (i == header.length() - 1) {
					sb.append(header.charAt(i));
					values.add(sb.toString());
				} else {
					sb.append(header.charAt(i));
				}

			}
			data = new String[values.size()];
			values.toArray(data);
		} else {
			data = header.split(" ");
		}
		return data;
	}

	/**
	 * Parse the header and set the values in the sequence
	 * @param header
	 * @param sequence
	 */
	@Override
	public void parseHeader(String header, S sequence) {
		//uniptrot
		// tr|Q0TET7|Q0TET7_ECOL5 Putative uncharacterized protein OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC) GN=ECP_2553 PE=4 SV=1
		sequence.setOriginalHeader(header);
		String[] data = getHeaderValues(header);

		if (data.length == 1) {
			sequence.setAccession(new AccessionID(data[0]));
		} else  if ("sp".equalsIgnoreCase(data[0]) || "tr".equalsIgnoreCase(data[0])) {
			if ("sp".equalsIgnoreCase(data[0])) {
				sequence.setAnnotationType(AnnotationType.CURATED);
			} else {
				sequence.setAnnotationType(AnnotationType.PREDICTED);
			}

			sequence.setAccession(new AccessionID(data[1], DataSource.UNIPROT));
			if (data.length > 2) {
				sequence.setDescription(data[2]);
			}

		} else if ("gi".equalsIgnoreCase(data[0])) {
			DataSource giSource = DataSource.UNKNOWN;
			if (data.length >= 3) {
				if ("gb".equalsIgnoreCase(data[2])) {
					giSource = DataSource.GENBANK;
				} else if ("emb".equalsIgnoreCase(data[2])) {
					giSource = DataSource.ENA;
				} else if ("dbj".equalsIgnoreCase(data[2])) {
					giSource = DataSource.DDBJ;
				}
				sequence.setAccession(new AccessionID(data[3], giSource));
			} else {
				sequence.setAccession(new AccessionID(header, giSource));
			}
		} else if ("pir".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[2], DataSource.NBRF));
		} else if ("prf".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[2], DataSource.PRF));
		} else if ("pdb".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[1] + ":" + data[2], DataSource.PDB1));
		} else if (data[0].startsWith("PDB")) {
			String[] pdbe = data[0].split(" ");
			String[] pdbaccession = pdbe[0].split(":");
			sequence.setAccession(new AccessionID(pdbaccession[1], DataSource.PDBe));
		} else if (data[0].indexOf(":") != -1 && data.length > 1 && "PDBID".equals(data[1])) {
			sequence.setAccession(new AccessionID(data[0], DataSource.PDB2));
		} else if ("pat".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[2], DataSource.PATENTS));
		} else if ("bbs".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[1], DataSource.GENINFO));
		} else if ("gnl".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[2], DataSource.GENERAL));
		} else if ("ref".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[1], DataSource.NCBI));
		} else if ("lcl".equalsIgnoreCase(data[0])) {
			sequence.setAccession(new AccessionID(data[1], DataSource.LOCAL));
		} else {
			sequence.setAccession(new AccessionID(data[0])); // avoid the common problem of picking up all the comments original header in getOriginalHeader
		}


	}

	
}
