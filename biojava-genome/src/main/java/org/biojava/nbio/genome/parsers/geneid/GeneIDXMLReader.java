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
package org.biojava.nbio.genome.parsers.geneid;

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.io.FastaWriterHelper;
import org.biojava.nbio.core.util.XMLHelper;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GeneIDXMLReader {

	private static final Logger logger = LoggerFactory.getLogger(GeneIDXMLReader.class);

	Document geneidDoc = null;

	public GeneIDXMLReader(String geneidXMLFile) throws Exception {
		logger.info("Start read of {}", geneidXMLFile);
		geneidDoc = XMLHelper.loadXML(geneidXMLFile);
		logger.info("Read finished");
	}

	public LinkedHashMap<String, ProteinSequence> getProteinSequences() throws Exception {
		LinkedHashMap<String, ProteinSequence> proteinSequenceList = new LinkedHashMap<String, ProteinSequence>();
		ArrayList<Element> elementList = XMLHelper.selectElements(geneidDoc.getDocumentElement(), "prediction/gene/protein");
		logger.info("{} hits", elementList.size());

		for (Element proteinElement : elementList) {
			Element geneElement = (Element) proteinElement.getParentNode();
			String sequence = proteinElement.getTextContent().replaceAll("\\W","");
			ProteinSequence proteinSequence = new ProteinSequence(sequence);
			String idGene = geneElement.getAttribute("idGene");
			proteinSequence.setAccession(new AccessionID(idGene));
			proteinSequenceList.put(idGene, proteinSequence);
		}

		return proteinSequenceList;
	}

	public LinkedHashMap<String, DNASequence> getDNACodingSequences() throws Exception {
		LinkedHashMap<String, DNASequence> dnaSequenceList = new LinkedHashMap<String, DNASequence>();
		ArrayList<Element> elementList = XMLHelper.selectElements(geneidDoc.getDocumentElement(), "prediction/gene/cDNA");
		logger.info("{} hits", elementList.size());

		for (Element dnaElement : elementList) {
			Element geneElement = (Element) dnaElement.getParentNode();
			String sequence = dnaElement.getTextContent().replaceAll("\\W","");
			DNASequence dnaSequence = new DNASequence(sequence);
			String idGene = geneElement.getAttribute("idGene");
			dnaSequence.setAccession(new AccessionID(idGene));
			dnaSequenceList.put(idGene, dnaSequence);
		}

		return dnaSequenceList;
	}

	public static void main(String[] args) {
		try {
			GeneIDXMLReader geneIDXMLReader = new GeneIDXMLReader("/Users/Scooter/scripps/dyadic/geneid/geneid/c1_geneid.xml");
			LinkedHashMap<String, ProteinSequence> proteinSequenceHashMap = geneIDXMLReader.getProteinSequences();
			FastaWriterHelper.writeProteinSequence(new File("/Users/Scooter/scripps/dyadic/geneid/geneid/c1_geneid.faa"), proteinSequenceHashMap.values());

			LinkedHashMap<String, DNASequence> dnaSequenceHashMap = geneIDXMLReader.getDNACodingSequences();
			FastaWriterHelper.writeNucleotideSequence(new File("/Users/Scooter/scripps/dyadic/geneid/geneid/c1_geneid.fna"), dnaSequenceHashMap.values());

		} catch (Exception e) {
			logger.error("Exception: ", e);
		}
	}
}
