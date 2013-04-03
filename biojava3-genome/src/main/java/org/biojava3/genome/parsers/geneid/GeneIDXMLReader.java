/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.genome.parsers.geneid;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.logging.Logger;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.io.FastaWriterHelper;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.io.FastaWriterHelper;
import org.biojava3.core.util.XMLHelper;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class GeneIDXMLReader {

    private static final Logger log = Logger.getLogger(GeneIDXMLReader.class.getName());
    Document geneidDoc = null;

    public GeneIDXMLReader(String geneidXMLFile) throws Exception {
        log.info("Start read of " + geneidXMLFile);
        geneidDoc = XMLHelper.loadXML(geneidXMLFile);
        log.info("Read finished");
    }

    public LinkedHashMap<String, ProteinSequence> getProteinSequences() throws Exception {
        LinkedHashMap<String, ProteinSequence> proteinSequenceList = new LinkedHashMap<String, ProteinSequence>();
        ArrayList<Element> elementList = XMLHelper.selectElements(geneidDoc.getDocumentElement(), "prediction/gene/protein");
        log.info(elementList.size() + " hits");

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
        log.info(elementList.size() + " hits");

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
            e.printStackTrace();
        }
    }
}
