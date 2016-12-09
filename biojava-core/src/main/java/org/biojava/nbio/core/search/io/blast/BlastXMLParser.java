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
package org.biojava.nbio.core.search.io.blast;


import org.biojava.nbio.core.search.io.Hit;
import org.biojava.nbio.core.search.io.Hsp;
import org.biojava.nbio.core.search.io.Result;
import org.biojava.nbio.core.search.io.ResultFactory;
import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.xpath.XPathException;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.core.util.XMLHelper;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

/**
 * Re-designed by Paolo Pavan on the footprint of: 
 * org.biojava.nbio.genome.query.BlastXMLQuery by Scooter Willis <willishf at gmail dot com>
 * 
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 *
 * @author Paolo Pavan
 */
public class BlastXMLParser <S extends Sequence<C>, C extends Compound> implements ResultFactory  <S,C>{
    private static final org.slf4j.Logger logger = LoggerFactory.getLogger(Hsp.class);
    Document blastDoc = null;
    private File targetFile;
    private List<S> queryReferences, databaseReferences;
    private Map<String,S> queryReferencesMap, databaseReferencesMap;
    
    public BlastXMLParser() {
        
    }
    @Override
    public void setFile(File f){
        targetFile = f;
    }
    
    private void readFile(String blastFile) throws IOException, ParseException{
        logger.info("Start reading " + blastFile);
        try {
            blastDoc = XMLHelper.loadXML(blastFile);
        } catch (SAXException ex) {
            logger.error("A parsing error has occurred while reading XML blast file");
            throw new ParseException(ex.getMessage(),0);
        } catch (ParserConfigurationException ex) {
            logger.error("Internal XML parser non properly configured");
            throw new ParseException(ex.getMessage(),0);
        }
        logger.info("Read finished");
    }
    
    @Override
    //@SuppressWarnings("unchecked")
    //@SuppressWarnings("rawtypes")
    public List<Result> createObjects(double maxEScore) throws IOException, ParseException {
        if (targetFile == null) throw new IllegalStateException("File to be parsed not specified.");
        
        // getAbsolutePath throws SecurityException
        readFile(targetFile.getAbsolutePath());
        // create mappings between sequences and blast id
        mapIds();
        
        ArrayList<Result> resultsCollection;
        ArrayList<Hit> hitsCollection;
        ArrayList<Hsp<S,C>> hspsCollection;
        
        try {
            // select top level elements
            String program = XMLHelper.selectSingleElement(blastDoc.getDocumentElement(),"BlastOutput_program").getTextContent();
            String version = XMLHelper.selectSingleElement(blastDoc.getDocumentElement(),"BlastOutput_version").getTextContent();
            String reference = XMLHelper.selectSingleElement(blastDoc.getDocumentElement(),"BlastOutput_reference").getTextContent();
            String dbFile = XMLHelper.selectSingleElement(blastDoc.getDocumentElement(),"BlastOutput_db").getTextContent();
            
            // try to guess the sequence type of the blast run
            // the most reliable way I can see is to exploit the program blast used
            Boolean isDNA = null;
            if (program.endsWith("blastp")) isDNA = false;
            if (program.endsWith("blastn")) isDNA = true;

            logger.info("Query for hits in "+ targetFile);
            ArrayList<Element> IterationsList = XMLHelper.selectElements(blastDoc.getDocumentElement(), "BlastOutput_iterations/Iteration[Iteration_hits]");
            logger.info(IterationsList.size() + " results");

            resultsCollection = new ArrayList<Result>();
            for (Element element : IterationsList) {
                BlastResultBuilder resultBuilder = new BlastResultBuilder();
                // will add BlastOutput* key sections in the result object
                resultBuilder
                    .setProgram(program)
                    .setVersion(version)
                    .setReference(reference)
                    .setDbFile(dbFile);

                // Iteration* section keys:
                resultBuilder
                    .setIterationNumber(new Integer(XMLHelper.selectSingleElement(element,"Iteration_iter-num").getTextContent()))
                    .setQueryID(XMLHelper.selectSingleElement(element,"Iteration_query-ID").getTextContent())
                    .setQueryDef(XMLHelper.selectSingleElement(element, "Iteration_query-def").getTextContent())
                    .setQueryLength(new Integer(XMLHelper.selectSingleElement(element,"Iteration_query-len").getTextContent()));
                
                if (queryReferences != null) resultBuilder.setQuerySequence(queryReferencesMap.get(
                        XMLHelper.selectSingleElement(element,"Iteration_query-ID").getTextContent()
                ));



                Element iterationHitsElement = XMLHelper.selectSingleElement(element, "Iteration_hits");
                ArrayList<Element> hitList = XMLHelper.selectElements(iterationHitsElement, "Hit");

                hitsCollection = new ArrayList<Hit>();
                for (Element hitElement : hitList) {
                    // deliberately not parametrized
                    BlastHitBuilder blastHitBuilder = new BlastHitBuilder();
                    blastHitBuilder
                        .setHitNum(new Integer(XMLHelper.selectSingleElement(hitElement, "Hit_num").getTextContent()))
                        .setHitId(XMLHelper.selectSingleElement(hitElement, "Hit_id").getTextContent())
                        .setHitDef(XMLHelper.selectSingleElement(hitElement, "Hit_def").getTextContent())
                        .setHitAccession(XMLHelper.selectSingleElement(hitElement, "Hit_accession").getTextContent())
                        .setHitLen(new Integer(XMLHelper.selectSingleElement(hitElement, "Hit_len").getTextContent()));
                    
                    // deliberately not parametrized (see above)
                    if (databaseReferences != null) blastHitBuilder.setHitSequence(databaseReferencesMap.get(
                        XMLHelper.selectSingleElement(hitElement, "Hit_id").getTextContent()
                    ));

                    Element hithspsElement = XMLHelper.selectSingleElement(hitElement, "Hit_hsps");
                    ArrayList<Element> hspList = XMLHelper.selectElements(hithspsElement, "Hsp");

                    hspsCollection = new ArrayList<Hsp<S,C>>();
                    for (Element hspElement : hspList) {
                        Double evalue = new Double(XMLHelper.selectSingleElement(hspElement, "Hsp_evalue").getTextContent());

                        // add the new hsp only if it pass the specified threshold. It can save lot of memory and some parsing time
                        if (evalue <= maxEScore) {
                            int hspNum = new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_num").getTextContent());
                            String hspQuerySeq = XMLHelper.selectSingleElement(hspElement, "Hsp_qseq").getTextContent();
                            String hspHitSeq = XMLHelper.selectSingleElement(hspElement, "Hsp_hseq").getTextContent();
                            BlastHspBuilder<S,C> blastHspBuilder = new BlastHspBuilder<S,C>();
                            blastHspBuilder
                                .setHspNum(hspNum)
                                .setHspBitScore(new Double(XMLHelper.selectSingleElement(hspElement, "Hsp_bit-score").getTextContent()))
                                .setHspScore(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_score").getTextContent()))
                                .setHspEvalue(evalue)
                                .setHspQueryFrom(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_query-from").getTextContent()))
                                .setHspQueryTo(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_query-to").getTextContent()))
                                .setHspHitFrom(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_hit-from").getTextContent()))
                                .setHspHitTo(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_hit-to").getTextContent()))
                                .setHspQueryFrame(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_query-frame").getTextContent()))
                                .setHspHitFrame(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_hit-frame").getTextContent()))
                                .setHspIdentity(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_identity").getTextContent()))
                                .setHspPositive(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_positive").getTextContent()))
                                .setHspGaps(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_gaps").getTextContent()))
                                .setHspAlignLen(new Integer(XMLHelper.selectSingleElement(hspElement, "Hsp_align-len").getTextContent()))
                                .setHspQseq(hspQuerySeq)
                                .setHspHseq(hspHitSeq)
                                .setHspIdentityString(XMLHelper.selectSingleElement(hspElement, "Hsp_midline").getTextContent());
                            
                            S s;
                            CompoundSet c;
                            String querySequenceString = hspQuerySeq.replace("-", "");
                            try {
                                if (isDNA) {
                                    c = DNACompoundSet.getDNACompoundSet();
                                    s = (S)new DNASequence(hspQuerySeq, c);
                                } else {
                                    c = AminoAcidCompoundSet.getAminoAcidCompoundSet();
                                    s = (S)new ProteinSequence(hspQuerySeq, c);
                                }
                                
                                blastHspBuilder.setSequence(s);
                                //blastHspBuilder.setCompound(c); // ???????
                            } catch (CompoundNotFoundException e){
                                logger.error("Compoundset not found error for HSP num "+hspNum
                                        +". It will not be set by blast parser and I will try to continue.");
                            }

                            hspsCollection.add(blastHspBuilder.createBlastHsp());
                        }
                    }
                    // finally set the computed hsp collection and create Hit object
                    // (deliberately not parametrized, see above)
                    blastHitBuilder.setHsps(hspsCollection);
                    hitsCollection.add(blastHitBuilder.createBlastHit());
                }
                // finally set the computed Hit collection to the result
                resultBuilder.setHits(hitsCollection);
                resultsCollection.add(resultBuilder.createBlastResult());
            }
        } catch (XPathException e) {
            throw new ParseException(e.getMessage(),0);
        }
        logger.info("Parsing of "+targetFile+" finished.");
        
        return resultsCollection;
    }
    
    @Override
    public List<String> getFileExtensions(){
        ArrayList<String> extensions = new ArrayList<String>(1);
        extensions.add("blastxml");
        return extensions;
    }

    @Override
    public void setQueryReferences(List<S> sequences) {
        queryReferences = sequences;
    }

    @Override
    public void setDatabaseReferences(List<S> sequences) {
        databaseReferences = sequences;
    }
    
    /**
     * fill the map association between sequences an a unique id
     */
    private void mapIds() {
        if (queryReferences != null) {
            queryReferencesMap = new HashMap<String,S>(queryReferences.size());
            for (int counter=0; counter < queryReferences.size() ; counter ++){
                String id = "Query_"+(counter+1);
                queryReferencesMap.put(id, queryReferences.get(counter));
            }  
        }
        
        if (databaseReferences != null) {
            databaseReferencesMap = new HashMap<String,S>(databaseReferences.size());
            for (int counter=0; counter < databaseReferences.size() ; counter ++){
                // this is strange: while Query_id are 1 based, Hit (database) id are 0 based
                String id = "gnl|BL_ORD_ID|"+(counter);
                databaseReferencesMap.put(id, databaseReferences.get(counter));
            }  
        }
    }

    @Override
    public void storeObjects(List<Result> results) throws IOException, ParseException {
        throw new UnsupportedOperationException("This parser does not support writing yet.");
    }
}


class BlastHsp<S extends Sequence<C>, C extends Compound> extends org.biojava.nbio.core.search.io.Hsp<S,C> {
    public BlastHsp(int hspNum, double hspBitScore, int hspScore, double hspEvalue, int hspQueryFrom, int hspQueryTo, int hspHitFrom, int hspHitTo, int hspQueryFrame, int hspHitFrame, int hspIdentity, int hspPositive, int hspGaps, int hspAlignLen, String hspQseq, String hspHseq, String hspIdentityString, Double percentageIdentity, Integer mismatchCount, S s, C c) {
        super(hspNum, hspBitScore, hspScore, hspEvalue, hspQueryFrom, hspQueryTo, hspHitFrom, hspHitTo, hspQueryFrame, hspHitFrame, hspIdentity, hspPositive, hspGaps, hspAlignLen, hspQseq, hspHseq, hspIdentityString, percentageIdentity, mismatchCount, s, c);
    }
    
}

class BlastHit extends org.biojava.nbio.core.search.io.Hit {
    public BlastHit(int hitNum, String hitId, String hitDef, String hitAccession, int hitLen, List<Hsp<?,?>> hitHsps, Sequence<?> hitSequence) {
        super(hitNum, hitId, hitDef, hitAccession, hitLen, hitHsps, hitSequence);
    }
    
}