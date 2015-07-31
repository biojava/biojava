/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io.blast;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Scanner;
import java.util.logging.Logger;
import org.biojava.nbio.core.exceptions.ParserException;
import org.biojava.nbio.core.search.io.Hit;
import org.biojava.nbio.core.search.io.Hsp;
import org.biojava.nbio.core.search.io.Result;
import org.biojava.nbio.core.search.io.ResultFactory;

/**
 *
 * @author Paolo Pavan, Genomnia srl
 * https://it.linkedin.com/pub/paolo-pavan/6/15a/956
 */
public class BlastTabularParser implements ResultFactory {
    private final String blastReference = 
            "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), A greedy algorithm for aligning DNA sequences&quot;, J Comput Biol 2000; 7(1-2):203-14.";
    /**
     * Tries to define a different level of consistency during parsing
     * LITERAL is intended a strict parsing much tight to the report.
     * IMPROVED consistency tries to import data much tight to the data model
     * (I hope you got the idea, if not have a look to the code.
     * I'm not very sure I will leave to the user the possibility to choose)
     */
    private enum PARSING_CONSISTENCY {
        IMPROVED,
        LITERAL
    }
    private static final Logger log = Logger.getLogger(BlastTabularParser.class.getName());
    
    private File targetFile;
    private int fileLinesCount;
    private PARSING_CONSISTENCY parsingConsistency = PARSING_CONSISTENCY.IMPROVED;
    
    
    // data imported private:
    int queryIdNumber = 0;
    HashMap<String,String> queryIdMapping = new HashMap();
    String programName=null, queryName = null, databaseFile = null;
    private String queryId      ;
    private String subjectId    ;
    private String percIdentity ;
    private String alnLength    ;
    private String mismatchCount;
    private String gapOpenCount ;
    private String queryStart   ;
    private String queryEnd     ;
    private String subjectStart ;
    private String subjectEnd   ;
    private String evalue       ;
    private String bitScore     ;
    

    @Override
    public List<String> getFileExtensions() {
        List l = new ArrayList();
        l.add("blasttabular");
        l.add("blasttxt");
        return l;
    }

    @Override
    public void setFile(File f) {
        targetFile = f;
    }

    @Override
    public List<Result> createObjects(double maxEScore) throws Exception {
        List<Result> results = new ArrayList();
        
        log.info("Query for hits");
        LineNumberReader  lnr = new LineNumberReader(new FileReader(targetFile));
        lnr.skip(Long.MAX_VALUE);
        fileLinesCount = lnr.getLineNumber();
        log.info(fileLinesCount + " hits approximately in all results");
        lnr.close();
        
        FileInputStream fileInputStream = new FileInputStream(targetFile);
        Scanner scanner = new Scanner(fileInputStream);
        
        String line = fetchData(scanner);
        while (scanner.hasNext()){
            try {
                BlastResultBuilder resultBuilder = new BlastResultBuilder();
                resultBuilder
                        .setQueryID(queryId)
                        .setDbFile(databaseFile)
                        .setProgram(programName)
                        .setQueryDef(queryName)
                        .setReference(blastReference);
                
                List<Hit> hits = new ArrayList();
                
                String currentQueryId = queryId;
                while (currentQueryId.equals(queryId) && scanner.hasNext()){
                    BlastHitBuilder hitBuilder = new BlastHitBuilder();
                    
                    List<Hsp> hsps = new ArrayList();
                    
                    String currentSubjectId=subjectId;
                    while (currentSubjectId.equals(subjectId) && scanner.hasNext()){
                        if (new Double(evalue) > maxEScore) {
                            line = fetchData(scanner);
                            continue;
                        }
                        BlastHspBuilder hspBuilder = new BlastHspBuilder();
                        hspBuilder
                            .setHspAlignLen(new Integer(alnLength))
                            .setHspGaps(new Integer(gapOpenCount))
                            .setHspQueryFrom(new Integer(queryStart))
                            .setHspQueryTo(new Integer(queryEnd))
                            .setHspHitFrom(new Integer(subjectStart))
                            .setHspHitTo(new Integer(subjectEnd))
                            .setHspEvalue(new Double(evalue))
                            .setHspBitScore(new Double(bitScore))
                            .setPercentageIdentity(new Double(percIdentity)/100)
                            .setMismatchCount(new Integer(mismatchCount));
                        hsps.add(hspBuilder.createBlastHsp());
                        line = fetchData(scanner);
                    }
                    hits.add(hitBuilder.setHsps(hsps).createBlastHit());
                }
                results.add(resultBuilder.setHits(hits).createBlastResult());
            } catch (NumberFormatException e) {
                throw new ParserException("Invalid numeric value met in:\n"+line);
            }
        }
        return results;
    }
    
    private String fetchData(Scanner scanner){
        String line;
        String[] split;
        
        line = scanner.nextLine();    
        while (line.startsWith("#")){
            // blast tabular with header options contains some more informations
            if (line.matches("#\\s.?BLAST.+")) programName = line.replace("#\\s","");
            if (line.startsWith("# Query:")) queryName = line.replace("# Query: ","");
            if (line.startsWith("# Database:")) databaseFile = line.replace("# Database: ","");
            
            // needed because blast report can end with a comment...
            if (!scanner.hasNext()) return null;
            line = scanner.nextLine(); 
        }
        
        // Here, programName != null checks if there was a header in the file
        boolean headerFound = programName != null;
        
        split = line.split("\\t");
        queryId      =split[0];
        subjectId    =split[1];
        percIdentity =split[2];
        alnLength    =split[3];
        mismatchCount=split[4];
        gapOpenCount =split[5];
        queryStart   =split[6];
        queryEnd     =split[7];
        subjectStart =split[8];
        subjectEnd   =split[9];
        evalue       =split[10];
        bitScore     =split[11];
        
        // blast tabular reports only the first word of the query name. 
        // If it was specified in the header it is better to use that definition
        if (parsingConsistency == PARSING_CONSISTENCY.IMPROVED && headerFound) {
            if (queryIdMapping.get(queryId)==null) {
                queryIdNumber ++;
                queryIdMapping.put(queryId,"Query_" + queryIdNumber);
            }
            // If a complete definition of the query name was readed, than we can use
            // a queryID schema that is consistent with blast xml report
            queryId = queryIdMapping.get(queryId);
        }
        if (!headerFound) queryName = queryId;
        
        return line;
    }

    @Override
    public void storeObjects(List<Result> results) throws Exception {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    /**
     * Intended for use with run module. 
     * Although possible does not make a lot of sense to have it with limited 
     * information in report
     * @param sequences 
     */
    @Override
    public void setQueryReferences(List sequences) {
        throw new UnsupportedOperationException("Not supported for this parser.");
    }
    /**
     * Intended for use with run module. 
     * Although possible does not make a lot of sense to have it with limited 
     * information in report
     * @param sequences 
     */
    @Override
    public void setDatabaseReferences(List sequences) {
        throw new UnsupportedOperationException("Not supported for this parser.");
    }

}
