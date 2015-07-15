/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 *  This class models a search result.
 * You will find one of this for every query sequence specified in the run.
 * @author pavanpa
 */
public abstract class Result implements Iterable<Hit>{
    private String program;
    private String version;
    private String reference;
    private String dbFile;

    private HashMap<String,String> programSpecificParameters;
    
    private int iterationNumber;
    private String queryID;
    private String queryDef;
    private int queryLength;
    private Sequence querySequence;
    private ArrayList<Hit> hits;
    private int hitCounter = -1;
  
    /*
    public Hit nextHit() {
        return(hits.get(hitCounter++));
    }
    */

    public Result(String program, String version, String reference, String dbFile, HashMap<String, String> programSpecificParameters, int iterationNumber, String queryID, String queryDef, int queryLength, ArrayList<Hit> hits, Sequence querySequence) {
        this.program = program;
        this.version = version;
        this.reference = reference;
        this.dbFile = dbFile;
        this.programSpecificParameters = programSpecificParameters;
        this.iterationNumber = iterationNumber;
        this.queryID = queryID;
        this.queryDef = queryDef;
        this.queryLength = queryLength;
        this.hits = hits;
        this.querySequence = querySequence;
    }

    public int getIterationNumber() {
        return iterationNumber;
    }

    public String getQueryID() {
        return queryID;
    }

    public String getQueryDef() {
        return queryDef;
    }

    public int getQueryLength() {
        return queryLength;
    }

    public int getHitCounter() {
        return hitCounter;
    }
    
    public String getProgram() {
        return program;
    }

    public String getVersion() {
        return version;
    }

    public String getReference() {
        return reference;
    }

    public String getDbFile() {
        return dbFile;
    }

    public Set<String> getProgramSpecificParametersList() {
        return programSpecificParameters.keySet();
    }
    
    public String getProgramSpecificParameter(String key) {
        return programSpecificParameters.get(key);
    }
    /**
     * returns the reference th the original and whole sequence used to query the database.
     * Available only if the ResultFactory implements setQueryReferences and
     * it was used before the parsing with SearchIO
     * @return Sequence object
     */
    public Sequence getQuerySequence() {
        return querySequence;
    }
    
    @Override
    public Iterator<Hit> iterator() {
        return new Iterator<Hit>() {
            int currentResult = 0;
            @Override
            public boolean hasNext() {
                return currentResult < hits.size();
            }

            @Override
            public Hit next() {
                return hits.get(currentResult++);
            }

            @Override
            public void remove() {
                throw new UnsupportedOperationException("The remove operation is not supported by this iterator");
            }
        };
    }
}
