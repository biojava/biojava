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
package org.biojava.nbio.core.search.io;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * This class models a search result.
 * You will find one of this for every query sequence specified in the run.
 * 
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
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
    private List<Hit> hits;
    private int hitCounter = -1;

    public Result(String program, String version, String reference, String dbFile, HashMap<String, String> programSpecificParameters, int iterationNumber, String queryID, String queryDef, int queryLength, List<Hit> hits, Sequence querySequence) {
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

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + (this.queryID != null ? this.queryID.hashCode() : 0);
        hash = 29 * hash + (this.queryDef != null ? this.queryDef.hashCode() : 0);
        hash = 29 * hash + (this.hits != null ? this.hits.hashCode() : 0);
        return hash;
    }
    /**
     * Experimental.
     * Wants to return an hashcode designed to allow conceptual comparisons of search results.
     * Wants to implement conceptual comparisons of search results.
     * Fields unrelated to search are deliberately not considered.
     * @return 
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Result other = (Result) obj;
        if ((this.queryID == null) ? (other.queryID != null) : !this.queryID.equals(other.queryID)) {
            return false;
        }
        if ((this.queryDef == null) ? (other.queryDef != null) : !this.queryDef.equals(other.queryDef)) {
            return false;
        }
        if (this.hits != other.hits && (this.hits == null || !this.hits.equals(other.hits))) {
            return false;
        }
        return true;
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
     * returns the reference to the original and whole sequence used to query the database.
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
