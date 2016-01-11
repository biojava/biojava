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
import java.util.HashMap;
import java.util.List;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 */

public class BlastResultBuilder {
    private String program;
    private String version;
    private String reference;
    private String dbFile;
    private HashMap<String, String> programSpecificParameters;
    private int iterationNumber;
    private String queryID;
    private String queryDef;
    private int queryLength;
    private Sequence querySequence;
    private List<Hit> hits;

    public BlastResultBuilder() {
    }

    public BlastResultBuilder setProgram(String program) {
        this.program = program;
        return this;
    }

    public BlastResultBuilder setVersion(String version) {
        this.version = version;
        return this;
    }

    public BlastResultBuilder setReference(String reference) {
        this.reference = reference;
        return this;
    }

    public BlastResultBuilder setDbFile(String dbFile) {
        this.dbFile = dbFile;
        return this;
    }

    public BlastResultBuilder setProgramSpecificParameters(HashMap<String, String> programSpecificParameters) {
        this.programSpecificParameters = programSpecificParameters;
        return this;
    }

    public BlastResultBuilder setIterationNumber(int iterationNumber) {
        this.iterationNumber = iterationNumber;
        return this;
    }

    public BlastResultBuilder setQueryID(String queryID) {
        this.queryID = queryID;
        return this;
    }

    public BlastResultBuilder setQueryDef(String queryDef) {
        this.queryDef = queryDef;
        return this;
    }

    public BlastResultBuilder setQueryLength(int queryLength) {
        this.queryLength = queryLength;
        return this;
    }

    public BlastResultBuilder setHits(List<Hit> hits) {
        this.hits = hits;
        return this;
    }
    
    public BlastResultBuilder setQuerySequence(Sequence s) {
        this.querySequence = s;
        return this;
    }

    public BlastResult createBlastResult() {
        return new BlastResult(program, version, reference, dbFile, programSpecificParameters, iterationNumber, queryID, queryDef, queryLength, hits, querySequence);
    }

}
