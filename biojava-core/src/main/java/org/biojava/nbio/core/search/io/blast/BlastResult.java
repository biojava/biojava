/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.biojava.nbio.core.search.io.blast;

import org.biojava.nbio.core.search.io.Hit;
import org.biojava.nbio.core.search.io.Result;
import java.util.ArrayList;
import java.util.HashMap;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 *  This class models a Blast/Blast plus result.
 * You will find one of this for every query sequence specified in the run.
 * @author pavanpa
 */
public class BlastResult extends Result{
    public BlastResult(String program, String version, String reference, String dbFile, HashMap<String, String> programSpecificParameters, int iterationNumber, String queryID, String queryDef, int queryLength, ArrayList<Hit> hits, Sequence querySequence) {
        super(program, version, reference, dbFile, programSpecificParameters, iterationNumber, queryID, queryDef, queryLength, hits, querySequence);
    }
    
}



