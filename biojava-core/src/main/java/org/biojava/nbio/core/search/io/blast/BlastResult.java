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
import org.biojava.nbio.core.search.io.Result;
import java.util.HashMap;
import java.util.List;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * This class models a Blast/Blast plus result.
 * You will find one of this for every query sequence specified in the run.
 * 
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info 
 * or discuss major changes.
 * https://github.com/paolopavan
 * 
 * @author Paolo Pavan
 * 
 */
public class BlastResult extends Result{
    public BlastResult(String program, String version, String reference, String dbFile, HashMap<String, String> programSpecificParameters, int iterationNumber, String queryID, String queryDef, int queryLength, List<Hit> hits, Sequence querySequence) {
        super(program, version, reference, dbFile, programSpecificParameters, iterationNumber, queryID, queryDef, queryLength, hits, querySequence);
    }
    
}



