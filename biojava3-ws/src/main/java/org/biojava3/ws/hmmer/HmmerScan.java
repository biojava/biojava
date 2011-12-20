package org.biojava3.ws.hmmer;

import java.io.IOException;

import java.util.SortedSet;

import org.biojava3.core.sequence.ProteinSequence;

/** Interface for performing Hmmscans on sequences.
 * 
 * @author Andreas Prlic
 *
 */
public interface HmmerScan {

	public  SortedSet<HmmerResult> scan(ProteinSequence sequence) throws IOException;
	
}
