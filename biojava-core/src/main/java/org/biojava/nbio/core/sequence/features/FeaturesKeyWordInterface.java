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
 * Created on 01-21-2010
 */

package org.biojava.nbio.core.sequence.features;

import java.util.List;

/**
 * Models the keywords that are annotated for a protein sequence at Uniprot. If a ProxySequenceReader
 * implements this interface then the sequence will call this method
 *
 * @author Scooter Willis 
 */
public interface FeaturesKeyWordInterface {

	/**
	 *
	 * @return
	 */
	public List<String> getKeyWords() ;
}
