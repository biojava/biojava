/*
 *                  BioJava development code
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
 * Created on Jan 18, 2008
 *
 * since 1.6
 */

package org.biojava.nbio.ontology.obo;

import org.biojava.nbio.ontology.Synonym;

/** an interface for events that occur during parsing of .obo files
 *
 * @author Andreas Prlic
 *
 */
public interface OboFileEventListener {

	/** starting to parse a new OBO file
	 *
	 *
	 */
	public void documentStart();

	/** end of parsing a new OBO file
	 *
	 *
	 */
	public void documentEnd();

	/** parsed a new OBO file header
	 *
	 *
	 */
	public void newOboFileHeader();

	/** parsed a new stanza in the file
	 *
	 * @param stanza
	 */
	public void newStanza(String stanza);

	/**found a new key in the file
	 *
	 * @param key
	 * @param value
	 */
	public void newKey(String key, String value );

	/** a new synonym has been found
	 *
	 * @param synonym
	 */
	public void newSynonym(Synonym synonym);
}
