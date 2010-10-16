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
 * created at Mar 4, 2008
 */
package org.biojava.bio.structure.io.mmcif;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;

/** Interface that needs to be implemented by an MMcifParser
 *
 * @author Andreas Prlic
 * @since 1.7
 */
public interface MMcifParser {

	/** Add a MMcifConsumer that listens to even being triggered by the parser and processes the data into a backend provided by the Consumer.
	 *
	 * @param consumer a consumer object.
	 */
	public void addMMcifConsumer(MMcifConsumer consumer);

	/** Remove all consumers from the parser.
	 *
	 */
	public void clearConsumers();

	/** remove a single consumer from the parser
	 *
	 * @param consumer
	 */
	public void removeMMcifConsumer(MMcifConsumer consumer);


	/** Start the actual parsing. The parser will trigger events that are defined by the MMcifConsumer class.
	 *
	 * @param buf a BufferedReader.
	 */
	public void parse(BufferedReader buf) throws IOException;

	/** Start the actual parsing. The parser will trigger events that are defined by the MMcifConsumer class.
	 *
	 * @param inStream InputStream to parse from.
	 */
	public void parse(InputStream inStream) throws IOException;


}
