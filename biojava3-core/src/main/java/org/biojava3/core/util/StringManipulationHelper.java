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
 * Created on Sep 14, 2011
 * Author: Amr AL-Hossary 
 *
 */
package org.biojava3.core.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;


/**
 * A utility class for common {@link String} manipulation tasks.
 * All functions are static methods.
 * 
 * @author Amr AL-Hossary
 */
public class StringManipulationHelper  {

	/**
	 * we are using Unix endline here, since this is used for testing XML and it
	 * is part of the XML recommendations: <a href
	 * ="http://www.w3.org/TR/REC-xml/#sec-line-ends"
	 * >http://www.w3.org/TR/REC-xml/#sec-line-ends</a>
	 */
	private static final String UNIX_NEWLINE = "\n";

	private StringManipulationHelper() {
		// to prevent instantiation
	}

	

	

	public static String convertStreamToString(InputStream stream) {
		BufferedReader reader = new BufferedReader(new InputStreamReader(stream));
		StringBuilder sb = new StringBuilder();

		String line = null;
		try {
			while ((line = reader.readLine()) != null) {

				sb.append(line + UNIX_NEWLINE);
			}
		} catch (IOException e) {
			// e.printStackTrace();
		} finally {
			try {
				stream.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return sb.toString();
	}

}
