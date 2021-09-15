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
 * Created on Mar 19, 2014
 * Author: andreas
 *
 */

package demo;

import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.io.OboParser;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import java.util.Iterator;
import java.util.Set;

public class ParseGO {

	/**
	 * Parses Biosapiens OBO file and prints out name/description
	 * pairs
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		try {

			OboParser parser = new OboParser();
			InputStream inStream = OboParser.class.getResourceAsStream("/ontology/biosapiens.obo");

			BufferedReader oboFile = new BufferedReader(new InputStreamReader(inStream));

			Ontology ontology = parser.parseOBO(oboFile, "BioSapiens", "the BioSapiens ontology");

			Set<Term> keys = ontology.getTerms();
			Iterator<Term> iter = keys.iterator();
			while (iter.hasNext()) {
				Term t = iter.next();
				System.out.println(t.getName() + " " + t.getDescription());
			}
		} catch (Exception e) {
			System.err.println("Exception: " + e);
		}
	}
}
