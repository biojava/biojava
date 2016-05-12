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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Iterator;
import java.util.Set;

public class ParseGO {

	private static final Logger logger = LoggerFactory.getLogger(ParseGO.class);

	public static void main(String[] args){

		String u = "http://sourceforge.net/p/song/svn/HEAD/tree/trunk/subsets/biosapiens.obo?format=raw";

		try {
			URL url = new URL(u);

			OboParser parser = new OboParser();
			InputStream inStream = url.openStream();

			BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );

			Ontology ontology = parser.parseOBO(oboFile, "BioSapiens", "the BioSapiens ontology");

			Set<Term> keys = ontology.getTerms();
			Iterator<Term> iter = keys.iterator();
			while (iter.hasNext()){
				Term t = iter.next();
				logger.info("{} [{}]", t.getName(), t.getDescription());
			}
		} catch (Exception e){
			logger.error("Exception: ", e);
		}
	}
}
