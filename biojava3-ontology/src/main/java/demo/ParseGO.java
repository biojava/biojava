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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.Iterator;
import java.util.Set;

import org.biojava3.ontology.Ontology;
import org.biojava3.ontology.Term;
import org.biojava3.ontology.io.OboParser;

public class ParseGO {
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
				System.out.println(t.getName() + " [" + t.getDescription() + "]");
			}


		} catch (Exception e){
			e.printStackTrace();
		}

	}
}

