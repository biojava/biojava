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
 */

package org.biojava.nbio.ontology.io;

import org.biojava.nbio.ontology.*;
import org.biojava.nbio.ontology.obo.OboFileEventListener;
import org.biojava.nbio.ontology.obo.OboFileHandler;
import org.biojava.nbio.ontology.obo.OboFileParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.text.ParseException;




/** Parses an OBO file.
 *
 * @author Andreas Prlic
 * @since 1.7
 *
 * <h2>Example</h2>
 * <pre>
 * OboParser parser = new OboParser();
		InputStream inStream = this.getClass().getResourceAsStream("/files/ontology/biosapiens.obo");

		BufferedReader oboFile = new BufferedReader ( new InputStreamReader ( inStream ) );
		try {
			Ontology ontology = parser.parseOBO(oboFile, "BioSapiens", "the BioSapiens ontology");

			Set keys = ontology.getTerms();
			Iterator iter = keys.iterator();
			while (iter.hasNext()){
				System.out.println(iter.next());
			}

		} catch (Exception e){
			e.printStackTrace();
		}
 * </pre>
 *
 */
public class OboParser {

	/** Parse a OBO file and return its content as a BioJava Ontology object
	 *
	 * @param oboFile the file to be parsed
	 * @param ontoName
	 * @param ontoDescription

	 * @return the ontology represented as a BioJava ontology file
	 * @throws ParseException
	 * @throws IOException
	 */
	public Ontology parseOBO(
			BufferedReader oboFile,
			String ontoName,
			String ontoDescription
			)
					throws ParseException, IOException {

		try {
			OntologyFactory factory = OntoTools.getDefaultFactory();
			Ontology ontology = factory.createOntology(ontoName, ontoDescription);

			OboFileParser parser = new OboFileParser();

			OboFileEventListener handler = new OboFileHandler(ontology);

			parser.addOboFileEventListener(handler);
			parser.parseOBO(oboFile);

			return ontology;


		} catch (AlreadyExistsException ex) {
			throw new RuntimeException( "Duplication in ontology");
		} catch (OntologyException ex) {
			throw new RuntimeException(ex);
		}

	}
}
