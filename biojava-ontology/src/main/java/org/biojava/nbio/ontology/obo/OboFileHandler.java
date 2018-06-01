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

package org.biojava.nbio.ontology.obo;

import org.biojava.nbio.ontology.AlreadyExistsException;
import org.biojava.nbio.ontology.Ontology;
import org.biojava.nbio.ontology.Synonym;
import org.biojava.nbio.ontology.Term;
import org.biojava.nbio.ontology.utils.Annotation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;

/** A  file handler for .obo files
 *
 * @author Andreas Prlic
 *
 */
public class OboFileHandler implements OboFileEventListener {

	private static final Logger logger = LoggerFactory.getLogger(OboFileEventListener.class);

	Ontology ontology;
	List<Term> termStack ;

	public static final String TERM        = "Term";
	public static final String TYPEDEF     = "Typedef";
	public static final String ONTOLOGY    = "ontologys";
	public static final String ID_KEY      = "id";
	public static final String SYNONYM     = "synonym";
	public static final String EXACT_SYNONYM  = "exact_synonym";
	public static final String BROAD_SYNONYM  = "broad_synonym";
	public static final String NARROW_SYNONYM = "narrow_synonym";
	public static final String REL_SYNONYM = "related_synonym";
	public static final String NAME        = "name";
	public static final String DEF         = "def";
	public static final String XREF_ANALOG = "xref_analog";
	public static final String COMMENT     = "comment";
	public static final String IS_A        = "is_a";
	public static final String IS_OBSOLETE = "is_obsolete";
	public static final String RELATIONSHIP = "relationship";
	public static final String DISJOINT_FROM = "disjoint_from";
	public static final String SUBSET       = "subset";
	public static final String INTERSECTION_OF = "intersection_of";
	public static final String NAMESPACE = "namespace";
	public static final String REPLACED_BY = "replaced_by";


	public static final String ALT_ID      = "alt_id";

	boolean isTerm ;

	private Term currentTerm;

	public OboFileHandler(Ontology ontology){
		this.ontology = ontology ;

		//Term isa = onto.importTerm(OntoTools.IS_A, null);
		//Term partof = onto.importTerm(OntoTools.PART_OF, null);;

	}

	@Override
	public void documentEnd() {
		// TODO Auto-generated method stub

	}

	@Override
	public void documentStart() {
		termStack = new ArrayList<Term>();
	}

	@Override
	public void newOboFileHeader() {
		// TODO Auto-generated method stub
	}

	@Override
	public void newStanza(String stanza) {
		//logger.info("got a new stanza: {}", stanza);
		if ( stanza.equals(TERM)){
			isTerm = true;
			currentTerm = null;
		} else {
			isTerm = false;
		}

	}

	@Override
	public void newKey(String key, String value) {
		if (isTerm) {

			if ( key.equals(ID_KEY)) {
				if ( ontology.containsTerm(key)){
					currentTerm = ontology.getTerm(key);
				} else {
					try {
						if (  ontology.containsTerm(value)) {
							currentTerm = ontology.getTerm(value);
						} else {
							currentTerm = ontology.createTerm(value);
						}
					} catch (AlreadyExistsException ex) {
						logger.error("Exception: ", ex);
					}

				}
				return;
			}
			if (currentTerm == null) {
				logger.warn("did not find ID for Term! ");
				return;
			}
			if (key.equals(NAMESPACE)){
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(NAMESPACE, value);
			}
			else if (key.equals(NAME)){
				currentTerm.setDescription(value);
			} else if (key.equals(DEF)){
				//TODO
				// set definition
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(DEF, value);
			} else if (key.equals(XREF_ANALOG)){
				// set xref analog
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(XREF_ANALOG, value);
			} else if (key.equals(IS_OBSOLETE)) {
				// ignore obsolete Terms...
				//logger.info("obsolete: {}", currentTerm);
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(IS_OBSOLETE, new Boolean(true));

			} else if (key.equals(IS_A) ||
					key.equals(RELATIONSHIP) ||
					key.equals(DISJOINT_FROM) ||
					key.equals(INTERSECTION_OF) ||
					key.equals(SUBSET)) {
				try {
					Term object = (ontology.containsTerm(value) ?
							ontology.getTerm(value): ontology.createTerm(value));
					Term predicate = (ontology.containsTerm(key) ? ontology.getTerm(key) : ontology.createTerm(key));
					ontology.createTriple(currentTerm, object, predicate, currentTerm + " " + predicate + " " + object, key+"-relationship");
				} catch (AlreadyExistsException ex) {
				}

			} else if (key.equals(COMMENT)){
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(COMMENT, value);
			} else if (key.equals(ALT_ID)){
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(ALT_ID, value);
			}
			else if (key.equals(REPLACED_BY)) {
				Annotation anno = currentTerm.getAnnotation();
				anno.setProperty(REPLACED_BY, value);
			}

			else {
				//logger.info("unknown key {}", key);
			}


		} else {
			//logger.info("not a term and ignoring: {}->{}", key, value);
		}

	}

	@Override
	public void newSynonym(Synonym synonym) {
		if (isTerm) {
			currentTerm.addSynonym(synonym);
		}
	}

}