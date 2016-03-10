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
 */
package org.biojava.nbio.ontology;

/**
 * A factory for Ontology instances.
 *
 * @author Matthew Pocock
 */

public interface OntologyFactory {
	/**
	 * Creates a new Ontology
	 *
	 * @param name  the name to give the ontology
	 * @param description the description for the ontology
	 * @return an Ontology
	 * @throws NullPointerException if either name or description are null
	 * @throws OntologyException if the ontology could not be created
	 */
	public Ontology createOntology(String name, String description)
	throws OntologyException;
}
