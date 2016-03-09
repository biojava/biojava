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

import java.util.Set;

/**
 * This is an interface for optimizing ontology operators.
 *
 * <p>
 * Some ontology implementations will be able to compute some derived properties
 * very quickly because of how they store their data. This is likely to out-
 * perform generic implementations of algorithms using the Ontology interface
 * to get the same result. Ontology instances provide an instance of
 * OntologyOps, publishing optimizations of some common operations. The reasoner
 * may then choose to call OntologyOps methods on the Ontology instance rather
 * than using its fall-back implementations.
 * </p>
 *
 * @author Matthew Pocock
 * @since 1.4
 */
public interface OntologyOps {
	/**
	 * Get the set of all remote terms.
	 *
	 * <p>
	 * We do not currently specify whether this set is mutable or not, and if it
	 * will reflect modifications to the optimised ontolgies.
	 * </p>
	 *
	 * @return a Set containing all remote terms in the ontology
	 */
	public Set getRemoteTerms();
}
