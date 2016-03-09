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
 *
 *
 * @author Matthew Pocock
 */
public interface Variable
extends Term {
	public static class Impl
	extends Term.Impl
	implements Variable
	{
		private static final long serialVersionUID = 1L;
		public Impl(Ontology ontology, String name, String description) {
			super(ontology, name, description);
		}
		public Impl(Ontology ontology, String name, String description, Object[] synonyms) {
			super(ontology, name, description, synonyms);
		}
	}
}
