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

import org.biojava.nbio.ontology.utils.Annotation;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;



/**
 * A term in another ontology.
 *
 * <p>
 * This is how you allow one ontology to refer to terms in another one. Since
 * these ontologies are designed to be modular and self-contained, it is
 * expected that you would not copy terms from one ontology into another. The
 * best-practice way to represent terms from another ontology in your one is to
 * use RemoteTerm instances. Ontology has a method importTerm that does this
 * for you. By default, imported terms will have names composed from the source
 * ontology and the imported term name. However, this should be over-rideable.
 * </p>
 *
 * <p>
 * The imported term will have the same name as the original term. They are
 * implicitly identical to each other. The most common use of imports will be
 * to slurp in the "core" ontology so that operations such as <code>is-a</code>
 * and <code>has-a</code> are available.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.4
 */

public interface RemoteTerm extends Term {
	/**
	 * Return the imported term
	 * @return the term
	 */

	public Term getRemoteTerm();

	/**
	 * Simple in-memory implementation of a remote ontology term.
	 *
	 * This can be used to implement Ontology.importTerm
	 */

	public final static class Impl
	extends AbstractTerm
	implements RemoteTerm, java.io.Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 922700041939183676L;
		private final Ontology ontology;
		private final Term remoteTerm;
		private final String name;
		private Set synonyms;

		public Impl(Ontology ontology, Term remoteTerm, String name) {
			this(ontology, remoteTerm, name, null);
		}

		public Impl(Ontology ontology, Term remoteTerm, String name, Object[] synonyms) {
			if (ontology == null) {
				throw new NullPointerException("Ontology must not be null");
			}
			if (remoteTerm == null) {
				throw new NullPointerException("RemoteTerm must not be null");
			}
			if(name == null) {
			  name = remoteTerm.getOntology().getName() + "." + remoteTerm.getName();
			}

			this.ontology = ontology;
			this.remoteTerm = remoteTerm;
			this.name = name;

			this.synonyms = new TreeSet();
			if (synonyms!=null) this.synonyms.addAll(Arrays.asList(synonyms));
		}

		@Override
		public void addSynonym(Object synonym) {
			this.synonyms.add(synonym);
		}

		@Override
		public void removeSynonym(Object synonym) {
			this.synonyms.remove(synonym);
		}

		@Override
		public Object[] getSynonyms() {
			return this.synonyms.toArray();
		}

		@Override
		public String getName() {
			return getOntology().getName() + ":" + remoteTerm.getName();
		}

		@Override
		public String getDescription() {
			return remoteTerm.getDescription();
		}

		@Override
		public Ontology getOntology() {
			return ontology;
		}

		@Override
		public Term getRemoteTerm() {
			return remoteTerm;
		}

		@Override
		public String toString() {
			return name;
		}

		@Override
		public Annotation getAnnotation() {
			return remoteTerm.getAnnotation();
		}
	}
}
