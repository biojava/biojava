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
 * A triple in an ontology.  This is two terms and a relationship between
 * them, similar to RDF and other similar logic systems.
 *
 * <p>
 * For documentation purposes, a Triple may provide a name. However, a Triple
 * may also be named as "(subject, object, predicate)" if no specific name is
 * provided.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.4
 * @see org.biojavax.ontology.ComparableTriple
 */

public interface Triple
extends Term {
	/**
	 * Return the subject term of this triple
	 * @return the subject term
	 */

	public Term getSubject();

	/**
	 * Return the object term of this triple.
	 * @return the object term
	 */

	public Term getObject();

	/**
	 * Return a Term which defines the type of relationship between the subject and object terms.
	 * @return the predicate
	 */

	public Term getPredicate();

	/**
	 * The hashcode for a Triple.
	 *
	 * <p>This <em>must</em> be implemented as:
	 * <pre>
	 * return getSubject().hashCode() +
	 * 31 * getObject().hashCode() +
	 * 31 * 31 * getPredicate().hashCode();
	 * </pre>
	 * If you do not implement hashcode in this way then you have no guarantee
	 * that your Triple objects will be found in an ontology and that they will
	 * not be duplicated.
	 * </p>
	 */
	@Override
	public int hashCode();

	/**
	 * Check to see if an object is an equivalent Triple.
	 *
	 * <p>
	 * Two triples are equivalent if they have the same subject, object and
	 * predicate fields.
	 * <pre>
	 * if (! (o instanceof Triple)) {
	 *     return false;
	 * }
	 * Triple to = (Triple) o;
	 * return to.getSubject() == getSubject() &&
	 *        to.getObject() == getObject() &&
	 *        to.getPredicate() == getPredicate();
	 * </pre>
	 * If you do not implement equals in this way then you have no guarantee
	 * that your Triple objects will be found in an ontology and that they will
	 * not be duplicated.
	 * </p>
	 */
	@Override
	public boolean equals(Object obj);

	/**
	 * Basic in-memory implementation of a Triple in an ontology
	 *
	 * This can be used to implement Ontology.createTriple
	 * @see org.biojavax.ontology.SimpleComparableTriple
	 */

	public static final class Impl

	implements Triple, java.io.Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 3807331980372839221L;
		private final Term subject;
		private final Term object;
		private final Term predicate;
		private /*final*/ String name;
		private /*final*/ String description;
		private Set<Object> synonyms;

		public Impl(Term subject, Term object, Term predicate) {
			this(subject, object, predicate, null, null, null);
		}

		public Impl(Term subject, Term object, Term predicate, Object[] synonyms) {
			this(subject, object, predicate, null, null, synonyms);
		}

		public Impl(Term subject,
				Term object,
				Term predicate,
				String name,
				String description) {
			this(subject,object,predicate,name,description,null);
		}

		public Impl(Term subject,
				Term object,
				Term predicate,
				String name,
				String description,
				Object[] synonyms)
		{
			if (subject == null) {
				throw new NullPointerException("Subject must not be null");
			}
			if (object == null) {
				throw new NullPointerException("Object must not be null");
			}
			if (predicate == null) {
				throw new NullPointerException("predicate must not be null");
			}

			if(
					subject.getOntology() != object.getOntology() ||
					subject.getOntology() != predicate.getOntology()
			) {
				throw new IllegalArgumentException(
						"All terms must be from the same ontology: " +
						subject.getOntology().getName() + ", " +
						object.getOntology().getName() + ", " +
						predicate.getOntology().getName());
			}

			if(description == null) {
				description = "";
			}

			this.subject = subject;
			this.object = object;
			this.predicate = predicate;
			this.name = name;
			this.description = description;

			this.synonyms = new TreeSet<Object>();
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
			if(name == null) {
				name = predicate + "(" + subject + ", " + object + ")";
			}
			return name;
		}

		@Override
		public String getDescription() {
			return description;
		}
		@Override
		public void setDescription(String desc){
			this.description = desc;
		}

		@Override
		public Ontology getOntology() {
			return subject.getOntology();
		}

		@Override
		public Term getSubject() {
			return subject;
		}

		@Override
		public Term getObject() {
			return object;
		}

		@Override
		public Term getPredicate() {
			return predicate;
		}

		@Override
		public Annotation getAnnotation() {
			return Annotation.EMPTY_ANNOTATION;
		}

		/**
		 * Two triples are equal if all their fields are identical.
		 */

		@Override
		public boolean equals(Object o) {
			if (! (o instanceof Triple)) {
				return false;
			}
			Triple to = (Triple) o;
			return to.getSubject().equals(getSubject()) &&
			to.getObject().equals(getObject()) &&
			to.getPredicate().equals(getPredicate());
		}

		@Override
		public int hashCode() {
			return getSubject().hashCode() +
			31 * getObject().hashCode() +
			31 * 31 * getPredicate().hashCode();
		}

		@Override
		public String toString() {
			if (getName().length() > 0)
				return getName();
			return subject + " " + predicate + " " + object;
		}
	}
}
