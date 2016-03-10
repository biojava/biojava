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

import org.biojava.nbio.ontology.utils.Annotatable;
import org.biojava.nbio.ontology.utils.Annotation;
import org.biojava.nbio.ontology.utils.SmallAnnotation;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;



/**
 * A term in an ontology.  This has an {@link org.biojava.nbio.Annotation Annotation}
 * which can be used for storing additional human-displayable information.  It
 * is strongly recommended that the Annotation is not used for any machine-readable
 * data -- this should be represented by relations in the ontology instead.
 *
 * <p>
 * Terms are things that represent things. They are the same sort of thing as a
 * Java object or a prolog atom. A sub-set of terms are themselves relations.
 * This means that they are used to describe associations between pairs of terms.
 * Since all terms can be described, it is possible (and indeed encouraged) to
 * describe relations. As a minimum, you should consider saying if they are
 * identity or partial order relations, or if they are transitive, reflexive,
 * symmetrical, anti-symmetrical or anything else you know about them. This gives
 * the inference engine some chance of working out what is going on.
 * </p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 * @since 1.4
 * @see org.biojavax.ontology.ComparableTerm
 */

public interface Term extends Annotatable {
	/**
	 * ChangeType which indicates that this term's ontology has been
	 * altered
	 */


	/**
	 * Return the name of this term.
	 * @return the name of the term
	 */

	public String getName();

	/**
	 * Return a human-readable description of this term, or the empty string if
	 * none is available.
	 * @return the description of the term
	 */

	public String getDescription();

	/** set the description of the term;
	 *
	 * @param description
	 *
	 */
	public void setDescription(String description);

	/**
	 * Return the ontology in which this term exists.
	 * @return the ontology
	 */

	public Ontology getOntology();

	/**
	 * Return the synonyms for this term.
	 * @return the synonyms
	 */

	public Object[] getSynonyms();

	/**
	 * Add a synonym for this term.
	 * @param synonym the synonym
	 */

	public void addSynonym(Object synonym);

	/**
	 * Remove a synonym for this term.
	 * @param synonym
	 */

	public void removeSynonym(Object synonym);

	/**
	 * Simple in-memory implementation of an ontology term.
	 * @see org.biojavax.ontology.SimpleComparableTerm
	 * This can be used to implement Ontology.createTerm
	 */

	public static class Impl
	extends AbstractTerm
	implements Term, java.io.Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 6561668917514377417L;

		private final String name;

		private final Ontology ontology;
		private Annotation annotation;
		private Set<Object> synonyms;

		public Impl(Ontology ontology, String name) {
			this(ontology,name,null,null);
		}

		public Impl(Ontology ontology, String name, String description) {
			this(ontology,name,description,null);
		}

		public Impl(Ontology ontology, String name, String description, Object[] synonyms) {
			if (name == null) {
				throw new NullPointerException("Name must not be null");
			}
			// by AP - description can change from now on...
			//if (description == null) {
			//    throw new NullPointerException("Description must not be null");
			//}
			if (ontology == null) {
				throw new NullPointerException("Ontology must not be null");
			}

			this.name = name;
			this.description = description;
			this.ontology = ontology;

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
			return name;
		}

		public void setAnnotation(Annotation annotation) {
			this.annotation = annotation;
		}

		public void setSynonyms(Set<Object> synonyms) {
			this.synonyms = synonyms;
		}

		@Override
		public String getDescription() {
			return description;
		}

		@Override
		public Ontology getOntology() {
			return ontology;
		}

		@Override
		public String toString() {
			return name;
		}

		@Override
		public Annotation getAnnotation() {
			if (annotation == null) {
				annotation = new SmallAnnotation();
			}
			return annotation;
		}

	  @Override
	public int hashCode() {
		int value = 17;
		if(getName() != null)
		  value *= 31 * getName().hashCode();
		return 17 * value;
	  }

	  @Override
	public boolean equals(Object obj)
	  {
		if(obj == this) return true;
		if(!(obj instanceof Term)) return false;

		Term that = (Term) obj;

		return this.getOntology() == that.getOntology() &&
				this.getName() == that.getName();
	  }
	}
}
