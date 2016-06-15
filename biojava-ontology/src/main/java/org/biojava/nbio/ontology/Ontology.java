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

import org.biojava.nbio.ontology.utils.AssertionFailure;

import java.util.*;




/**
 * An ontology.
 *
 * <p>This is just a set of Term objects, and a set of
 * Triple objects describing relationships between these terms.
 * This class does not itself contain any reasoning functionality. Ontology is
 * a collection of facts, or axioms.</p>
 *
 * @author Thomas Down
 * @author Matthew Pocock
 *
 * @since 1.4
 * @see org.biojavax.ontology.ComparableOntology
 */

public interface Ontology  {


	/**
	 * Return the name of this ontology
	 * @return the name of the ontology
	 */

	public String getName();

	/** Set the name for this ontology
	 *
	 * @param name - the name
	 *
	 */
	public void setName(String name);

	/**
	 * Return a human-readable description of this ontology, or the empty
	 * string if none is available
	 * @return the description of the term
	 */

	public String getDescription();

	/** set the description of this ontology
	 *
	 * @param description
	 */
	public void setDescription(String description);



	/**
	 * Return all the terms in this ontology
	 * @return a Set of all Terms of the ontology.
	 */

	public Set<Term> getTerms();

	/**
	 * Fetch the term with the specified name.
	 * @param name the name of the term
	 *
	 * @return The term named <code>name</code>
	 * @throws NoSuchElementException if no term exists with that name
	 */

	public Term getTerm(String name) throws NoSuchElementException;

	/**
	 * Return all triples from this ontology which match the supplied
	 * pattern.  If any of the parameters of this method are <code>null</code>,
	 * they are treated as wildcards.
	 *
	 * @param subject The subject to search for, or <code>null</code>
	 * @param object The object to search for, or <code>null</code>
	 * @param predicate The relationship to search for, or <code>null</code>.
	 * @return a Set of triples
	 */

	public Set<Triple> getTriples(Term subject, Term object, Term predicate);

	/**
	 * Return the associated OntologyOps.
	 *
	 * This method should be implemented by ontology
	 * implementors to allow OntoTools
	 * to get optimized access to some usefull ontology operations. It is not
	 * intended that users will ever invoke this. A sensible dumb implementation
	 * of this would return a per-ontology instance of DefaultOps.
	 *
	 * @return the OntologyOps instance associated with this instance.
	 */

	public OntologyOps getOps();

	/**
	 * Create a new term in this ontology.
	 *
	 * @param name The name of the term (must be unique))
	 * @throws IllegalArgumentException if either <code>name</code> or
	 *         <code>description</code> is <code>null</code>, or violates
	 *         some other constraint of this implementation.
	 * @throws AlreadyExistsException if a term of this name already exists
	 * @return The newly created term.
	 * @throws ChangeVetoException
	 */

	public Term createTerm(String name)
			throws
			AlreadyExistsException,

			IllegalArgumentException;

	/**
	 * Create a new term in this ontology.
	 *
	 * @param name The name of the term (must be unique)
	 * @param description A human-readable description (may be empty)
	 * @throws IllegalArgumentException if either <code>name</code> or
	 *         <code>description</code> is <code>null</code>, or violates
	 *         some other constraint of this implementation.
	 * @throws AlreadyExistsException if a term of this name already exists
	 * @return The newly created term.
	 * @throws ChangeVetoException
	 */

	public Term createTerm(String name, String description)
			throws
			AlreadyExistsException,

			IllegalArgumentException;

	/**
	 * Create a new term in this ontology.
	 *
	 * @param name The name of the term (must be unique)
	 * @param description A human-readable description (may be empty)
	 * @param synonyms Some synonyms for this term.
	 * @throws IllegalArgumentException if either <code>name</code> or
	 *         <code>description</code> is <code>null</code>, or violates
	 *         some other constraint of this implementation.
	 * @throws AlreadyExistsException if a term of this name already exists
	 * @return The newly created term.
	 * @throws ChangeVetoException
	 */

	public Term createTerm(String name, String description, Object[] synonyms)
			throws
			AlreadyExistsException,

			IllegalArgumentException;

	/**
	 * Create a new term in this ontology that is used as a variable.
	 *
	 * @param name The name of the term (must be unique)
	 * @param description A human-readable description (may be empty)
	 * @throws IllegalArgumentException if either <code>name</code> or
	 *         <code>description</code> is <code>null</code>, or violates
	 *         some other constraint of this implementation.
	 * @throws AlreadyExistsException if a term of this name already exists
	 * @return The newly created term.
	 * @throws ChangeVetoException
	 */

	public Variable createVariable(String name, String description)
			throws
			AlreadyExistsException,

			IllegalArgumentException;

	/**
	 * Create a view of a term from another ontology.  If the requested term
	 * has already been imported under that name, this method returns the existing
	 * RemoteTerm object. If the term that is being imported is itself a
	 * RemoteTerm instance then first unwrap the term back to the orriginal
	 * term it represents and then produce a RemoteTerm from that. If the term
	 * being imported orriginated from this ontology, then return that term
	 * unaltered.
	 *
	 * @param t  the Term to import
	 * @param localName  the local name to import it under, optionally null
	 * @return a Term
	 * @throws ChangeVetoException
	 * @throws IllegalArgumentException
	 */

	public Term importTerm(Term t, String localName)
			throws

			IllegalArgumentException;

	/**
	 * Creates a new Triple.
	 *
	 * @param subject     the subject Term
	 * @param object      the object Term
	 * @param predicate    the predicate Term
	 * @param name        the name of the triple, or null
	 * @param description the description of the triple, or null
	 * @return  a new Triple over these three terms
	 * @throws AlreadyExistsException if a triple already exists with the same
	 *      subject, object and predicate, regardless of the name and description
	 * @throws ChangeVetoException
	 * @throws NullPointerException if subject, object or predicate are null
	 * @throws IllegalArgumentException if subject, object or predicate are not all
	 *      from the same ontology
	 */
	public Triple createTriple(Term subject, Term object, Term predicate, String name, String description)
			throws
			AlreadyExistsException
			;

	/**
	 * See if a triple exists in this ontology
	 * @param subject
	 * @param object
	 * @param predicate
	 * @return true if contained
	 */

	public boolean containsTriple(Term subject, Term object, Term predicate);

	/**
	 * Remove a term from an ontology, together with all triples which refer to it.
	 * @param t
	 * @throws ChangeVetoException
	 */

	public void deleteTerm(Term t) ;

	/**
	 * Determines if this ontology currently contains a term named <code>name</code>
	 * @param name
	 * @return true is contained
	 */

	public boolean containsTerm(String name);

	/**
	 * A basic in-memory implementation of an ontology
	 *
	 * @author Thomas Down
	 * @author Matthew Pocock
	 * @since 1.3
	 */

	// AP: I am setting name and description to public changeable fields
	// e.g during parsing of an .obo file we don't know them when the ontology is instanciated
	public final class Impl

	implements Ontology, java.io.Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = -8064461497813727957L;
		private final Map<String,Term> terms;
		private final Set<Triple> triples;
		private final Map<Term, Set<Triple>> subjectTriples;
		private final Map<Term, Set<Triple>> objectTriples;
		private final Map<Term, Set<Triple>> relationTriples;
		private final Map<Term,RemoteTerm> remoteTerms;
		private final Set<Term> localRemoteTerms;

		private /*final*/ String name;
		private /*final*/ String description;
		private final OntologyOps ops;

		{
			terms            = new HashMap<String, Term>();
			triples          = new HashSet<Triple>();
			subjectTriples   = new HashMap<Term, Set<Triple>>();
			objectTriples    = new HashMap<Term, Set<Triple>>();
			relationTriples  = new HashMap<Term, Set<Triple>>();
			remoteTerms      = new HashMap<Term, RemoteTerm>();
			localRemoteTerms = new HashSet<Term>();
		}

		public Impl(String name, String description) {
			this.name = name;
			this.description = description;
			ops = new DefaultOps() {
				/**
				 *
				 */
				private static final long serialVersionUID = -2135777733685713181L;

				@Override
				public Set<Term> getRemoteTerms() {
					return localRemoteTerms;
				}
			};
		}

		@Override
		public String getName() {
			return name;
		}

		@Override
		public String getDescription() {
			return description;
		}


		@Override
		public void setDescription(String description){
			this.description = description;
		}

		@Override
		public Set<Term> getTerms() {
			return new HashSet<Term>(terms.values());
		}

		@Override
		public Term getTerm(String name)
				throws NoSuchElementException
				{
			Term t = terms.get(name);
			if (t == null) {
				throw new NoSuchElementException("No term named '" + name + "'");
			} else {
				return terms.get(name);
			}
				}

		@Override
		public Set<Triple> getTriples(Term subject, Term object, Term predicate) {
			if(subject != null && subject.getOntology() != this) {
				throw new IllegalArgumentException("Subject is not in this ontology: " + subject + " " + this);
			}

			if(object != null && object.getOntology() != this) {
				throw new IllegalArgumentException("Object is not in this ontology: " + object + " " + this);
			}

			if(predicate != null && predicate.getOntology() != this) {
				throw new IllegalArgumentException("Predicate is not in this ontology: " + predicate + " " + this);
			}

			if (subject != null) {
				return filterTriples( subjectTriples.get(subject), null, object, predicate);
			} else if (object != null) {
				return filterTriples( objectTriples.get(object), subject, null, predicate);
			} else if (predicate != null) {
				return filterTriples( relationTriples.get(predicate), subject, object, null);
			} else {
				return filterTriples(triples, subject, object, predicate);
			}
		}

		private Set<Triple> filterTriples(Set<Triple> base, Term subject, Term object, Term predicate) {
			if (base == null) {
				return Collections.EMPTY_SET;
			} else if (subject == null && object == null && predicate == null) {
				return Collections.unmodifiableSet(new HashSet<Triple>(base));
			}

			Set<Triple> retval = new HashSet<Triple>();
			for (Iterator<Triple> i = base.iterator(); i.hasNext(); ) {
				Triple t = i.next();
				if (subject != null && t.getSubject() != subject) {
					continue;
				}
				if (object != null && t.getObject() != object) {
					continue;
				}
				if (predicate != null && t.getPredicate() != predicate) {
					continue;
				}
				retval.add(t);
			}
			return retval;
		}

		private void addTerm(Term t)
				throws AlreadyExistsException, IllegalArgumentException
				{
			if (terms.containsKey(t.getName())) {
				throw new AlreadyExistsException("Ontology " + getName() + " already contains " + t.toString());
			}


			terms.put(t.getName(), t);

				}

		@Override
		public Term createTerm(String name)
				throws AlreadyExistsException, IllegalArgumentException
				{
			Term t = new Term.Impl(this, name);
			addTerm(t);
			return t;
				}

		@Override
		public Term createTerm(String name, String description)
				throws AlreadyExistsException, IllegalArgumentException
				{
			Term t = new Term.Impl(this, name, description);
			addTerm(t);
			return t;
				}

		@Override
		public Term createTerm(String name, String description, Object[] synonyms)
				throws AlreadyExistsException, IllegalArgumentException
				{
			Term t = new Term.Impl(this, name, description, synonyms);
			addTerm(t);
			return t;
				}

		@Override
		public Variable createVariable(String name, String description)
				throws
				AlreadyExistsException,

				IllegalArgumentException {
			Variable var = new Variable.Impl(this, name, description);
			addTerm(var);
			return var;
		}

		public OntologyTerm createOntologyTerm(Ontology o)
				throws AlreadyExistsException
				{
			OntologyTerm ot = new OntologyTerm.Impl(this, o);
			addTerm(ot);
			return ot;
				}


		@Override
		public Term importTerm(Term t, String name)
				throws IllegalArgumentException
				{
			// unpack any potential indirection - belt & braces
			while(t instanceof RemoteTerm) {
				t = ((RemoteTerm) t).getRemoteTerm();
			}

			if(t.getOntology() == this) {
				return t;
			}

			RemoteTerm rt = remoteTerms.get(t);
			if (rt == null) {
				rt = new RemoteTerm.Impl(this, t, name);
				try {
					addTerm(rt);
				} catch (AlreadyExistsException e) {
					throw new AssertionFailure("This term can not exist", e);
				}
				if(name == null) {
					remoteTerms.put(t, rt);
				}
				localRemoteTerms.add(rt);
			}
			return rt;
				}

		@Override
		public void deleteTerm(Term t)

		{
			String name = t.getName();
			if (terms.get(name) != t) {
				return; // Should this be an exception?
			}

			terms.remove(name);
			if(t instanceof Triple) {
				removeTriple((Triple) t);
			}

		}

		@Override
		public boolean containsTerm(String name) {
			return terms.containsKey(name);
		}

		private boolean containsTerm(Term t) {
			return (terms.get(t.getName()) == t);
		}

		@Override
		public boolean containsTriple(Term subject, Term object, Term predicate) {
			if(!(subject.getOntology() == this)) return false;
			if(!(object.getOntology() == this)) return false;
			if(!(predicate.getOntology() == this)) return false;

			return triples.contains(new Triple.Impl(subject, object, predicate));
		}

		@Override
		public Triple createTriple(Term subject,
				Term object,
				Term predicate,
				String name,
				String description)
						throws
						AlreadyExistsException,
						IllegalArgumentException,
						NullPointerException,
						IllegalArgumentException
						{
			Triple t = new Triple.Impl(subject, object, predicate, name, description);
			if (!containsTerm(subject)) {
				throw new IllegalArgumentException("Ontology " + getName() + " doesn't contain " + subject);
			}
			if (!containsTerm(predicate)) {
				throw new IllegalArgumentException("Ontology " + getName() + " doesn't contain " + predicate);
			}
			if (!containsTerm(object)) {
				throw new IllegalArgumentException("Ontology " + getName() + " doesn't contain " + object);
			}
			if (triples.contains(t)) {
				throw new AlreadyExistsException("Ontology " + getName() + " already contains " + t.toString());
			}


			addTerm(t);
			addTriple(t);

			return t;
						}

		private void addTriple(Triple t) {
			triples.add(t);
			pushTriple(subjectTriples, t.getSubject(), t);
			pushTriple(objectTriples, t.getObject(), t);
			pushTriple(relationTriples, t.getPredicate(), t);
		}

		private void pushTriple(Map<Term,Set<Triple>> m, Term key, Triple t) {
			Set<Triple> s = m.get(key);
			if (s == null) {
				s = new HashSet<Triple>();
				m.put(key, s);
			}
			s.add(t);
		}

		private void removeTriple(Triple t) {
			triples.remove(t);
			pullTriple(subjectTriples, t.getSubject(), t);
			pullTriple(objectTriples, t.getObject(), t);
			pullTriple(relationTriples, t.getPredicate(), t);
		}

		private void pullTriple(Map<Term, Set<Triple>> subjectTriples2, Term key, Triple t) {
			Set<Triple> s = subjectTriples2.get(key);
			if (s != null) {
				s.remove(t);
			}
		}

		@Override
		public OntologyOps getOps() {
			return ops;
		}

		@Override
		public String toString() {
			return "ontology: " + getName();
		}

		@Override
		public void setName(String name) {
			this.name=name;

		}
	}
}

