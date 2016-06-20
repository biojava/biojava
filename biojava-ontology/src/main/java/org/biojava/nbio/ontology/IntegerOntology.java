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
import org.biojava.nbio.ontology.utils.WeakValueHashMap;

import java.util.*;



/**
 *
 *
 * @author Matthew Pocock
 */
public class IntegerOntology

implements Ontology {
	private final Map termCache;

	IntegerOntology() {
		termCache = new WeakValueHashMap();
	}

	@Override
	public String getName() {
		return "core.integer";
	}

	@Override
	public String getDescription() {
		return "Ontology containing all integers";
	}

	@Override
	public void setDescription(String description){
	}

	@Override
	public Set getTerms() {
		return new AbstractSet() {
			@Override
			public boolean contains(Object o) {
				return o instanceof IntTerm;
			}

			@Override
			public int size() {
				return Integer.MAX_VALUE;
			}

			@Override
			public Iterator iterator() {
				return new Iterator() {
					int i = 0;

					@Override
					public boolean hasNext() {
						return i > 0;
					}

					@Override
					public Object next() {
                        if(!hasNext()){
                            throw new NoSuchElementException();
                        }
                        return resolveInt(i++);
					}

					@Override
					public void remove() {
						throw new UnsupportedOperationException();
					}
				};
			}
		};
	}

	@Override
	public Term getTerm(String s) throws NoSuchElementException {
		int val = Integer.parseInt(s);
		return resolveInt(val);
	}

	@Override
	public Set getTriples(Term subject, Term object, Term predicate) {
		return Collections.EMPTY_SET;
	}

	@Override
	public OntologyOps getOps() {
		return new DefaultOps() {
			@Override
			public Set getRemoteTerms() {
				return Collections.EMPTY_SET;
			}
		};
	}

	@Override
	public Term createTerm(String name) throws AlreadyExistsException,  IllegalArgumentException {
		throw new IllegalArgumentException(getName() + " is immutable");
	}

	@Override
	public Term createTerm(String name, String description)
			throws
			AlreadyExistsException,

			IllegalArgumentException
			{
		throw new IllegalArgumentException(getName() + " is immutable");
			}

	@Override
	public Term createTerm(String name, String description, Object[] synonyms)
			throws
			AlreadyExistsException,

			IllegalArgumentException
			{
		throw new IllegalArgumentException(getName() + " is immutable");
			}

	@Override
	public Variable createVariable(String name, String description)
			throws
			AlreadyExistsException,

			IllegalArgumentException
			{
		throw new IllegalArgumentException(getName() + " is immutable");
			}

	@Override
	public Term importTerm(Term t, String name)


	{
		throw new IllegalArgumentException(getName() + " is immutable");
	}

	@Override
	public Triple createTriple(Term subject, Term object, Term predicate, String name, String description)
			throws
			AlreadyExistsException
			{
		throw new IllegalArgumentException(getName() + " is immutable");
			}

	@Override
	public boolean containsTriple(Term subject, Term object, Term predicate) {
		return false;
	}

	@Override
	public void deleteTerm(Term t)  {
		throw new RuntimeException(getName() + " is immutable");
	}

	@Override
	public boolean containsTerm(String name) {
		// uglee hack - perhaps we should use a regex?
		try {
			Integer.parseInt(name);
		} catch (NumberFormatException e) {
			return false;
		}

		return true;
	}

	public IntTerm resolveInt(int val) {
		Integer i = new Integer(val);
		IntTerm term = (IntTerm) termCache.get(i);

		if(term == null) {
			term = new IntTerm(val);
			termCache.put(i, term);
		}

		return term;
	}

	public final class IntTerm

	implements Term {
		private final int val;
		private Set synonyms;

		public IntTerm(int val) {
			this(val, null);
		}

		public IntTerm(int val, Object[] synonyms) {
			this.val = val;

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

		public int intValue() {
			return val;
		}

		@Override
		public String getName() {
			return String.valueOf(val);
		}

		@Override
		public String getDescription() {
			return "The integer " + getName();
		}

		@Override
		public void setDescription(String description){

		}

		@Override
		public Ontology getOntology() {
			return IntegerOntology.this;
		}

		@Override
		public Annotation getAnnotation() {
			return Annotation.EMPTY_ANNOTATION;
		}
	}

	@Override
	public void setName(String name) {
		//ignore

	}


}
