package org.biojava3.ontology;

import java.util.AbstractSet;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.TreeSet;

import org.biojava3.ontology.utils.Annotation;
import org.biojava3.ontology.utils.ChangeVetoException;
import org.biojava3.ontology.utils.Unchangeable;
import org.biojava3.ontology.utils.WeakValueHashMap;



/**
 *
 *
 * @author Matthew Pocock
 */
public class IntegerOntology
extends Unchangeable
implements Ontology {
	private final Map termCache;

	IntegerOntology() {
		termCache = new WeakValueHashMap();
	}

	public String getName() {
		return "core.integer";
	}

	public String getDescription() {
		return "Ontology containing all integers";
	}
	
	public void setDescription(String description){		
	}

	public Set getTerms() {
		return new AbstractSet() {
			public boolean contains(Object o) {
				return o instanceof IntTerm;
			}

			public int size() {
				return Integer.MAX_VALUE;
			}

			public Iterator iterator() {
				return new Iterator() {
					int i = 0;

					public boolean hasNext() {
						return i > 0;
					}

					public Object next() {
						return resolveInt(i++);
					}

					public void remove() {
						throw new UnsupportedOperationException();
					}
				};
			}
		};
	}

	public Term getTerm(String s) throws NoSuchElementException {
		int val = Integer.parseInt(s);
		return resolveInt(val);
	}

	public Set getTriples(Term subject, Term object, Term predicate) {
		return Collections.EMPTY_SET;
	}

	public OntologyOps getOps() {
		return new DefaultOps() {
			public Set getRemoteTerms() {
				return Collections.EMPTY_SET;
			}
		};
	}

	public Term createTerm(String name) throws AlreadyExistsException, ChangeVetoException, IllegalArgumentException {
		throw new ChangeVetoException(getName() + " is immutable");
	}

	public Term createTerm(String name, String description)
	throws
	AlreadyExistsException,
	ChangeVetoException,
	IllegalArgumentException
	{
		throw new ChangeVetoException(getName() + " is immutable");
	}

	public Term createTerm(String name, String description, Object[] synonyms)
	throws
	AlreadyExistsException,
	ChangeVetoException,
	IllegalArgumentException
	{
		throw new ChangeVetoException(getName() + " is immutable");
	}

	public Variable createVariable(String name, String description)
	throws
	AlreadyExistsException,
	ChangeVetoException,
	IllegalArgumentException
	{
		throw new ChangeVetoException(getName() + " is immutable");
	}

	public Term importTerm(Term t, String name)
	throws
	ChangeVetoException
	{
		throw new ChangeVetoException(getName() + " is immutable");
	}

	public Triple createTriple(Term subject, Term object, Term predicate, String name, String description)
	throws
	AlreadyExistsException,
	ChangeVetoException {
		throw new ChangeVetoException(getName() + " is immutable");
	}

	public boolean containsTriple(Term subject, Term object, Term predicate) {
		return false;
	}

	public void deleteTerm(Term t) throws ChangeVetoException {
		throw new ChangeVetoException(getName() + " is immutable");
	}

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
	extends Unchangeable
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

		public void addSynonym(Object synonym) {
			this.synonyms.add(synonym);
		}

		public void removeSynonym(Object synonym) {
			this.synonyms.remove(synonym);
		}

		public Object[] getSynonyms() {
			return this.synonyms.toArray();
		}

		public int intValue() {
			return val;
		}

		public String getName() {
			return String.valueOf(val);
		}

		public String getDescription() {
			return "The integer " + getName();
		}

		public void setDescription(String description){

		}

		public Ontology getOntology() {
			return IntegerOntology.this;
		}

		public Annotation getAnnotation() {
			return Annotation.EMPTY_ANNOTATION;
		}
	}

	public void setName(String name) {
		//ignore
		
	}


}
