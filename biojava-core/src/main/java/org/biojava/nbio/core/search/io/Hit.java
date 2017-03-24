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
package org.biojava.nbio.core.search.io;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

import org.biojava.nbio.core.sequence.template.Compound;
import org.biojava.nbio.core.sequence.template.Sequence;

/**
 * This class models a search Hit.
 * You will retrieve a list of this using iterator of a Result
 * Designed by Paolo Pavan.
 * You may want to find my contacts on Github and LinkedIn for code info
 * or discuss major changes.
 * https://github.com/paolopavan
 *
 * @author Paolo Pavan
 */

public abstract class Hit<S extends Sequence<C>,C extends Compound> implements Iterable<Hsp<S,C>>{
	private final int hitNum;
	private final String hitId;
	private final String hitDef;
	private final String hitAccession;
	/**
	 * the length of the hit sequence
	 */
	private final int hitLen;
	private final List<Hsp<S,C>> hsps;
	private S hitSequence;



	public Hit(int hitNum, String hitId, String hitDef, String hitAccession, int hitLen, List<Hsp<S,C>> hsps, S hitSequence) {
		this.hitNum = hitNum;
		this.hitId = hitId;
		this.hitDef = hitDef;
		this.hitAccession = hitAccession;
		this.hitLen = hitLen;
		this.hsps = hsps;
		this.hitSequence = hitSequence;
	}

	@Override
	public int hashCode() {
		int hash = 3;
		hash = 89 * hash + this.hitLen;
		hash = 89 * hash + (this.hsps != null ? this.hsps.hashCode() : 0);
		return hash;
	}
	 /**
	 * Implements conceptual comparisons of search results.
	 * Fields unrelated to search are deliberately not considered.
	 * @return
	 */
	@Override
	public boolean equals(Object obj) {
		if (obj == null) {
			return false;
		}
		if (getClass() != obj.getClass()) {
			return false;
		}
		final Hit<?,?> other = (Hit<?,?>) obj;
		if (this.hitLen != other.hitLen) {
			return false;
		}
		if( this.hsps == null ) {
			return other.hsps == null;
		} else {
			return this.hsps.equals(other.hsps);
		}
	}

	public int getHitNum() {
		return hitNum;
	}

	public String getHitId() {
		return hitId;
	}

	public String getHitDef() {
		return hitDef;
	}

	public String getHitAccession() {
		return hitAccession;
	}

	public int getHitLen() {
		return hitLen;
	}

	/**
	 * returns the reference to the original and whole sequence hit in the database.
	 * Available only if the ResultFactory implements setHitReferences and
	 * it was used before the parsing with SearchIO
	 * @return Sequence object
	 */
	public S getHitSequence() {
		return hitSequence;
	}

	@Override
	public Iterator<Hsp<S,C>> iterator() {
		return new Iterator<Hsp<S,C>>() {
			int current = 0;
			@Override
			public boolean hasNext() {
				return current < hsps.size();
			}

			@Override
			public Hsp<S,C> next() {
                if(!hasNext()){
                    throw new NoSuchElementException();
                }
                return hsps.get(current++);
			}

			@Override
			public void remove() {
				throw new UnsupportedOperationException("The remove operation is not supported by this iterator");
			}
		};
	}
}
