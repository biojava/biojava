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
 * Created on 24.05.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.nbio.structure;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Iterator;
import java.util.NoSuchElementException;


/** an iterator over all atoms of a structure / group.
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */

public class AtomIterator implements Iterator<Atom> {

	private final static Logger logger = LoggerFactory.getLogger(AtomIterator.class);

	private Group     group         ;
	private int current_atom_pos    ;
	private GroupIterator groupiter ;

	/**
	 * Constructs an AtomIterator object over all models
	 *
	 * @param struct  a Structure object
	 */
	public AtomIterator(Structure struct) {
		current_atom_pos = -1 ;

		groupiter = new GroupIterator(struct) ;
		if ( groupiter.hasNext() ) {
			group = groupiter.next() ;
		}
		else
			group = null ;
	}
	
	/**
	 * Constructs an AtomIterator object over a single model
	 *
	 * @param struct  a Structure object
	 */
	public AtomIterator(Structure struct,int modelNr) {
		current_atom_pos = -1 ;

		groupiter = new GroupIterator(struct,modelNr) ;
		if ( groupiter.hasNext() ) {
			group = groupiter.next() ;
		}
		else
			group = null ;
	}

	/** Get the  chain that contains the current atom.
	 *
	 * @return a Chain object
	 */
	public Chain getCurrentChain(){
		return groupiter.getCurrentChain();
	}


	/** Get the model number of the model containing the current atom.
	 *
	 * @return the number of the model
	 */
	public int getCurrentModel(){
		return groupiter.getCurrentModel();
	}

	/**
	 * Constructs an AtomIterator object.
	 *
	 * @param g  a Group object
	 */
	public AtomIterator(Group g) {
		group = g ;
		current_atom_pos = -1 ;
		groupiter = null ;
	}

	/** Is there a next atom ?
	 * @return true if there is an atom after the current one
	 * */
	@Override
	public boolean hasNext() {

		// trying to iterate over an empty structure...

		if ( group == null)
			return false;

		// if there is another group ...
		if ( current_atom_pos < group.size()-1 ) {
			return true ;
		} else {
			// search through the next groups if they contain an atom
			if (groupiter != null) {
				GroupIterator tmp = (GroupIterator) groupiter.clone() ;
				while (tmp.hasNext()) {
					Group tmpg = tmp.next() ;

					if ( tmpg.size() > 0 ) {
						return true ;
					}

				}
			} else {
				// just an iterator over one group ...
				return false ;
			}
		}
		return false ;
	}

	/** Return next atom.
	 *
	 * @return the next Atom
	 * @throws NoSuchElementException if there is no atom after the current one
	 */
	@Override
	public Atom next()
	throws NoSuchElementException
	{
		current_atom_pos++ ;
		if ( current_atom_pos >= group.size() ) {
			if ( groupiter == null ) {
				throw new NoSuchElementException("no more atoms found in group!");

			}
			if ( groupiter.hasNext() ) {
				group = groupiter.next() ;
				current_atom_pos = -1 ;
				return next();
			} else {
				throw new NoSuchElementException("no more atoms found in structure!");
			}
		}

		Atom a ;


		a = group.getAtom(current_atom_pos);
		if ( a == null) {
			logger.error("current_atom_pos {} group {} size: {}", current_atom_pos, group, group.size());

			throw new NoSuchElementException("error wile trying to retrieve atom");
		}

		return a ;

	}

	/** does nothing. */
	@Override
	public void remove() {
	}

}

