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
 * Created on 09.05.2004
 * @author Andreas Prlic
 *
 */

package org.biojava.nbio.structure;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;


/**
 * An iterator over all groups of a structure.
 * @author Andreas Prlic
 * @since 1.4
 * @version %I% %G%
 */

public class GroupIterator implements Iterator<Group> {

	private Structure structure   ;
	private int current_model_pos ;
	private int current_chain_pos ;
	private int current_group_pos ;
	private boolean fixed_model   ;


	/**
	 * Constructs a GroupIterator object over all models
	 *
	 * @param struct  a Structure object
	 */
	public GroupIterator (Structure struct) {
		structure = struct     ;
		current_model_pos = 0  ;
		current_chain_pos = 0  ;
		current_group_pos = -1 ;
		fixed_model = false    ;
	}
	/**
	 * Constructs a GroupIterator object over a specific model
	 *
	 * @param struct  a Structure object
	 */
	public GroupIterator (Structure struct, int modelNr) {
		structure = struct     ;
		current_model_pos = modelNr;
		current_chain_pos = 0  ;
		current_group_pos = -1 ;
		fixed_model = true     ;
	}


	/** needed to do a copy of iterator ... */
	private Structure getStructure() { return structure         ;}
	private int  getModelPos()       { return current_model_pos ;}
	private void setModelPos(int pos){ current_model_pos = pos  ;}
	private int  getChainPos()       { return current_chain_pos ;}
	private void setChainPos(int pos){ current_chain_pos = pos  ;}
	private int  getGroupPos()       { return current_group_pos ;}
	private void setGroupPos(int pos){ current_group_pos = pos  ;}

	/**  Creates and returns a copy of this object. */
	@Override
	public Object clone () {

		GroupIterator gr = new GroupIterator(this.getStructure()) ;
		gr.setModelPos(this.getModelPos());
		gr.setChainPos(this.getChainPos());
		gr.setGroupPos(this.getGroupPos());
		gr.fixed_model = this.fixed_model;
		return gr ;

	}


	/** is there a group after the current one in the structure? */
	@Override
	public boolean hasNext() {
		return hasSubGroup(current_model_pos,current_chain_pos , current_group_pos +1) ;
	}

	/** recursive method to determine if there is a next group. Helper
	 * method for hasNext().
	 * @see #hasNext
	 */
	private boolean hasSubGroup(int tmp_model,int tmp_chain,int tmp_group) {

		if (tmp_model >= structure.nrModels()) {
			return false;
		}

		List<Chain> model = structure.getModel(tmp_model);

		if (tmp_chain >= model.size()) {
			if(fixed_model)
				return false;
			return hasSubGroup(tmp_model + 1, 0, 0);
		}

		Chain chain = model.get(tmp_chain);

		// start search at beginning of next chain.
		return tmp_group < chain.getAtomLength() || hasSubGroup(tmp_model, tmp_chain + 1, 0);

	}

	/** Get the model number of the current model.
	 *
	 * @return the number of the model
	 */
	public int getCurrentModel(){

		return current_model_pos;
	}

	/** Get the current Chain. Returns null if we are at the end of the iteration.
	 *
	 * @return the Chain of the current position
	 */
	public Chain getCurrentChain(){
		if ( current_model_pos >= structure.nrModels()){
			return null;
		}

		List<Chain> model =  structure.getModel(current_model_pos);

		if ( current_chain_pos >= model.size() ){
			return null;
		}

		return model.get(current_chain_pos);

	}

	/** get next Group.
	 * @return next Group
	 * @throws NoSuchElementException ...
	 */
	@Override
	public Group next()
	throws NoSuchElementException
	{

		return getNextGroup(current_model_pos,current_chain_pos,current_group_pos+1);
	}

	/** recursive method to retrieve the next group. Helper
	 * method for gext().
	 * @see #next
	 */
	private Group getNextGroup(int tmp_model,int tmp_chain,int tmp_group)
	throws NoSuchElementException
	{

		if ( tmp_model >= structure.nrModels()){
			throw new NoSuchElementException("arrived at end of structure!");
		}

		List<Chain> model = structure.getModel(tmp_model);

		if ( tmp_chain >= model.size() ){
			if(fixed_model)
				throw new NoSuchElementException("arrived at end of model!");
			return getNextGroup(tmp_model+1,0,0);
		}

		Chain     chain = model.get(tmp_chain);

		if (tmp_group  >= chain.getAtomLength()){
			// start search at beginning of next chain.
			return getNextGroup(tmp_model,tmp_chain+1,0);
		} else {
			current_model_pos = tmp_model;
			current_chain_pos = tmp_chain;
			current_group_pos = tmp_group;
			return chain.getAtomGroup(current_group_pos);
		}

	}

	/** does nothing .
	 */
	@Override
	public void remove() {
		throw new UnsupportedOperationException("Cannot call remove() for GroupIterator");
	}


}

