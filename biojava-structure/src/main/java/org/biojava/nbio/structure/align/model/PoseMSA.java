package org.biojava.nbio.structure.align.model;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * An implementation of a Pose for a {@link MultipleAlignment}.
 * 
 * @author Aleix Lafita
 * 
 */
public class PoseMSA extends PoseAbstractImpl{

	private static final long serialVersionUID = 378522458943613367L;
	private MultipleAlignment parent;

	/**
	 * Constructor.
	 * @param multAln the parent MultipleAlignment of the PoseMSA instance. Cannot be null.
	 * @return PoseMSA a PoseMSA instance linked to its parent MultipleAlignment.
	 * @throws StructureAlignmentException if parent is null.
	 */
	public PoseMSA(MultipleAlignment multAln) throws StructureAlignmentException {
		super();
		
		if (multAln == null) throw new StructureAlignmentException("Parent of Pose is null.");
		parent = multAln;
	}
	
	/**
	 * Copy constructor.
	 * @param p PoseMSA object to be copied.
	 * @return PoseMSA an identical copy of the input PoseMSA object.
	 */
	public PoseMSA(PoseMSA p) {
		
		parent = p.parent;
		
		rotation = null;
		if (p.rotation!=null){
			//Make a deep copy of everything
			rotation = new ArrayList<Matrix>();
			for (int i=0; i<p.rotation.size(); i++)
				rotation.add((Matrix) p.rotation.get(i).clone());
		}
		
		translation = null;
		if (p.translation!=null){
			//Make a deep copy of everything
			translation = new ArrayList<Atom>();
			for (int i=0; i<p.translation.size(); i++)
				translation.add((Atom) p.translation.get(i).clone());
		}
		
		backDistMatrix = null;
		if (p.backDistMatrix!=null){
			//Make a deep copy of everything
			backDistMatrix = new ArrayList<List<Matrix>>();
			for (int k=0; k<p.backDistMatrix.size(); k++){
				backDistMatrix.add(new ArrayList<Matrix>());
				for (int i=0; i<p.backDistMatrix.get(k).size(); i++)
					backDistMatrix.get(k).add((Matrix) p.backDistMatrix.get(k).get(i).clone());
			}
		}
		
		listRMSD = new ArrayList<Double>(p.listRMSD);
		listTMscore = new ArrayList<Double>(p.listTMscore);
		
		rmsd = p.rmsd;
		tmScore = p.tmScore;
	}
	
	@Override
	public PoseMSA clone(){
		return new PoseMSA(this);
	}

	@Override
	public String toString() {
		return "PoseMSA [parent=" + parent + super.toString() + "]";
	}

	@Override
	public MultipleAlignment getParent() {
		return parent;
	}

	@Override
	protected List<Atom[]> getAtomArrays() throws StructureAlignmentException {
		return parent.getAtomArrays();
	}

	@Override
	protected List<Block> getBlocks() throws StructureAlignmentException {
		List<Block> blocks = new ArrayList<Block>();
		//Loop through all the BlockSet of the alignment and add the Blocks
		for (int bs=0; bs<parent.getBlockSetNum(); bs++){
			for (int b=0; b<parent.getBlockSets().get(bs).getBlockNum(); b++){
				blocks.add(parent.getBlockSets().get(bs).getBlocks().get(b));
			}
		}
		return blocks;
	}
	
}
