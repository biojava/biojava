package org.biojava.nbio.structure.align.model;

import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * An implementation of a Pose for a {@link BlockSet}.
 * 
 * @author Aleix Lafita
 * 
 */
public class PoseBS extends PoseAbstractImpl{

	private static final long serialVersionUID = -8309408466325388360L;
	private BlockSet parent;
	
	/**
	 * Constructor.
	 * @param blockSet the parent BlockSet of the PoseBS instance. Cannot be null.
	 * @return PoseBS a PoseBS instance linked to its parent BlockSet.
	 * @throws StructureAlignmentException if the parent is null.
	 */
	public PoseBS(BlockSet blockSet) throws StructureAlignmentException {
		super();
		
		if (blockSet == null) throw new StructureAlignmentException("Parent of Pose is null.");
		parent = blockSet;
	}
	
	/**
	 * Copy constructor.
	 * @param p PoseBS object to be copied.
	 * @return PoseBS an identical copy of the input PoseBS object.
	 */
	public PoseBS(PoseBS p) {
		
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
	public Object clone(){
		return new PoseBS(this);
	}

	@Override
	public String toString() {
		return "PoseBS [parent=" + parent + super.toString() + "]";
	}

	@Override
	public BlockSet getParent() {
		return parent;
	}

	@Override
	protected List<Atom[]> getAtomArrays() throws StructureAlignmentException {
		return parent.getMultipleAlignment().getAtomArrays();
	}

	@Override
	protected List<Block> getBlocks() throws StructureAlignmentException {
		return parent.getBlocks();
	}
	
}
