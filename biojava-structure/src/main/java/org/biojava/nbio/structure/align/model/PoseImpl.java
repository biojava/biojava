package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A general implementation of a Pose.
 * 
 * @author Aleix Lafita
 * 
 */
public class PoseImpl implements Serializable, Pose{

	private static final long serialVersionUID = -8309408466325388360L;
	
	private BlockSet parent;

	private List<Matrix> rotation;					//rotation Matrix for every structure to calculate the 3D superimposition.
	private List<Atom> translation;					//translation Vectors for the atoms of every structure.
	private List<List<Matrix>> backDistMatrix;    	//Background distances of every structure (dot-plot). 
														//TODO In theory there would be N^2 distance matrices, one for each pairwise comparison. Is that feasible? Or too much space and not scalable?
	private double rmsd;							//RMSD of the 3D superposition
	private double tmScore;							//score of the 3D superposition

	/**
	 * Constructor.
	 * @param blockSet the parent BlockSet of the PoseImpl instance.
	 * @return PoseImpl a BlockImpl instance linked to its parent BlockSet.
	 */
	public PoseImpl(BlockSet blockSet) {
		
		parent = blockSet;
		rotation = null;
		translation = null;
		backDistMatrix = null;
		rmsd = 0;
		tmScore = 0;
	}
	
	/**
	 * Copy constructor.
	 * @param p PoseImpl object to be copied.
	 * @return PoseImpl an identical copy of the input PoseImpl object.
	 */
	public PoseImpl(PoseImpl p) {
		
		parent = p.parent;
		
		rotation = null;
		if (p.getRotation()!=null){
			//Make a deep copy of everything
			rotation = new ArrayList<Matrix>();
			for (int i=0; i<p.getRotation().size(); i++)
				rotation.add((Matrix) p.getRotation().get(i).clone());
		}
		
		translation = null;
		if (p.getTranslation()!=null){
			//Make a deep copy of everything
			translation = new ArrayList<Atom>();
			for (int i=0; i<p.getRotation().size(); i++)
				translation.add((Atom) p.getTranslation().get(i).clone());
		}
		
		backDistMatrix = null;
		if (p.getBackDistMatrix()!=null){
			//Make a deep copy of everything
			backDistMatrix = new ArrayList<List<Matrix>>();
			for (int k=0; k<p.getBackDistMatrix().size(); k++){
				backDistMatrix.add(new ArrayList<Matrix>());
				for (int i=0; i<p.getBackDistMatrix().get(k).size(); i++)
					backDistMatrix.get(k).add((Matrix) p.getBackDistMatrix().get(k).get(i).clone());
			}
		}
	}
	
	@Override
	public Object clone(){
		
		return new PoseImpl(this);
	}

	@Override
	public String toString() {
		return "PoseImpl [parent=" + parent + ", rotation=" + rotation
				+ ", translation=" + translation + ", backDistMatrix="
				+ backDistMatrix + "]";
	}

	@Override
	public void setBlockSet(BlockSet parent) {
		this.parent = parent;
	}

	@Override
	public BlockSet getBlockSet() {
		return parent;
	}

	@Override
	public List<Matrix> getRotation() {
		return rotation;
	}

	@Override
	public void setRotation(List<Matrix> rotation) {
		this.rotation = rotation;
	}

	@Override
	public List<Atom> getTranslation() {
		return translation;
	}

	@Override
	public void setTranslation(List<Atom> translation) {
		this.translation = translation;
	}

	@Override
	public List<List<Matrix>> getBackDistMatrix() {
		return backDistMatrix;
	}

	@Override
	public void setBackDistMatrix(List<List<Matrix>> backDistMatrix) {
		this.backDistMatrix = backDistMatrix;
	}

	@Override
	public int size() {
		//Get the size from the variables that can contain the information
		if (parent != null) return parent.size();
		else if (backDistMatrix != null) return backDistMatrix.size();
		else if (rotation != null) return rotation.size();
		else if (translation != null) return translation.size();
		else return 0;
	}

	@Override
	public double getRMSD() {
		return rmsd;
	}

	@Override
	public void updateRMSD() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double getTMscore() {
		return tmScore;
	}

	@Override
	public void updateTMscore() {
		// TODO Auto-generated method stub
		
	}
	
}
