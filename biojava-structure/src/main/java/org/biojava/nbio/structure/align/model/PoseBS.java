package org.biojava.nbio.structure.align.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.SVDSuperimposer;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * An implementation of a Pose for a BlockSet.
 * 
 * @author Aleix Lafita
 * 
 */
public class PoseBS implements Serializable, Pose{

	private static final long serialVersionUID = -8309408466325388360L;
	
	private BlockSet parent;

	private List<Matrix> rotation;					//rotation Matrix for every structure to calculate the 3D superimposition.
	private List<Atom> translation;					//translation Vectors for the atoms of every structure.
	private List<List<Matrix>> backDistMatrix;    	//Background distances of every structure (dot-plot).
	private double rmsd;							//RMSD of the 3D superposition
	private double tmScore;							//score of the 3D superposition

	/**
	 * Constructor.
	 * @param blockSet the parent BlockSet of the PoseBS instance. Cannot be null.
	 * @return PoseBS a PoseBS instance linked to its parent BlockSet.
	 * @throws StructureAlignmentException if the parent is null.
	 */
	public PoseBS(BlockSet blockSet) throws StructureAlignmentException {
		
		if (blockSet == null) throw new StructureAlignmentException("Parent of Pose is null.");
		parent = blockSet;
		
		rotation = null;
		translation = null;
		backDistMatrix = null;
		
		rmsd = -1;
		tmScore = -1;
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
		
		rmsd = p.rmsd;
		tmScore = p.tmScore;
	}
	
	@Override
	public Object clone(){
		return new PoseBS(this);
	}

	@Override
	public String toString() {
		return "PoseBS [parent=" + parent + ", rotation=" + rotation
				+ ", translation=" + translation + ", backDistMatrix="
				+ backDistMatrix + ", rmsd=" + rmsd + ", tmScore=" + tmScore
				+ "]";
	}

	@Override
	public BlockSet getParent() {
		return parent;
	}

	@Override
	public List<Matrix> getRotation() throws StructureAlignmentException {
		if (rotation == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return rotation;
	}

	@Override
	public List<Atom> getTranslation() throws StructureAlignmentException {
		if (translation == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return translation;
	}

	@Override
	public List<List<Matrix>> getBackDistMatrix() throws StructureAlignmentException {
		if (rotation == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return backDistMatrix;
	}

	@Override
	public int size() throws StructureAlignmentException {
		return parent.size(); //Get the size from the parent
	}

	@Override
	public double getRMSD() throws StructureAlignmentException {
		if (rmsd == -1) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return rmsd;
	}

	@Override
	public double getTMscore() throws StructureAlignmentException {
		if (tmScore == -1) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return tmScore;
	}
	
	@Override
	public void updatePose(PoseMethod method) throws StructureException, StructureAlignmentException{
		
		//Initialize or replace the rotation and translation variables
		rotation = new ArrayList<Matrix>();
		translation = new ArrayList<Atom>();
		
		switch (method) {
		case REFERENCE:
			//We suppose the first molecule as reference and superimpose everything to it
			for (int i=0; i<size(); i++){
				List<Atom> atomSet1 = new ArrayList<Atom>();
				List<Atom> atomSet2 = new ArrayList<Atom>();
				for (int k=0; k<parent.getBlockNum(); k++){
					for (int j=0; j<parent.getBlocks().get(k).length(); j++){
						Integer pos1 = parent.getBlocks().get(k).getAlignRes().get(0).get(j);
						Integer pos2 = parent.getBlocks().get(k).getAlignRes().get(i).get(j);
						
						if (pos1==null || pos2==null) continue;
						atomSet1.add((Atom) parent.getMultipleAlignment().getAtomArrays().get(0)[pos1]);
						atomSet2.add((Atom) parent.getMultipleAlignment().getAtomArrays().get(i)[pos2].clone());
					}
				}
				Atom[] array1 = atomSet1.toArray(new Atom[atomSet1.size()]);
				Atom[] array2 = atomSet2.toArray(new Atom[atomSet2.size()]);
				
				//From the superimposer we obtain the rotation and translation information
				SVDSuperimposer svd = new SVDSuperimposer(array1, array2);
				rotation.add(svd.getRotation());
				translation.add(svd.getTranslation());
			}
			break;
		case MEDIAN:
			//TODO implement this type of superimposition
			System.out.println("Not yet implemented!");
			break;
		case CONSENSUS:
			//TODO implement this type of superimposition
			System.out.println("Not yet implemented!");
			break;
		}
		
		//Once translation and rotation are set, update the other Cache variables (in that order)
		updateBackDistMatrix();
		updateRMSD();
		updateTMscore();
	}

	private void updateTMscore() {
		// TODO Auto-generated method stub
		
	}

	private void updateRMSD() {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public void updateBackDistMatrix() throws StructureAlignmentException, StructureException {
		
		//Reset or initialize the backDistanceMatrices
		backDistMatrix = new ArrayList<List<Matrix>>();
		
		//First rotate the atoms
		List<Atom[]> rotatedAtoms = getRotatedAtoms();
		//Loop thorugh all pairwise structure comparisons
		for (int s1=0; s1<size(); s1++){
			int n = rotatedAtoms.get(s1).length;
			backDistMatrix.add(new ArrayList<Matrix>());
			
			for (int s2=0; s2<size(); s2++){
				int m = rotatedAtoms.get(s2).length;
				backDistMatrix.get(s1).add(new Matrix(n, m));
				if (s1==s2) continue;
				
				//Loop through all residue pairwise combinations and complete the Matrix
				for (int i=0; i<n; i++){
					for (int j=0; j<m; j++){
						double distance = Calc.getDistance(rotatedAtoms.get(s1)[n], rotatedAtoms.get(s2)[m]);
						backDistMatrix.get(s1).get(s2).set(i, j, distance);
					}
				}
			}
		}
	}

	@Override
	public List<Atom[]> getRotatedAtoms() throws StructureAlignmentException, StructureException {
		
		if (rotation == null || translation == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		
		List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();
		//Rotate the atom coordinates of all the structures
		for (int i=0; i<size(); i++){
			Atom[] rotCA = StructureTools.cloneAtomArray(parent.getMultipleAlignment().getAtomArrays().get(i));
			Calc.rotate(rotCA, rotation.get(i));
			Calc.shift(rotCA, translation.get(i));
			rotatedAtoms.add(rotCA);
		}
		return rotatedAtoms;
	}
	
}
