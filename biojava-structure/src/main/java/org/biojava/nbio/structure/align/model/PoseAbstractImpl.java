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
 * An abstract implementation of a Pose to generalize the shared code of PoseBS ({@link BlockSet}) and PoseMSA ({@link MultipleAlignment}).
 * This Pose abstract class does not have a parent.
 * 
 * @author Aleix Lafita
 * 
 */
public abstract class PoseAbstractImpl implements Serializable, Pose{

	private static final long serialVersionUID = 4630634095792288834L;

	protected List<Matrix> rotation;					//rotation Matrix for every structure to calculate the 3D superimposition.
	protected List<Atom> translation;					//translation Vectors for the atoms of every structure.
	protected List<List<Matrix>> backDistMatrix;    	//Background distances of every structure (dot-plot).
	protected List<Double> listRMSD;					//RMSD of every individual structure
	protected List<Double> listTMscore;				//TMscore of every individual structure
	protected double rmsd;							//RMSD of the 3D superposition
	protected double tmScore;							//TMscore of the 3D superposition

	/**
	 * Constructor.
	 */
	public PoseAbstractImpl() {
		
		rotation = null;
		translation = null;
		backDistMatrix = null;
		listRMSD = null;
		listTMscore = null;
		rmsd = -1;
		tmScore = -1;
	}
	
	@Override
	public abstract Object clone();

	@Override
	public String toString() {
		return "rotation=" + rotation + ", translation="
				+ translation + ", backDistMatrix=" + backDistMatrix
				+ ", listRMSD=" + listRMSD + ", listTMscore=" + listTMscore
				+ ", rmsd=" + rmsd + ", tmScore=" + tmScore;
	}

	@Override
	public abstract Object getParent();

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
		if (getAtomArrays() == null) throw new StructureAlignmentException("Empty MultipleAlignmentEnsemble: atomArrays == null");
		else return getAtomArrays().size();
	}

	@Override
	public double getRMSD() throws StructureAlignmentException {
		if (rmsd == -1) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return rmsd;
	}
	
	@Override
	public double getRMSD(int structureNr) throws StructureAlignmentException {
		if (listRMSD == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return listRMSD.get(structureNr);
	}

	@Override
	public double getTMscore() throws StructureAlignmentException {
		if (tmScore == -1) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return tmScore;
	}
	
	@Override
	public double getTMscore(int structureNr) throws StructureAlignmentException {
		if (listTMscore == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		else return listTMscore.get(structureNr);
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
				
				//Loop through all the Blocks that define the aligned positions
				for (int k=0; k<getBlocks().size(); k++){
					//Loop through all the columns of each Block (aligned residues)
					for (int j=0; j<getBlocks().get(k).length(); j++){
						Integer pos1 = getBlocks().get(k).getAlignRes().get(0).get(j);
						Integer pos2 = getBlocks().get(k).getAlignRes().get(i).get(j);
						
						if (pos1==null || pos2==null) continue;
						atomSet1.add((Atom) getAtomArrays().get(0)[pos1]);
						atomSet2.add((Atom) getAtomArrays().get(i)[pos2].clone());
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
		
		//Once translation and rotation are set, update the other Cache variables
		updateBackDistMatrix();
		updateRMSDandScore();
	}

	@Override
	public List<Atom[]> getRotatedAtoms() throws StructureAlignmentException, StructureException {
		
		if (rotation == null || translation == null) throw new StructureAlignmentException("Empty Pose: updatePose() first.");
		
		List<Atom[]> rotatedAtoms = new ArrayList<Atom[]>();
		//Rotate the atom coordinates of all the structures
		for (int i=0; i<size(); i++){
			Atom[] rotCA = StructureTools.cloneAtomArray(getAtomArrays().get(i));
			Calc.rotate(rotCA, rotation.get(i));
			Calc.shift(rotCA, translation.get(i));
			rotatedAtoms.add(rotCA);
		}
		return rotatedAtoms;
	}
	
	/**
	 * Helper method.
	 * Calculates and sets the background distance Matrices of all structural pairwise comparisons.
	 * The number of comparisons is quatratic in number of structures and quadratic in number of residues.
	 * The distances can be used either for RMSD and TMscore calculation or to calculate algorithm scores to
	 * optimize the alignment.
	 * @throws StructureAlignmentException
	 * @throws StructureException
	 */
	private void updateBackDistMatrix() throws StructureAlignmentException, StructureException {
		
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
						double distance = Calc.getDistance(rotatedAtoms.get(s1)[i], rotatedAtoms.get(s2)[j]);
						backDistMatrix.get(s1).get(s2).set(i, j, distance);
					}
				} //finish Matrix
			}
		} //finish all pairwise comparisons
	}
	
	/**
	 * Helper method.
	 * Calculates and sets the global and individual RMSD and TMscore of the alignment.
	 * Needs the updated background distances matrices, since it assumes all distances to be already calculated.
	 * @throws StructureAlignmentException
	 * @throws StructureException
	 */
	private void updateRMSDandScore() throws StructureAlignmentException, StructureException{
		
		int size = size();  //better performance because it is used a lot.
		rmsd =0;
		tmScore=0;
		listRMSD = new ArrayList<Double>();
		listTMscore = new ArrayList<Double>();
		
		//Loop for every structure (s1) and look at the distances to every other aligned atom of the other structures (s2)
		for (int s1=0; s1<size; s1++){
			int lenRMSD = 0;  //the number of aligned residues to this structure s1.
			int lenScore = 0; //the total number of possible residues aligned of this structure s1.
			double sumRMSD = 0.0;
			double sumScore = 0.0;
			for (int s2=0; s2<size; s2++){
				//Constants for the TM score calculation (taken from SVDSuperimposer)
				int Lmin = Math.min(getAtomArrays().get(s1).length,getAtomArrays().get(s2).length);
				double d0 = 1.24 * Math.cbrt(Lmin - 15.) - 1.8;
				double d0sq = d0*d0;
				lenScore += Lmin;
				
				//Loop through all the Blocks that define the aligned positions
				for (int k=0; k<getBlocks().size(); k++){
					//Loop through all the columns of each Block (aligned residues)
					for (int j=0; j<getBlocks().get(k).length(); j++){
						Integer pos1 = getBlocks().get(k).getAlignRes().get(s1).get(j);
						Integer pos2 = getBlocks().get(k).getAlignRes().get(s2).get(j);
						
						if (pos1==null || pos2==null) continue;
						double d = backDistMatrix.get(s1).get(s2).get(pos1, pos2);
						sumRMSD += (d*d);
						sumScore += 1./(1+d*d/d0sq);   //Formulas taken from SVDSuperimposer class
						lenRMSD++;
					}
				} //end Block comparison
			}
			listRMSD.add(Math.sqrt(sumRMSD/ lenRMSD));
			listTMscore.add(sumScore/lenScore);
		} //end all-to-all comparisons
		
		//The global RMSD is the average RMSD of all the structures 
		//TODO should it be weighted average depending on the aligned residues of each structure? (only important if there are a lot of GAPS between structures)
		for (int i=0; i<size; i++){
			rmsd += listRMSD.get(i);
			tmScore += listTMscore.get(i);
		}
		rmsd /= size;
		tmScore /= size;
	}
	
	/**
	 * Helper method.
	 * Returns the List of Atoms of the structures, without rotation.
	 * The Atoms returned are not a copy of the originals.
	 * @return List of Atom arrays
	 * @throws StructureAlignmentException
	 */
	protected abstract List<Atom[]> getAtomArrays() throws StructureAlignmentException;
	
	/**
	 * Helper method.
	 * Returns the List of alignment Blocks associated with this Pose.
	 * @return List of Blocks
	 * @throws StructureAlignmentException
	 */
	protected abstract List<Block> getBlocks() throws StructureAlignmentException;
	
}
