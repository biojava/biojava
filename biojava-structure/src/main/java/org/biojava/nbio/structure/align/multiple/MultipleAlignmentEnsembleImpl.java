package org.biojava.nbio.structure.align.multiple;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A general implementation of a {@link MultipleAlignmentEnsemble}.
 * 
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentEnsembleImpl extends AbstractScoresCache implements MultipleAlignmentEnsemble, Serializable, Cloneable {

	private static final long serialVersionUID = -5732485866623431898L;
	
	//Creation Properties
	String algorithmName;
	String version;
	Long ioTime;
	Long calculationTime;							//running time of the algorithm
	
	//Structure Identifiers
	private List<String> structureNames;  			//names of the structures in PDB or SCOP format
	private List<Atom[]> atomArrays;      			//arrays of atoms for every structure in the alignment
	private List<Matrix> distanceMatrix; 			//A list of n (l*l)-matrices that store the distance between every pair of residues for every protein
													//n=number structures; l=length of the protein
	private List<MultipleAlignment> multipleAlignments;
	
	
	/**
	 * Default Constructor. Empty ensemble, no structures assigned.
	 * @return MultipleAlignmentEnsemble an empty MultipleAlignmentEnsemble instance.
	 */
	public MultipleAlignmentEnsembleImpl(){
		
		algorithmName = null;
		version = "1.0";
		ioTime = null;
		calculationTime = null;
		
		structureNames = null;
		atomArrays = null;
		distanceMatrix = null;
		multipleAlignments = null;
	}
	
	/**
	 * Constructor using structure identifiers.
	 * @param structureNames List of Structure names, that can be parsed by AtomCache.
	 * @return MultipleAlignmentEnsemble a MultipleAlignmentEnsemble instance with the structures.
	 */
	public MultipleAlignmentEnsembleImpl(List<String> structureNames){
		
		this();
		setStructureNames(structureNames);
	}
	
	/**
	 * Copy constructor.
	 * @param e MultipleAlignmentEnsembleImpl to copy.
	 * @return MultipleAlignmentEnsembleImpl identical copy of the input MultipleAlignmentEnsembleImpl.
	 */
	public MultipleAlignmentEnsembleImpl(MultipleAlignmentEnsembleImpl e){
		
		algorithmName = e.algorithmName;
		version = e.version;
		ioTime = e.ioTime;
		calculationTime = e.calculationTime;
		
		atomArrays = null;
		if (e.atomArrays != null){
			//Make a deep copy of everything
			atomArrays = new ArrayList<Atom[]>();
			for (Atom[] array:e.atomArrays){
				Atom[] newArray = new Atom[array.length];
				for (int i=0; i<array.length; i++){
					newArray[i] = (Atom) array[i].clone();
				}
				atomArrays.add(newArray);
			}
		}
		
		distanceMatrix = null;
		if (e.distanceMatrix!=null){
			//Make a deep copy of everything
			distanceMatrix = new ArrayList<Matrix>();
			for (Matrix mat:e.distanceMatrix){
				distanceMatrix.add((Matrix) mat.clone());
			}
		}
		
		multipleAlignments = null;
		if (e.multipleAlignments!=null){
			//Make a deep copy of everything
			multipleAlignments = new ArrayList<MultipleAlignment>();
			for (MultipleAlignment msa:e.multipleAlignments){
				MultipleAlignment newMSA = msa.clone();
				newMSA.setEnsemble(this);  //This automatically adds the newMSA to the multipleAlignments list
			}
		}
		
		structureNames = new ArrayList<String>(e.structureNames);
	}
	
	/**
	 * Constructor from an AFPChain instance. Creates an equivalent pairwise alignment.
	 * @param ensemble parent MultipleAlignmentEnsemble.
	 * @return MultipleAlignment a MultipleAlignment instance part of an MultipleAlignmentEnsemble.
	 * @throws StructureAlignmentException 
	 * @throws StructureException 
	 */
	public MultipleAlignmentEnsembleImpl(AFPChain afpChain, Atom[] ca1, Atom[] ca2) throws StructureAlignmentException, StructureException {
		
		//Copy all the creation and algorithm information
		this();
		setAtomArrays(Arrays.asList(ca1,ca2));
		setStructureNames(Arrays.asList(afpChain.getName1(),afpChain.getName2()));
		setAlgorithmName(afpChain.getAlgorithmName());
		setVersion(afpChain.getVersion());
		setCalculationTime(afpChain.getCalculationTime());
		
		MultipleAlignment alignment = new MultipleAlignmentImpl(this);
		setMultipleAlignments(Arrays.asList((MultipleAlignment) alignment));
		
		//Convert the rotation and translation to a Matrix4D and copy it to the MultipleAlignment
		Matrix4d ident = new Matrix4d();
		ident.setIdentity();
		alignment.setTransformations(Arrays.asList(ident, Calc.getTransformation(afpChain.getBlockRotationMatrix()[0], afpChain.getBlockShiftVector()[0])));
		
		//Create a BlockSet for every block in AFPChain and set its transformation
		List<Block>blocks = new ArrayList<Block>(afpChain.getBlockNum());
		for (int bs=0; bs<afpChain.getBlockNum(); bs++){
			BlockSet blockSet = new BlockSetImpl(alignment);
			Block block = new BlockImpl(blockSet);
			block.getAlignRes().add(new ArrayList<Integer>()); //add the two chains
			block.getAlignRes().add(new ArrayList<Integer>());
			blocks.add(block);
			//Set the transformation (convert as before the rotation and translation to a 4D matrix)
			blockSet.setTransformations(Arrays.asList(ident, Calc.getTransformation(afpChain.getBlockRotationMatrix()[bs], afpChain.getBlockShiftVector()[bs])));
			
			for (int i=0; i<afpChain.getOptLen()[bs]; i++){
				block.getAlignRes().get(0).add(afpChain.getOptAln()[bs][0][i]);
				block.getAlignRes().get(1).add(afpChain.getOptAln()[bs][1][i]);
			}
		}
	}

	
	@Override
	public MultipleAlignmentEnsembleImpl clone() {
		return new MultipleAlignmentEnsembleImpl(this);
	}
	
	@Override
	public String getAlgorithmName() {
		return algorithmName;
	}

	@Override
	public void setAlgorithmName(String algorithmName) {
		this.algorithmName = algorithmName;
	}

	@Override
	public String getVersion() {
		return version;
	}

	@Override
	public void setVersion(String version) {
		this.version = version;
	}

	@Override
	public Long getIoTime() {
		return ioTime;
	}

	@Override
	public void setIoTime(Long millis) {
		this.ioTime = millis;
	}
	
	@Override
	public Long getCalculationTime() {
		return calculationTime;
	}

	@Override
	public void setCalculationTime(Long millis) {
		this.calculationTime = millis;
	}

	@Override
	public List<String> getStructureNames() {
		return structureNames;
	}

	@Override
	public void setStructureNames(List<String> structureNames) {
		this.structureNames = structureNames;
	}

	@Override
	public List<Atom[]> getAtomArrays() throws StructureAlignmentException {
		if (atomArrays == null)
			try {
				updateAtomArrays();
			} catch (IOException e) {
				throw new StructureAlignmentException(e.getMessage(),e);
			} catch (StructureException e) {
				throw new StructureAlignmentException(e.getMessage(),e);
			}
		return atomArrays;
	}

	@Override
	public void setAtomArrays(List<Atom[]> atomArrays) {
		this.atomArrays = atomArrays;
	}
	
	/**
	 * Force the atom arrays to regenerate based on {@link #getStructureNames()}
	 * @throws StructureAlignmentException
	 * @throws IOException
	 * @throws StructureException
	 */
	public void updateAtomArrays() throws IOException, StructureException{
		AtomCache cache = new AtomCache();
		atomArrays = new ArrayList<Atom[]>();
		for (String name : getStructureNames() ){
			Atom[] array = cache.getRepresentativeAtoms(name);
			atomArrays.add(array);
		}
		//TODO update superposition & other properties
	}

	@Override
	public int getAlignmentNum() {
		return multipleAlignments.size();
	}

	@Override
	public List<Matrix> getDistanceMatrix() throws StructureAlignmentException {
		if (distanceMatrix == null) updateDistanceMatrix();
		return distanceMatrix;
	}

	/**
	 * Force recalculation of the distance matrices
	 * @throws StructureAlignmentException
	 */
	public void updateDistanceMatrix() throws StructureAlignmentException {
		
		//Reset the distance Matrix variable
		distanceMatrix = new ArrayList<Matrix>();
		
		for (int s=0; s<size(); s++){
			int n = atomArrays.get(s).length;  //length of the structure
			Matrix distMat = new Matrix(n,n);
			
			//Calculate all distances between every pair of atoms and set the entries
			for (int a1=0; a1<n; a1++){
				for (int a2=0; a2<n; a2++){
					double dist = Calc.getDistance(atomArrays.get(s)[a1], atomArrays.get(s)[a2]);
					distMat.set(a1, a2, dist);
				}
			}
			distanceMatrix.add(distMat);
		}
	}

	@Override
	public List<MultipleAlignment> getMultipleAlignments() {
		if (multipleAlignments == null) multipleAlignments = new ArrayList<MultipleAlignment>();
		return multipleAlignments;
	}

	@Override
	public void setMultipleAlignments(List<MultipleAlignment> multipleAlignments) {
		this.multipleAlignments = multipleAlignments;
	}
	
	/**
	 * Add a new multiple alignment to the end of the ensemble and set its
	 * ensemble to this.
	 * @param alignment
	 */
	@Override
	public void addMultipleAlignment( MultipleAlignment alignment) {
		multipleAlignments.add(alignment);
		alignment.setEnsemble(this);
	}


	@Override
	public int size() throws StructureAlignmentException {
		if (structureNames != null) return structureNames.size();
		else if (atomArrays != null) return atomArrays.size();
		else throw new StructureAlignmentException("Empty MultipleAlignmentEnsemble: structureNames == null");
	}
	
	/**
	 * Clear scores and distance matrix. Recursively clears member alignments.
	 */
	@Override
	public void clear() {
		super.clear();
		distanceMatrix = null;
		for(MultipleAlignment a : getMultipleAlignments()) {
			a.clear();
		}
	}

}
