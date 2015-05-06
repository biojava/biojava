package org.biojava.nbio.structure.align.model;

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A MultipleAlignmentEnsemble is a collection of {@link MultipleAlignments} that share the same Atom arrays (structures) and 
 * creation properties (algorithm, version, creation time, etc.).
 * This class deals with the multiple alternatives returned by the MSTA algorithms, so that only one object is returned
 * with more than one alignment option (depending on the RMSD/length trade-offs made, for example).
 * 
 * @author Aleix Lafita
 *
 */
public interface MultipleAlignmentEnsemble extends Cloneable{

	/**
	 * Creates and returns an identical copy of this object.
	 * @return EnsembleMSTA identical copy of this object.
	 */
	public Object clone();
	
	/**
	 * Returns the name of the multiple structure alignment algorithm that created the MultipleAlignment objects.
	 * @return String name of the algorithm.
	 * @see #setAlgorithmName(String)
	 */
	public String getAlgorithmName();

	/**
	 * Set the name of the multiple structure alignment algorithm that created the MultipleAlignment objects.
	 * @param algorithmName name of the algorithm.
	 * @see #getAlgorithmName()
	 */
	public void setAlgorithmName(String algorithmName);

	/**
	 * Returns the version of the algorithm used to generate the MultipleAlignment objects.
	 * @return String version of the algorithm.
	 * @see #setVersion(String)
	 */
	public String getVersion();

	/**
	 * Sets the version of the algorithm used to generate the MultipleAlignment objects.
	 * @param version the version of the algorithm.
	 * @see #getVersion()
	 */
	public void setVersion(String version);

	/**
	 * Returns the creation time of this object, in milliseconds.
	 * @return long creation time.
	 */
	public long getIoTime();

	/**
	 * Returns the running time of the structure alignment calculation, in milliseconds.
	 * @return long running time of the calculation.
	 * @see #setCalculationTime(long)
	 */
	public long getCalculationTime();

	/**
	 * Set the running time of the structure alignment calculation, in milliseconds.
	 * @param calculationTime running time of the calculation.
	 * @see #getCalculationTime()
	 */
	public void setCalculationTime(long calculationTime);

	/**
	 * Returns the structure alignment object ID.
	 * @return long structure alignment ID.
	 * @see #setId(long)
	 */
	public long getId();

	/**
	 * Sets the structure alignment object ID.
	 * @param id structure alignment object ID.
	 * @see #getId()
	 */
	public void setId(long id);
	
	/**
	 * Returns a List containing the names of the structures aligned (i.e.: PDB code, SCOP domain, etc.).
	 * The names are structure identifiers of the structures.
	 * They are in the same order as in the alignment Blocks (same index number for same structure).
	 * @return List of String names of the structures
	 * @see #setStructureNames(List)
	 * @see #getAtomArrays()
	 */
	public List<String> getStructureNames();

	/**
	 * Set the List containing the names of the structures aligned (i.e.: PDB code, SCOP domain, etc.).
	 * The names are structure identifiers of the structures.
	 * @param structureNames names of the structures, structure identifiers
	 * @see #getStructureNames()
	 * @see #setAtomArrays(List)
	 */
	public void setStructureNames(List<String> structureNames);

	/**
	 * Returns the List of Atom arrays. Every structure has an Atom array associated.
	 * The Atom arrays are only stored as a cache, and must be deleted when the alignment is serialized or stored.
	 * @return List of Atom[].
	 * @see #setAtomArrays(List)
	 */
	public List<Atom[]> getAtomArrays();

	/**
	 * Sets the List of Atom arrays. Every structure has an Atom array associated.
	 * The Atom arrays are only stored as a cache, and must be deleted when the alignment is serialized or stored.
	 * @param atomArrays the List of Atom[].
	 * @see #getAtomArrays()
	 * @see #setStructureNames(List)
	 */
	public void setAtomArrays(List<Atom[]> atomArrays);
	
	/**
	 * Downloads and sets the List of Atom arrays from the Structure identifiers.
	 * The Atom arrays are only stored as a cache, and must be deleted when the alignment is serialized or stored.
	 * @param atomArrays the List of Atom[].
	 * @throws StructureAlignmentException
	 * @throws StructureException 
	 * @throws IOException 
	 * @see #setAtomArrays(List)
	 * @see #setStructureNames(List)
	 */
	public void updateAtomArrays() throws StructureAlignmentException, IOException, StructureException;
	
	/**
	 * Return the number of alternative alignments stored in the EnsembleMSTA object.
	 * @return int number of alternative alignments.
	 * @see #size()
	 */
	public int getAlignmentNum();

	/**
	 * Returns the List containing the interatomic distance Matrix of each structure.
	 * @return List of Matrix interatomic distance matrices.
	 * @throws StructureAlignmentException 
	 * @see #updateDistanceMatrix()
	 */
	public List<Matrix> getDistanceMatrix() throws StructureAlignmentException;

	/**
	 * Calculates and sets the cache List containing the interatomic distance Matrix of each structure.
	 * @throws StructureAlignmentException 
	 * @see #getDistanceMatrix()
	 */
	public void updateDistanceMatrix() throws StructureAlignmentException;
	
	/**
	 * Returns the List of MultipleAlignments in the MultipleAlignmentEnsemble object.
	 * @return List of MultipleAlignment in the MultipleAlignmentEnsemble.
	 * @see #setMultipleAlignments()
	 * @see #getOptimalMultipleAlignment()
	 */
	public List<MultipleAlignment> getMultipleAlignments();
	
	/**
	 * Set the List of MultipleAlignments in the MultipleAlignmentEnsemble object.
	 * @param List of MultipleAlignments that are part of the ensemble.
	 * @see #getMultipleAlignments()
	 * @see #getOptimalMultipleAlignment()
	 */
	public void setMultipleAlignments(List<MultipleAlignment> multipleAlignments);
	
	/**
	 * Returns the optimal MultipleAlignment of the ensemble, stored at the first position (index 0) of the List.
	 * @return MultipleAlignment optimal MSTA of the MultipleAlignmentEnsemble.
	 * @throws StructureAlignmentException if the MultipleAlignmentEnsemble is empty.
	 * @see #getMultipleAlignments()
	 * @see #setMultipleAlignments()
	 */
	public MultipleAlignment getOptimalMultipleAlignment() throws StructureAlignmentException;
	
	/**
	 * Returns the number of aligned structures in the MultipleAlignments.
	 * @return int number of aligned structures.
	 * @throws StructureAlignmentException if atomArrays is null.
	 * @see #getStructureNames()
	 * @see #getAtomArrays()
	 */
	public int size() throws StructureAlignmentException;
	
}
