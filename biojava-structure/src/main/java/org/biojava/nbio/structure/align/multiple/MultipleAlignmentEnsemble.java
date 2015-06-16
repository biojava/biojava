package org.biojava.nbio.structure.align.multiple;

import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.jama.Matrix;

/**
 * A MultipleAlignmentEnsemble is a collection of {@link MultipleAlignment}s that share the same structures (Atoms) 
 * and creation properties (algorithm, version, creation time, etc.).
 * <p>
 * This class is the top level of the hierarchy and allows the storage of a set of alignment alternatives created by
 * a multiple structure alignment algorithm, so that only one object is returned with more than one alignment option.
 * 
 * @author Aleix Lafita
 *
 */
public interface MultipleAlignmentEnsemble extends ScoresCache {

	/**
	 * Creates and returns an identical copy of this ensemble, including a deep
	 * clone of all constituent alignments.
	 * @return MultipleAlignmentEnsemble identical copy of this object.
	 */
	public MultipleAlignmentEnsemble clone();
	
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
	 * Returns a List containing the names of the structures aligned (i.e.: PDB code, SCOP domain, etc.).<p>
	 * The names are structure identifiers of the structures.
	 * They are in the same order as in the alignment Blocks (same index number for same structure).
	 * @return List of String names of the structures
	 * @see #setStructureNames(List)
	 * @see #getAtomArrays()
	 */
	public List<String> getStructureNames();

	/**
	 * Set the List containing the names of the structures aligned (i.e.: PDB code, SCOP domain, etc.).<p>
	 * The names are structure identifiers of the structures.
	 * @param structureNames names of the structures, structure identifiers
	 * @see #getStructureNames()
	 * @see #setAtomArrays(List)
	 */
	public void setStructureNames(List<String> structureNames);

	/**
	 * Get an array of representative atoms for each structure (CA atoms for proteins).
	 * <p>
	 * Atoms should be unrotated. Thus, to obtain a superimposed set of structures,
	 * each atom array should be cloned and then rotated according to the
	 * transformation matrix.
	 * <p>
	 * If atoms have not previously been set using {@link #setAtomArrays(List)},
	 * attempts to load representative atoms based on {@link #getStructureNames()}.
	 * If it fails to load the Atoms it gives a NullPointerException before returning null.
	 * @return List of Atom[].
	 * @see #setAtomArrays(List)
	 */
	public List<Atom[]> getAtomArrays();

	/**
	 * Sets the List of Atom arrays. Every structure has an Atom array associated.
	 * Note that this should be called in conjunction with {@link #setStructureNames(List)}
	 * <p>
	 * Setting the atom arrays to null will cause them to be automatically
	 * regenerated based on {@link #getStructureNames()} during the next call to
	 * {@link #getAtomArrays()}
	 * @param atomArrays the List of representative (C-alpha) atom arrays
	 * @see #getAtomArrays()
	 * @see #setStructureNames(List)
	 */
	public void setAtomArrays(List<Atom[]> atomArrays);

	/**
	 * Return the number of alternative alignments stored in the Ensemble.
	 * @return int number of alternative alignments.
	 * @see #size()
	 */
	public int getAlignmentNum();

	/**
	 * Returns the List containing the interatomic distance Matrix of each structure.
	 * @return List of Matrix interatomic distance matrices.
	 * @see #updateDistanceMatrix()
	 */
	public List<Matrix> getDistanceMatrix();

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
	 * Add a new multiple alignment to the end of the ensemble and set its
	 * ensemble to this.
	 * @param alignment
	 */
	public void addMultipleAlignment( MultipleAlignment alignment);
	
	/**
	 * Returns the number of aligned structures in the MultipleAlignments.
	 * @return int number of aligned structures.
	 * @see #getStructureNames()
	 * @see #getAtomArrays()
	 */
	public int size();

	/**
	 * Returns the io time for this object, in milliseconds.
	 * @return long creation time, or null if unset
	 */
	public Long getIoTime();
	/**
	 * Set the IO time to load this object
	 * @param millis
	 */
	public void setIoTime(Long millis);

	/**public
	 * Returns the running time of the structure alignment calculation, in milliseconds.
	 * @return long running time of the calculation, or null if unset
	 * @see #getIoTime()L
	 */
	public Long getCalculationTime();
	/**
	 * Set the time needed to calculate this alignment
	 * @param millis
	 */
	public void setCalculationTime(Long millis);

	/**
	 * Clear scores and other properties which depend on the specific alignment.
	 * <p>
	 * This can free memory and ensures consistency for cached variables.
	 */
	public void clear();

}
