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
 */
package org.biojava.nbio.structure.symmetry.internal;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.symmetry.internal.CESymmParameters.SymmetryType;

/**
 * Data Structure that stores all the symmetry axis that describe
 * the symmetry of a structure. Generalizes to all types of symmetry,
 * the classic ones (Cn, Dn) and any hierarchical or local symmetries.
 * <p>
 * Hierarchical symmetry can be visualized as a tree, where each level
 * has a fixed branching factor. Each level of the tree is associated
 * with a transformation operator, whose order determines the degree of
 * nodes at that level of the tree. Leaves of the tree implicitly
 * represent aligned repeats (indexed 0 to n-1), so care must be taken to
 * keep external references to the repeats (e.g. rows of a
 * {@link MultipleAlignment} in the same order implied by the tree.
 * <p>
 * Each node of the tree specifies an alignment between those repeats
 * below each of its children. It is also associated with a symmetry axis,
 * which is calculated based on the associated operator as well as any parent
 * operators.
 * It also stores the parts of the structure (symmetric units) involved
 * in each axis, in addition to the way to calculate them.
 * <p>
 * This is intended to provide a general axis support for the multiple
 * repeat alignment optimization and the axis display in Jmol. This
 * object is related to a MultipleAlignment object that defines the
 * symmetric units.
 *
 * @author Aleix Lafita
 * @since 4.2.0
 *
 */
public class SymmetryAxes {
	/*
	 * Implementation note: The tree is a nice explanation and a good image
	 * for developing algorithms, but it is not constructed explicitly.
	 * Instead, we just store the axis and order for each level and reconstruct
	 * which operators apply to a particular leaf based on that leaf's index.
	 */
	/**
	 * List of all symmetry axis. They are sorted from higher to lower
	 * in the symmetry hierarchy, where higher means that they apply
	 * more globally and lower means that they apply to a local region
	 * of the higher axis division.
	 * 
	 * (Operator for each level in the hierarchy)
	 */
	private final List<Matrix4d> axes;

	/**
	 * Degree (branching factor) for each level in the hierarchy. This should
	 * be the order of the corresponding transformation operator
	 */
	private final List<Integer> degrees;

	private final List<SymmetryType> symmTypes;
	/**
	 * Constructor.
	 * Initializes variables only.
	 */
	public SymmetryAxes(){
		axes = new ArrayList<>();
		degrees = new ArrayList<>();
		symmTypes = new ArrayList<>();
	}

	/**
	 * Adds a new axis of symmetry.
	 * The repeats that participate in this axis and their superposition
	 * relation should also be indicated.
	 *
	 * @param axis the new axis of symmetry found
	 * @param superposition repeats participating and superposition relation
	 * @param repeats number of times the transformation is applied to every
	 * 			repeat. index1=repeat, index2=times.
	 * @param division number of parts that this axis divides the structure in
	 *
	 * @throws IllegalArgumentException if the repeat relation is in the
	 * 			wrong format: should be double List of equal sizes.
	 * @deprecated Use {@link #addAxis(Matrix4d, int, SymmetryType)} instead.
	 *  Repeats and Superposition are now inferred automatically.
	 */
	@Deprecated
	public void addAxis(Matrix4d axis, List<List<Integer>> superposition,
			List<Integer> repeats, Integer division) {

		//Check correct format of repeat relations
		if (superposition.size() != 2){
			throw new IllegalArgumentException(
					"Wrong superposition format: should be double List.");
		} else if (superposition.get(0).size() != superposition.get(1).size()){
			throw new IllegalArgumentException(
					"Wrong superposition format: not equal List sizes.");
		}
		// Now ignores superposition & repeats except to guess symmetry type
		SymmetryType type;
		// Closed if superposition has a circular permutation
		List<Integer> superPos1 = superposition.get(1);
		if(superPos1.get(0) > superPos1.get(superPos1.size()-1)) {
			type = SymmetryType.CLOSED;
		} else {
			type = SymmetryType.OPEN;
		}
		this.addAxis(axis,division,type);
	}
	/**
	 * Adds a new axis of symmetry to the bottom level of the tree
	 *
	 * @param axis the new axis of symmetry found
	 * @param order number of parts that this axis divides the structure in
	 * @param type indicates whether the axis has OPEN or CLOSED symmetry
	 */
	public void addAxis(Matrix4d axis, int order, SymmetryType type) {
		if (order < 2) {
			throw new IllegalArgumentException("A symmetry axis should divide a structure in > 2 parts");
		}
		if(type != SymmetryType.OPEN && type != SymmetryType.CLOSED) {
			throw new IllegalArgumentException("Invalid symmetry type. Only OPEN and CLOSED are allowed");
		}

		axes.add(axis);
		degrees.add(order);
		symmTypes.add(type);
	}

	/**
	 * Return a list giving the number of times each axis must be applied
	 * to generate the given repeat.
	 * <P>
	 * For instance, for a D3 case <tt>getAxisCounts(4)</tt> would return [2,0],
	 * indicating that repeat 4 is generated by two applications of the 3-fold
	 * axis followed by 0 appications of the two-fold axis.
	 * 
	 * @param repeat Index of the desired repeat
	 * @return array of the same length as axes giving the number of times
	 *  to apply each axis.
	 */
	private int[] getAxisCounts(int repeat) {
		int[] counts = new int[degrees.size()];
		
		for(int i = counts.length-1; i >= 0; i--) {
			int d = degrees.get(i);
			counts[i] = repeat % d;
			repeat /= d;
		}
		assert repeat == 0 : "Invalid repeat index";
		return counts;
	}

	/**
	 * Inverse of {@link #getAxisCounts(int)}; Calculates the repeat for a
	 * particular number of applications of each axis
	 * @param counts Number of times to apply each axis
	 * @return Repeat index
	 */
	private int getRepeatIndex(int[] counts) {
		int repeat = 0;
		for(int i = 0; i< counts.length; i++) {
			repeat += counts[i]*degrees.get(i);
		}
		return repeat;
	}
	/**
	 * Updates an axis of symmetry, after the superposition changed.
	 *
	 * @param index old axis index
	 * @param newAxis
	 */
	public void updateAxis(Integer index, Matrix4d newAxis){
		axes.set(index, newAxis);
	}

	/**
	 * Return all elementary axes of symmetry of the structure, that is,
	 * the axes stored in the List as unique and from which all the symmetry
	 * axes are constructed.
	 *
	 * @return axes elementary axes of symmetry.
	 */
	public List<Matrix4d> getElementaryAxes(){
		return axes;
	}

	/**
	 * Returns two lists of the same length.
	 * The first gives a list of all repeat indices which are aligned
	 * at the specified level of symmetry (e.g. 0 through the degree of this level).
	 * The second list gives the corresponding repeats after applying the
	 * operator once.
	 *
	 * @param level the axis index
	 * @return the double List of repeat relations, or null if the
	 * 			level is invalid
	 */
	public List<List<Integer>> getRepeatRelation(int level){
		int m = getNumRepeats(level+1);//size of the children
		int d = degrees.get(level); // degree of this node
		int n = m*d; // number of repeats included
		if(symmTypes.get(level) == SymmetryType.OPEN) {
			n -= m; // leave off last child for open symm
		}
		List<Integer> repeats = new ArrayList<>(n);
		List<Integer> equiv = new ArrayList<>(n);
		for(int i=0;i<n;i++) {
			repeats.add(i);
			equiv.add( (i+m)%(m*d) );
		}
		return Arrays.asList(repeats,equiv);
	}


	/**
	 * Return the transformation that needs to be applied to a
	 * repeat in order to superimpose onto repeat 0.
	 *
	 * @param repeat the repeat index
	 * @return transformation matrix for the repeat
	 */
	public Matrix4d getRepeatTransform(int repeat){

		Matrix4d transform = new Matrix4d();
		transform.setIdentity();

		int[] counts = getAxisCounts(repeat);

		for(int t = counts.length-1; t>=0; t--) {
			if( counts[t] == 0 )
				continue;
			Matrix4d axis = new Matrix4d(axes.get(t));
			for(int i=0;i<counts[t];i++) {
				transform.mul(axis);
			}
		}
		return transform;
	}

	/**
	 * Return all symmetry axes of of the structure: the set of axes that
	 * describe all parts of the structure. This combines the elementary
	 * axes to generate all possible axes. The axes are returned in the repeat
	 * degrees.
	 * @return axes all symmetry axes of the structure.
	 */
	public List<Matrix4d> getSymmetryAxes(){

		List<Matrix4d> symmAxes = new ArrayList<Matrix4d>();

		Matrix4d prior = new Matrix4d();
		prior.setIdentity();
		
		getSymmetryAxes(symmAxes,prior,0);
		
		
		return symmAxes;
	}
	/**
	 * Recursive helper
	 * @param symmAxes output list
	 * @param prior transformation aligning the first repeat of this axis with the first overall
	 * @param level current level
	 */
	private void getSymmetryAxes(List<Matrix4d> symmAxes, Matrix4d prior, int level) {
		if(level >= degrees.size() ) {
			return;
		}
		
		Matrix4d elementary = axes.get(level);

		// Current axis:
		// elementary maps B -> A
		// prior maps I -> A and J -> B
		// want J -> I = J -> B -> A <- I= inv(prior) * elementary * prior
		Matrix4d invPrior = new Matrix4d(prior);
		invPrior.invert();
		Matrix4d currAxis = new Matrix4d(prior);
		currAxis.mul(elementary);
		Matrix4d newPrior = new Matrix4d(currAxis);//save intermediate for later
		currAxis.mul(invPrior);
		symmAxes.add(currAxis);
		//New prior is elementary^d*prior
		//Remember that all degrees are at least 2
		getSymmetryAxes(symmAxes,prior,level+1);
		getSymmetryAxes(symmAxes,newPrior,level+1);
		for(int d=2;d<degrees.get(level);d++) {
			newPrior.mul(elementary);
			getSymmetryAxes(symmAxes,newPrior,level+1);
		}
	}
	
//	public Matrix4d getSymmetryAxis(int level, int axisNum) {
//		if(level == 0) {
//			if( axisNum != 0 )
//				throw new IndexOutOfBoundsException("Axis number out of bounds");
//			return axes.get(0);
//		} else {
//			if( axisNum >= degrees.get(level-1) )
//				throw new IndexOutOfBoundsException("Axis number out of bounds");
//			// Convert axisNum into a count of 
//		
//	}
	/**
	 * Get the number of repeats. This is equal to the product of all degrees.
	 * @return Number of repeats (leaves of the tree).
	 */
	public int getNumRepeats() {
		return getNumRepeats(0);
	}

	/**
	 * Get the number of leaves from a node at the specified level. This is
	 * equal to the product of all degrees at or below the level.
	 * @param level level of the tree to cut at
	 * @return Number of repeats (leaves of the tree).
	 */
	private int getNumRepeats(int level) {
		int size = 1;
		// Return 1 for illegally high level
		if(level < degrees.size()) {
			for(int order : degrees.subList(level, degrees.size())) {
				size *= order;
			}
		}
		return size;
	}
	
	public SymmetryType getSymmetryType(int level) {
		return symmTypes.get(level);
	}
	
	public Matrix4d getElementaryAxis(int level) {
		return axes.get(level);
	}
	
	public int getNumLevels() {
		return degrees.size();
	}
}
