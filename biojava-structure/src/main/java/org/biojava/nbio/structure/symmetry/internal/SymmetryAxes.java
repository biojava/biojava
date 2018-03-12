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
import java.util.Iterator;
import java.util.List;

import javax.vecmath.Matrix4d;

import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.util.RotationAxis;
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
	 * Instead, we just store one elementary axis for each level and reconstruct
	 * which operators apply to a particular leaf based on that leaf's index.
	 */
	
	/**
	 * Represents an axis of symmetry
	 * @author Spencer Bliven
	 *
	 */
	public static class Axis {
		private Matrix4d operator;
		private int order;
		private SymmetryType symmType;
		private int level;
		//private int indexInLevel;
		private int firstRepeat;
		private RotationAxis rotAxis;

		public Axis(Matrix4d operator, int order, SymmetryType type, int level, int firstRepeat) {
			if (order < 2) {
				throw new IllegalArgumentException("A symmetry axis should divide a structure in > 2 parts");
			}
			if(type != SymmetryType.OPEN && type != SymmetryType.CLOSED) {
				throw new IllegalArgumentException("Invalid symmetry type. Only OPEN and CLOSED are allowed");
			}

			this.operator = operator;
			this.order = order;
			this.symmType = type;
			setLevel(level);
			setFirstRepeat(firstRepeat);
			rotAxis = null;
		}
		/**
		 * Get the transformation operator for this axis as an homogeneous matrix
		 * @return the transformation operator
		 */
		public Matrix4d getOperator() {
			return operator;
		}
		public void setOperator(Matrix4d op) {
			this.operator = op;
		}
		/**
		 * Get the order of this axis (closed symm) or the number of repeats
		 * (open symm)
		 * @return the order
		 */
		public int getOrder() {
			return order;
		}
		/**
		 * @return the symmType (OPEN or CLOSED only)
		 */
		public SymmetryType getSymmType() {
			return symmType;
		}
		
		/**
		 * Get the transformation operator as a rotation axis. For open
		 * symmetry this will have a non-zero screw component.
		 * @return a RotationAxis for this Axis
		 */
		public RotationAxis getRotationAxis() {
			if( rotAxis == null) {
				rotAxis = new RotationAxis(operator);
			}
			return rotAxis;
		}
		/**
		 * @return The level of this axis within it's parent hierarchy, or -1 if unset
		 */
		public int getLevel() {
			return level;
		}
		/**
		 * 
		 * @param level The level of this axis within it's parent hierarchy. Must be positive
		 */
		public void setLevel(int level) {
			if(level < 0) throw new IndexOutOfBoundsException("Level must be positive");
			this.level = level;
		}
//		/**
//		 * Each level can contain multiple equivalent axes. This index is
//		 * used to distinguish them.
//		 * @return the index of this axis relative to others at the same level
//		 */
//		public int getIndexInLevel() {
//			return indexInLevel;
//		}
//		/**
//		 * 
//		 * @param indexInLevel the index of this axis relative to others at the same level
//		 */
//		public void setIndexInLevel(int indexInLevel) {
//			if( indexInLevel < 0 || getOrder() <= indexInLevel )
//				throw new IndexOutOfBoundsException("Invalid index for order "+getOrder());
//			this.indexInLevel = indexInLevel;
//		}
		/**
		 * Get the index of the first repeat used by this axis
		 * @return the firstRepeat
		 */
		public int getFirstRepeat() {
			return firstRepeat;
		}
		/**
		 * @param firstRepeat the index of the first repeat used by this axis
		 */
		public void setFirstRepeat(int firstRepeat) {
			this.firstRepeat = firstRepeat;
		}
	}
	
	/**
	 * List of all symmetry axis. They are sorted from higher to lower
	 * in the symmetry hierarchy, where higher means that they apply
	 * more globally and lower means that they apply to a local region
	 * of the higher axis division.
	 */
	private final List<Axis> axes;

	/**
	 * Constructor.
	 * Initializes variables only.
	 */
	public SymmetryAxes(){
		axes = new ArrayList<>();
	}

	/**
	 * Adds a new axis of symmetry to the bottom level of the tree
	 *
	 * @param axis the new axis of symmetry found
	 * @param order number of parts that this axis divides the structure in
	 * @param type indicates whether the axis has OPEN or CLOSED symmetry
	 */
	public void addAxis(Matrix4d axis, int order, SymmetryType type) {
		axes.add(new Axis(axis,order,type,axes.size(),0));
	}

	/**
	 * Return a list giving the number of times each axis must be applied
	 * to generate the given repeat.
	 * <P>
	 * For instance, for a D3 case <tt>getAxisCounts(4)</tt> would return [2,0],
	 * indicating that repeat 4 is generated by two applications of the 3-fold
	 * axis followed by 0 applications of the two-fold axis.
	 * 
	 * @param repeat Index of the desired repeat
	 * @return array of the same length as axes giving the number of times
	 *  to apply each axis.
	 */
	private int[] getAxisCounts(int repeat) {
		int[] counts = new int[getNumLevels()];
		
		for(int i = counts.length-1; i >= 0; i--) {
			int d = axes.get(i).getOrder();
			counts[i] = repeat % d;
			repeat /= d;
		}
		assert repeat == 0 : "Invalid repeat index";
		return counts;
	}

//	/**
//	 * Inverse of {@link #getAxisCounts(int)}; Calculates the repeat for a
//	 * particular number of applications of each axis
//	 * @param counts Number of times to apply each axis
//	 * @return Repeat index
//	 */
//	private int getRepeatIndex(int[] counts) {
//		int repeat = 0;
//		for(int i = 0; i< counts.length; i++) {
//			repeat += counts[i]*axes.get(i).getOrder();
//		}
//		return repeat;
//	}
	/**
	 * Updates an axis of symmetry, after the superposition changed.
	 *
	 * @param index old axis index
	 * @param newAxis
	 */
	public void updateAxis(Integer index, Matrix4d newAxis){
		axes.get(index).setOperator(newAxis);
	}

	/**
	 * Return the operator for all elementary axes of symmetry of the structure, that is,
	 * the axes stored in the List as unique and from which all the symmetry
	 * axes are constructed.
	 *
	 * @return axes elementary axes of symmetry.
	 */
	public List<Matrix4d> getElementaryAxes(){
		List<Matrix4d> ops = new ArrayList<Matrix4d>(getNumLevels());
		for(Axis axis : axes) {
			ops.add(axis.getOperator());
		}
		return ops;
	}
	
	/**
	 * Return all elementary axes of symmetry of the structure, that is,
	 * the axes stored in the List as unique and from which all the symmetry
	 * axes are constructed.
	 *
	 * @return axes elementary axes of symmetry.
	 */
	public List<Axis> getElementaryAxesObjects() {
		return axes;
	}

	/**
	 * Get the indices of participating repeats in Cauchy two-line form.
	 * <p>
	 * Returns two lists of the same length.
	 * The first gives a list of all repeat indices which are aligned
	 * at the specified level of symmetry (e.g. 0 through the degree of this level).
	 * The second list gives the corresponding repeats after applying the
	 * operator once.
	 *
	 * @param level the axis index
	 * @return the double List of repeat relations, or null if the
	 * 			level is invalid
	 * @see #getRepeatsCyclicForm(int, int) for an equivalent specification with half the memory
	 */
	public List<List<Integer>> getRepeatRelation(int level){
		return getRepeatRelation(level,0);
	}

	public List<List<Integer>> getRepeatRelation(Axis axis){
		return getRepeatRelation(axis.getLevel(),axis.getFirstRepeat());
	}

	public List<List<Integer>> getRepeatRelation(int level, int firstRepeat) {
		Axis axis = axes.get(level);
		int m = getNumRepeats(level+1);//size of the children
		int d = axis.getOrder(); // degree of this node
		int n = m*d; // number of repeats included
		if(firstRepeat % n != 0)
			throw new IllegalArgumentException(String.format("Repeat %d cannot start a block at level %s of this tree",firstRepeat,level));
		if(axis.getSymmType() == SymmetryType.OPEN) {
			n -= m; // leave off last child for open symm
		}
		List<Integer> repeats = new ArrayList<>(n);
		List<Integer> equiv = new ArrayList<>(n);
		for(int i=0;i<n;i++) {
			repeats.add(i+firstRepeat);
			equiv.add( (i+m)%(m*d)+firstRepeat );
		}
		return Arrays.asList(repeats,equiv);

	}

	/**
	 * Get the indices of participating repeats in cyclic form.
	 * <p>
	 * Each inner list gives a set of equivalent repeats and should have length
	 * equal to the order of the axis' operator. 
	 * @param level
	 * @param firstRepeat
	 * @return
	 */
	public List<List<Integer>> getRepeatsCyclicForm(int level, int firstRepeat) {
		Axis axis = axes.get(level);
		int m = getNumRepeats(level+1);//size of the children
		int d = axis.getOrder(); // degree of this node
		int n = m*d; // number of repeats included
		if(firstRepeat % n != 0) {
			throw new IllegalArgumentException(String.format("Repeat %d cannot start a block at level %s of this tree",firstRepeat,level));
		}
		if(axis.getSymmType() == SymmetryType.OPEN) {
			n -= m; // leave off last child for open symm
		}
		
		List<List<Integer>> repeats = new ArrayList<>(m);
		for(int i=0;i<m;i++) {
			List<Integer> cycle = new ArrayList<>(d);
			for(int j=0;j<d;j++) {
				cycle.add(firstRepeat+i+j*m);
			}
			repeats.add(cycle);
		}
		return repeats;
	}
	public List<List<Integer>> getRepeatsCyclicForm(Axis axis) {
		return getRepeatsCyclicForm(axis.getLevel(),axis.getFirstRepeat());
	}
	public List<List<Integer>> getRepeatsCyclicForm(int level) {
		return getRepeatsCyclicForm(level,0);
	}
	public String getRepeatsCyclicForm(Axis axis, List<?> repeats) {
		if(repeats.size() != getNumRepeats()) {
			throw new IllegalArgumentException("Mismatch in the number of repeats");
		}
		return getRepeatsCyclicForm(getRepeatsCyclicForm(axis), repeats);
	}
	public static String getRepeatsCyclicForm(List<List<Integer>> cycleForm, List<?> repeats) {
		StringBuilder str = new StringBuilder();
		for(List<Integer> cycle : cycleForm) {
			str.append("(");
			Iterator<Integer> cycleIt = cycle.iterator();
			str.append(repeats.get(cycleIt.next())); //should be at least one
			while(cycleIt.hasNext()) {
				str.append(";")
				.append(repeats.get( cycleIt.next() ));
			}
			str.append(")");
		}
		return str.toString();
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
			Matrix4d axis = new Matrix4d(axes.get(t).getOperator());
			for(int i=0;i<counts[t];i++) {
				transform.mul(axis);
			}
		}
		return transform;
	}
	
	/**
	 * Return the transformation that needs to be applied to
	 * repeat x in order to superimpose onto repeat y.
	 *
	 * @param x the first repeat index (transformed)
	 * @param y the second repeat index (fixed)
	 * @return transformation matrix for the repeat x
	 */
	public Matrix4d getRepeatTransform(int x, int y){

		Matrix4d transform = new Matrix4d();
		transform.setIdentity();

		int[] iCounts = getAxisCounts(x);
		int[] jCounts = getAxisCounts(y);
		
		int[] counts = new int[iCounts.length];
		for (int k = 0; k < iCounts.length; k++)
			counts[k] = iCounts[k] - jCounts[k];
		
		for(int t = counts.length-1; t>=0; t--) {
			if(counts[t] == 0)
				continue;
			if (counts[t] > 0) {
				Matrix4d axis = new Matrix4d(axes.get(t).getOperator());
				for(int i=0;i<counts[t];i++)
					transform.mul(axis);
			} else if (counts[t] < 0) {
				Matrix4d axis = new Matrix4d(axes.get(t).getOperator());
				axis.invert();
				for(int i=0;i<counts[t];i++)
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
	public List<Axis> getSymmetryAxes(){

		List<Axis> symmAxes = new ArrayList<>();

		Matrix4d prior = new Matrix4d();
		prior.setIdentity();
		
		getSymmetryAxes(symmAxes,prior,0,0);
		
		
		return symmAxes;
	}
	/**
	 * Recursive helper
	 * @param symmAxes output list
	 * @param prior transformation aligning the first repeat of this axis with the first overall
	 * @param level current level
	 */
	private void getSymmetryAxes(List<Axis> symmAxes, Matrix4d prior, int level, int firstRepeat) {
		if(level >= getNumLevels() ) {
			return;
		}

		Axis elem = axes.get(level);
		Matrix4d elemOp = elem.getOperator();

		// Current axis:
		// elementary maps B -> A
		// prior maps I -> A and J -> B
		// want J -> I = J -> B -> A <- I= inv(prior) * elementary * prior
		Matrix4d currAxisOp = new Matrix4d(prior);
		currAxisOp.invert();
		currAxisOp.mul(elemOp);
		currAxisOp.mul(prior);
		Axis currAxis = new Axis(currAxisOp,elem.getOrder(),elem.getSymmType(),level,firstRepeat);
		symmAxes.add(currAxis);
		
		//Remember that all degrees are at least 2
		getSymmetryAxes(symmAxes,prior,level+1,firstRepeat);
		//New prior is elementary^d*prior
		Matrix4d newPrior = new Matrix4d(elemOp);
		newPrior.mul(prior);
		int childSize = getNumRepeats(level+1);
		getSymmetryAxes(symmAxes,newPrior,level+1,firstRepeat+childSize);
		for(int d=2;d<elem.getOrder();d++) {
			newPrior.mul(elemOp,newPrior);
			getSymmetryAxes(symmAxes,newPrior,level+1,firstRepeat+childSize*d);
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
		if(level < getNumLevels()) {
			for(Axis axis : axes.subList(level, getNumLevels())) {
				size *= axis.getOrder();
			}
		}
		return size;
	}
	
	/**
	 * Get the first repeat index of each axis of a specified level.
	 * @param level level of the tree to cut at
	 * @return List of first Repeats of each index, sorted in ascending order
	 */
	public List<Integer> getFirstRepeats(int level) {
		List<Integer> firstRepeats = new ArrayList<Integer>();
		int m = getNumRepeats(level+1); //size of the level
		int d = axes.get(level).getOrder(); //degree of this level
		int n = m*d; // number of repeats included in each axis
		for (int firstRepeat = 0; firstRepeat < getNumRepeats(); firstRepeat+=n)
			firstRepeats.add(firstRepeat);
		return firstRepeats;
	}

	public Axis getElementaryAxis(int level) {
		return axes.get(level);
	}

	public int getNumLevels() {
		return axes.size();
	}

}
