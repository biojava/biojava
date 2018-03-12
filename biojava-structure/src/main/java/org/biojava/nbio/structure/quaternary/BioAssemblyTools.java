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
package org.biojava.nbio.structure.quaternary;

import org.biojava.nbio.structure.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.vecmath.Point3d;





/**
 * @author Peter Rose
 *
 *
 */
public class BioAssemblyTools {


	/**
	 * Checks if the passed in expression is a unary operator expression
	 * Example: (1,2,3) or (1-60) are unary operator expressions
	 *          (1-60)(61-88) is a binary operator expression, representing
	 *          a cartesian product of the two parenthesised lists
	 *
	 * @param expression
	 * @return true if expression is a unary operator expression
	 */
	public static boolean isUnaryExpression(String expression) {
		int first = expression.indexOf("(");
		int last = expression.lastIndexOf("(");
		if (first < 0 || last < 0) {
			return true;
		}
		return ! (first == 0 && last > first);
	}

	public static List<String> parseUnaryOperatorExpression(String operatorExpression) throws IllegalArgumentException {
		return parseSubExpression(operatorExpression);
	}

	private static List<String> parseSubExpression(String expression) throws IllegalArgumentException {
		// remove parenthesis, if any
		String tmp = expression.replace("(", "");
		tmp = tmp.replace(")", "");

		// separate the operators
		List<String> components = null;
		try {
			components = Arrays.asList(tmp.split(","));
		} catch (Exception e) {
			throw new IllegalArgumentException("Invalid oper_expression: " + expression);
		}

		// expand ranges if present, i.e. 1-60 -> 1, 2, 3, ..., 60
		List<String> operators = new ArrayList<String>();
		for (String component : components) {
			if (component.contains("-")) {
				operators.addAll(expandRange(component));
			} else {
				operators.add(component);
			}
		}
		return operators;
	}

	/**
	 * Expands a range expression, i.e. (1-6) to a list 1,2,3,4,5,6
	 * @param expression the expression to be expanded
	 * @return list of items in range
	 * @throws IllegalArgumentException
	 */
	private static List<String> expandRange(String expression) throws IllegalArgumentException {
		int first = 0;
		int last = 0;
		try {
			String[] range = expression.split("-");
			first = Integer.parseInt(range[0]);
			last = Integer.parseInt(range[1]);
		} catch (Exception e) {
			throw new IllegalArgumentException("Invalid range specification in oper_expression: " + expression);
		}

		List<String> expandedExpression = new ArrayList<String>(last-first+1);
		for (int i = first; i <= last; i++) {
			expandedExpression.add(String.valueOf(i));
		}
		return expandedExpression;
	}

	public static List<OrderedPair<String>> parseBinaryOperatorExpression(String expression)
			throws IllegalArgumentException {
		// split operator expression, i.e. (1,2,3)(4,5) into two subexpressions
		String[] subExpressions = null;
		try {
			subExpressions = expression.split("\\)\\(");
		} catch (Exception e) {
			throw new IllegalArgumentException("Invalid oper_expression: " + expression);
		}
		if (subExpressions.length != 2) {
			throw new IllegalArgumentException("Invalid oper_expression: " + expression);
		}
		List<String> leftSide = parseSubExpression(subExpressions[0]);
		List<String> rightSide = parseSubExpression(subExpressions[1]);

		// form the cartesian product of the two lists
		CartesianProduct<String> product = new CartesianProduct<String>(leftSide, rightSide);
		return product.getOrderedPairs();
	}

	public static double[][]  getBiologicalMoleculeBounds(Structure asymStructure,List<BiologicalAssemblyTransformation> transformations) {
		final double[][] coordinateBounds = new double[2][3];
		coordinateBounds[0][0] = Double.MAX_VALUE;  // min x
		coordinateBounds[0][1] = Double.MAX_VALUE;  // min y
		coordinateBounds[0][2] = Double.MAX_VALUE;  // min z
		coordinateBounds[1][0] = Double.MIN_VALUE;  // max x
		coordinateBounds[1][1] = Double.MIN_VALUE;  // max y
		coordinateBounds[1][2] = Double.MIN_VALUE;  // max z

		double[] transformedCoords = new double[3];

		Atom[] atoms = StructureTools.getAllAtomArray(asymStructure);

		for ( Atom a : atoms) {

			Chain c = a.getGroup().getChain();
			String intChainID = c.getId();

			for (BiologicalAssemblyTransformation m: transformations) {
				if ( ! m.getChainId().equals(intChainID))
					continue;
				Point3d coords = a.getCoordsAsPoint3d();
				transformedCoords[0] = coords.x;
				transformedCoords[1] = coords.y;
				transformedCoords[2] = coords.z;

				if (transformedCoords[0] < coordinateBounds[0][0] ) {
					coordinateBounds[0][0] = transformedCoords[0];  // min x
				}

				if (transformedCoords[1] < coordinateBounds[0][1] ) {
					coordinateBounds[0][1] = transformedCoords[1];  // min y
				}

				if (transformedCoords[2] < coordinateBounds[0][2] ) {
					coordinateBounds[0][2] = transformedCoords[2];  // min z
				}

				if (transformedCoords[0] > coordinateBounds[1][0] ) {
					coordinateBounds[1][0] = transformedCoords[0];  // max x
				}

				if (transformedCoords[1] > coordinateBounds[1][1] ) {
					coordinateBounds[1][1] = transformedCoords[1];  // max y
				}

				if (transformedCoords[2] > coordinateBounds[1][2] ) {
					coordinateBounds[1][2] = transformedCoords[2];  // max z
				}
			}
		}
		return coordinateBounds;
	}
	public static double[][] getAtomCoordinateBounds(Structure s){

		Atom[] atoms = StructureTools.getAllAtomArray(s);
		int atomCount = atoms.length;
		final double[][] coordinateBounds = new double[2][3];
		if ( atomCount <= 0 ) {
			return coordinateBounds;
		}

		Atom a = atoms[0];

		coordinateBounds[0][0] = a.getX();  // min x
		coordinateBounds[0][1] = a.getY();  // min y
		coordinateBounds[0][2] = a.getZ();  // min z
		coordinateBounds[1][0] = a.getX();  // max x
		coordinateBounds[1][1] = a.getY();  // max y
		coordinateBounds[1][2] = a.getZ();  // max z

		for ( int i=1; i<atomCount; i++ )
		{
			a =atoms[i];

			if ( a.getX() < coordinateBounds[0][0] ) {
				coordinateBounds[0][0] = a.getX();  // min x
			}

			if ( a.getY() < coordinateBounds[0][1] ) {
				coordinateBounds[0][1] = a.getY();  // min y
			}

			if ( a.getZ() < coordinateBounds[0][2] ) {
				coordinateBounds[0][2] = a.getZ();  // min z
			}

			if ( a.getX() > coordinateBounds[1][0] ) {
				coordinateBounds[1][0] = a.getX();  // max x
			}

			if ( a.getY() > coordinateBounds[1][1] ) {
				coordinateBounds[1][1] = a.getY();  // max y
			}

			if ( a.getZ() > coordinateBounds[1][2] ) {
				coordinateBounds[1][2] = a.getZ();  // max z
			}
		}

		return coordinateBounds;
	}

	/**
	 * Returns the maximum extend of the structure in the x, y, or z direction.
	 * @param structure
	 * @return maximum extend
	 */
	public static double getMaximumExtend( final Structure structure ) {
		double[][] bounds = getAtomCoordinateBounds(structure);
		double xMax = Math.abs(bounds[0][0] - bounds[1][0]);
		double yMax = Math.abs(bounds[0][1] - bounds[1][1]);
		double zMax = Math.abs(bounds[0][2] - bounds[1][2]);
		return Math.max(xMax, Math.max(yMax, zMax));
	}

	/**
	 * Returns the maximum extend of the biological molecule in the x, y, or z direction.
	 * @param structure
	 * @return maximum extend
	 */
	public static double getBiologicalMoleculeMaximumExtend( final Structure structure,List<BiologicalAssemblyTransformation> transformations ) {
		double[][] bounds = getBiologicalMoleculeBounds(structure, transformations);
		double xMax = Math.abs(bounds[0][0] - bounds[1][0]);
		double yMax = Math.abs(bounds[0][1] - bounds[1][1]);
		double zMax = Math.abs(bounds[0][2] - bounds[1][2]);
		return Math.max(xMax, Math.max(yMax, zMax));
	}

	/**
	 * Returns the centroid of the biological molecule.
	 * @param structure
	 * @return centroid
	 * @throws IllegalArgumentException if structure is null
	 */

	public static double[] getBiologicalMoleculeCentroid( final Structure asymUnit,List<BiologicalAssemblyTransformation> transformations ) throws IllegalArgumentException {
		if ( asymUnit == null ) {
			throw new IllegalArgumentException( "null structure" );
		}

		Atom[] atoms = StructureTools.getAllAtomArray(asymUnit);
		int atomCount = atoms.length;
		double[] centroid = new double[3];

		if ( atomCount <= 0 ) {
			return centroid;
		}

		if ( transformations.size() == 0) {
			return Calc.getCentroid(atoms).getCoords();

		}



		int count = 0;
		double[] transformedCoordinate = new double[3];

		for (int i = 0; i < atomCount; i++)
		{
			Atom atom = atoms[i];
			Chain chain = atom.getGroup().getChain();
			String intChainID = chain.getId();


			for (BiologicalAssemblyTransformation m: transformations) {
				if (!  m.getChainId().equals(intChainID))
					continue;

				Point3d coords = atom.getCoordsAsPoint3d();
				transformedCoordinate[0] = coords.x;
				transformedCoordinate[1] = coords.y;
				transformedCoordinate[2] = coords.z;
				m.transformPoint(transformedCoordinate);
				centroid[0] += transformedCoordinate[0];
				centroid[1] += transformedCoordinate[1];
				centroid[2] += transformedCoordinate[2];
				count++;
			}
		}



		centroid[0] /= count;
		centroid[1] /= count;
		centroid[2] /= count;

		return centroid;
	}

	/** 
	 * Reduce a structure to a single-atom representation (e.g. CA atoms)
	 *
	 * @param orig
	 * @return
	 * @since Biojava 4.1.0
	 */
	public static Structure getReducedStructure(Structure orig){
		Structure s = new StructureImpl();
		s.setPDBHeader(orig.getPDBHeader());
		for ( Chain c : orig.getChains()){

			Chain c1 = new ChainImpl();
			c1.setId(c.getId());
			c1.setName(c.getName());
			s.addChain(c1);

			for (Group g : c.getAtomGroups()){

				Atom a = null;
				switch(g.getType()) {
				case AMINOACID:
					a = g.getAtom(StructureTools.CA_ATOM_NAME);
					break;
				case NUCLEOTIDE:
					a = g.getAtom(StructureTools.NUCLEOTIDE_REPRESENTATIVE);
					break;
				default:
					//omit group
				}
				if ( a != null){

					Group g1 = (Group)g.clone();
					g1.clearAtoms();
					g1.addAtom(a);
					c1.addGroup(g1);

				}

			}

		}
		return s;
	}

}
