package org.biojava.bio.structure.quaternary;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.biojava.bio.structure.jama.Matrix;
import org.biojava3.core.util.PrettyXMLWriter;
import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;



/**
 * Isolation issue -
 * 
 * The original version of this (renamed as GLTransformationMatrix) uses FloatBuffer objects obtained from
 * GL to store and operate on the matrices - doing so destroys the ability
 * to isolate the model, because the loaders use this to generate and store biological unit and non-crystallographic
 * transformations.
 * 
 * This version is a reimplementation that uses a simple array of floats.  These are almost always stored in
 * lists - note GLTransformationList has a utility to return a list of FloatBuffer objects for use in
 * gl rendering.
 * 
 * TODO: re-implement this as a subclass of Matrix3x and remove redundant functionality.
 *       Revert the array implementation to a Vector3x/Point3x implementation.
 *       All of the vertices and points should be carried consistently.
 *       24-Nov-08 - rickb
 * 
 * @see org.rcsb.mbt.glscene.jogl.GLTransformationMatrix
 * @see org.rcsb.vf.glscene.jogl.GLTransformationList
 * @see org.rcsb.mbt.model.geometry.ModelTransformationList
 * 
 * @author rickb
 *
 */
public class ModelTransformationMatrix implements Cloneable {
	public String id = null;

	public String ndbChainId = null;

	public String symmetryShorthand = null;

	public String code = null;

	public float values[];
	
	/**
	 * Default Constructor
	 */
	public ModelTransformationMatrix()
	{
		init();
	}

	/**
	 * Copy Constructor
	 * 
	 * @param src
	 */
	public ModelTransformationMatrix(final ModelTransformationMatrix src)
	{
		init();
		for (int ix = 0; ix < 16; ix++)
			values[ix] = src.values[ix];

		this.id = src.id;
		//this.cell = src.cell;
		this.ndbChainId = src.ndbChainId;
		this.symmetryShorthand = src.symmetryShorthand;
		this.code = src.code;
	}

	
	public void setTransformationMatrix(final Matrix matrix, final double[] vector) {
		this.init();

		synchronized(this.values)
		{
			// column-major order for OpenGl
			this.values[0] = (float) (matrix.get(0,0));
			this.values[1] = (float) (matrix.get(0,1));
			this.values[2] = (float) (matrix.get(0,2));
			this.values[3] = (0);
			this.values[4] = (float) (matrix.get(1,0));
			this.values[5] = (float) (matrix.get(1,1));
			this.values[6] = (float) (matrix.get(1,2));
			this.values[7] = (0);
			this.values[8] = (float) (matrix.get(2,0));
			this.values[9] = (float) (matrix.get(2,1));
			this.values[10] = (float) (matrix.get(2,2));
			this.values[11] = (0);
			this.values[12] = (float) (vector[0]);
			this.values[13] = (float) (vector[1]);
			this.values[14] = (float) (vector[2]);
			this.values[15] = (1);
		}
	}

	
	/**
	 * This function will change the contents of result, but will not change point.
	 */
	public void transformPoint(final double[] point, final double[] result) {
		result[0] = this.values[0] * point[0] + this.values[4] * point[1] + this.values[8] * point[2] + this.values[12];
		result[1] = this.values[1] * point[0] + this.values[5] * point[1] + this.values[9] * point[2] + this.values[13];
		result[2] = this.values[2] * point[0] + this.values[6] * point[1] + this.values[10] * point[2] + this.values[14];
	}



	/**
	 * The provided rotation matrix is:
	 * m00 m01 m02
	 * m10 m11 m12
	 * m20 m21 m22
	 * 
	 * And the provided translation vector is <v0 v1 v2>
	 */
	public void setTransformationMatrix(final float m00, final float m01, final float m02, final float m10, final float m11, final float m12, final float m20, final float m21, final float m22, final float v0, final float v1, final float v2) {
		this.init();

		synchronized(this.values) {
			// column-major order for OpenGl
			this.values[0] = (m00);
			this.values[1] = (m01);
			this.values[2] = (m02);
			this.values[3] = (0);
			this.values[4] = (m10);
			this.values[5] = (m11);
			this.values[6] = (m12);
			this.values[7] = (0);
			this.values[8] = (m20);
			this.values[9] = (m21);
			this.values[10] = (m22);
			this.values[11] = (0);
			this.values[12] = (v0);
			this.values[13] = (v1);
			this.values[14] = (v2);
			this.values[15] = (1);
		}
	}

	public void setIdentity() {
		//		 column-major order for OpenGl
		this.values[0] = (1);
		this.values[1] = (0);
		this.values[2] = (0);
		this.values[3] = (0);
		this.values[4] = (0);
		this.values[5] = (1);
		this.values[6] = (0);
		this.values[7] = (0);
		this.values[8] = (0);
		this.values[9] = (0);
		this.values[10] = (1);
		this.values[11] = (0);
		this.values[12] = (0);
		this.values[13] = (0);
		this.values[14] = (0);
		this.values[15] = (1);
	}



	/**
	 * Given a 4x4 array "matrix0", this function replaces it with the LU
	 * decomposition of a row-wise permutation of itself. The input parameters
	 * are "matrix0" and "dimen". The array "matrix0" is also an output
	 * parameter. The vector "row_perm[4]" is an output parameter that contains
	 * the row permutations resulting from partial pivoting. The output
	 * parameter "even_row_xchg" is 1 when the number of row exchanges is even,
	 * or -1 otherwise. Assumes data type is always double.
	 * 
	 * This function is similar to luDecomposition, except that it is tuned
	 * specifically for 4x4 matrices.
	 * 
	 * @return true if the matrix is nonsingular, or false otherwise.
	 */
	//
	// Reference: Press, Flannery, Teukolsky, Vetterling,
	// _Numerical_Recipes_in_C_, Cambridge University Press,
	// 1988, pp 40-45.
	//
	private boolean luDecomposition(final double[] matrix0, final int[] row_perm) {

		final double row_scale[] = new double[4];

		// Determine implicit scaling information by looping over rows
		{
			int i, j;
			int ptr, rs;
			double big, temp;

			ptr = 0;
			rs = 0;

			// For each row ...
			i = 4;
			while (i-- != 0) {
				big = 0.0;

				// For each column, find the largest element in the row
				j = 4;
				while (j-- != 0) {
					temp = matrix0[ptr++];
					temp = Math.abs(temp);
					if (temp > big) {
						big = temp;
					}
				}

				// Is the matrix singular?
				if (big == 0.0) {
					return false;
				}
				row_scale[rs++] = 1.0 / big;
			}
		}

		{
			int j;
			int mtx;

			mtx = 0;

			// For all columns, execute Crout's method
			for (j = 0; j < 4; j++) {
				int i, imax, k;
				int target, p1, p2;
				double sum, big, temp;

				// Determine elements of upper diagonal matrix U
				for (i = 0; i < j; i++) {
					target = mtx + (4 * i) + j;
					sum = matrix0[target];
					k = i;
					p1 = mtx + (4 * i);
					p2 = mtx + j;
					while (k-- != 0) {
						sum -= matrix0[p1] * matrix0[p2];
						p1++;
						p2 += 4;
					}
					matrix0[target] = sum;
				}

				// Search for largest pivot element and calculate
				// intermediate elements of lower diagonal matrix L.
				big = 0.0;
				imax = -1;
				for (i = j; i < 4; i++) {
					target = mtx + (4 * i) + j;
					sum = matrix0[target];
					k = j;
					p1 = mtx + (4 * i);
					p2 = mtx + j;
					while (k-- != 0) {
						sum -= matrix0[p1] * matrix0[p2];
						p1++;
						p2 += 4;
					}
					matrix0[target] = sum;

					// Is this the best pivot so far?
					if ((temp = row_scale[i] * Math.abs(sum)) >= big) {
						big = temp;
						imax = i;
					}
				}

				if (imax < 0) {
					new Exception().printStackTrace();
				}

				// Is a row exchange necessary?
				if (j != imax) {
					// Yes: exchange rows
					k = 4;
					p1 = mtx + (4 * imax);
					p2 = mtx + (4 * j);
					while (k-- != 0) {
						temp = matrix0[p1];
						matrix0[p1++] = matrix0[p2];
						matrix0[p2++] = temp;
					}

					// Record change in scale factor
					row_scale[imax] = row_scale[j];
				}

				// Record row permutation
				row_perm[j] = imax;

				// Is the matrix singular
				if (matrix0[(mtx + (4 * j) + j)] == 0.0) {
					return false;
				}

				// Divide elements of lower diagonal matrix L by pivot
				if (j != (4 - 1)) {
					temp = 1.0 / (matrix0[(mtx + (4 * j) + j)]);
					target = mtx + (4 * (j + 1)) + j;
					i = 3 - j;
					while (i-- != 0) {
						matrix0[target] *= temp;
						target += 4;
					}
				}
			}
		}

		return true;
	}

	/**
	 * From j3d.org's vecmath library. Solves a set of linear equations. The
	 * input parameters "matrix1", and "row_perm" come from luDecompostionD4x4
	 * and do not change here. The parameter "matrix2" is a set of column
	 * vectors assembled into a 4x4 matrix of floating-point values. The
	 * procedure takes each column of "matrix2" in turn and treats it as the
	 * right-hand side of the matrix equation Ax = LUx = b. The solution vector
	 * replaces the original column of the matrix.
	 * 
	 * If "matrix2" is the identity matrix, the procedure replaces its contents
	 * with the inverse of the matrix from which "matrix1" was originally
	 * derived.
	 */
	//
	// Reference: Press, Flannery, Teukolsky, Vetterling,
	// _Numerical_Recipes_in_C_, Cambridge University Press,
	// 1988, pp 44-45.
	//
	private void luBacksubstitution(final double[] matrix1, final int[] row_perm,
			final double[] matrix2) {

		int i, ii, ip, j, k;
		int rp;
		int cv, rv;

		// rp = row_perm;
		rp = 0;

		// For each column vector of matrix2 ...
		for (k = 0; k < 4; k++) {
			// cv = &(matrix2[0][k]);
			cv = k;
			ii = -1;

			// Forward substitution
			for (i = 0; i < 4; i++) {
				double sum;

				ip = row_perm[rp + i];
				sum = matrix2[cv + 4 * ip];
				matrix2[cv + 4 * ip] = matrix2[cv + 4 * i];
				if (ii >= 0) {
					// rv = &(matrix1[i][0]);
					rv = i * 4;
					for (j = ii; j <= i - 1; j++) {
						sum -= matrix1[rv + j] * matrix2[cv + 4 * j];
					}
				} else if (sum != 0.0) {
					ii = i;
				}
				matrix2[cv + 4 * i] = sum;
			}

			// Backsubstitution
			// rv = &(matrix1[3][0]);
			rv = 3 * 4;
			matrix2[cv + 4 * 3] /= matrix1[rv + 3];

			rv -= 4;
			matrix2[cv + 4 * 2] = (matrix2[cv + 4 * 2] - matrix1[rv + 3]
					* matrix2[cv + 4 * 3])
					/ matrix1[rv + 2];

			rv -= 4;
			matrix2[cv + 4 * 1] = (matrix2[cv + 4 * 1] - matrix1[rv + 2]
					* matrix2[cv + 4 * 2] - matrix1[rv + 3]
							* matrix2[cv + 4 * 3])
							/ matrix1[rv + 1];

			rv -= 4;
			matrix2[cv + 4 * 0] = (matrix2[cv + 4 * 0] - matrix1[rv + 1]
					* matrix2[cv + 4 * 1] - matrix1[rv + 2]
							* matrix2[cv + 4 * 2] - matrix1[rv + 3]
									* matrix2[cv + 4 * 3])
									/ matrix1[rv + 0];
		}
	}

	public ModelTransformationMatrix inverse3() {
		final ModelTransformationMatrix returned = new ModelTransformationMatrix();
		returned.init();

		final double temp[] = new double[16];
		final double result[] = new double[16];
		final int row_perm[] = new int[4];
		int i;
		// Copy source matrix to t1tmp
		temp[0] = this.values[0];
		temp[1] = this.values[4];
		temp[2] = this.values[8];
		temp[3] = this.values[12];

		temp[4] = this.values[1];
		temp[5] = this.values[5];
		temp[6] = this.values[9];
		temp[7] = this.values[13];

		temp[8] = this.values[2];
		temp[9] = this.values[6];
		temp[10] = this.values[10];
		temp[11] = this.values[14];

		temp[12] = this.values[3];
		temp[13] = this.values[7];
		temp[14] = this.values[11];
		temp[15] = this.values[15];

		// Calculate LU decomposition: Is the matrix singular?
		if (!luDecomposition(temp, row_perm)) {
			// Matrix has no inverse
			new Exception("Error: matrix has no inverse.").printStackTrace();
		}

		// Perform back substitution on the identity matrix
		for (i = 0; i < 16; i++) {
			result[i] = 0.0;
		}
		result[0] = 1.0;
		result[5] = 1.0;
		result[10] = 1.0;
		result[15] = 1.0;
		luBacksubstitution(temp, row_perm, result);

		returned.values[0] = (float) result[0];
		returned.values[4] = (float) result[1];
		returned.values[8] = (float) result[2];
		returned.values[12] = (float) result[3];

		returned.values[1] = (float) result[4];
		returned.values[5] = (float) result[5];
		returned.values[9] = (float) result[6];
		returned.values[13] = (float) result[7];

		returned.values[2] = (float) result[8];
		returned.values[6] = (float) result[9];
		returned.values[10] = (float) result[10];
		returned.values[14] = (float) result[11];

		returned.values[3] = (float) result[12];
		returned.values[7] = (float) result[13];
		returned.values[11] = (float) result[14];
		returned.values[15] = (float) result[15];

		return returned;
	}


	public void updateFullSymmetryDataWithInverseFractionalTransform(
			final ModelTransformationMatrix fractional,
			final ModelTransformationMatrix fractionalInverse) {
		ModelTransformationMatrix result = multiply4square_x_4square2(fractionalInverse,
				this);
		result.printMatrix("**fractional * symmetry**");
		result = multiply4square_x_4square2(result, fractional);
		this.values = result.values;

		// quick fix, to remove rounding errors.
		for(int i = 0; i < 16; i++) {
			final float val = this.values[i];
			if(val != 0 && val < 0.0000001 && val > -0.0000001) {
				this.values[i] = 0;
			}
		}

		this.printMatrix("**symmetry * inverse fractional**");
	}

	public static ModelTransformationMatrix multiply4square_x_4square2(
			final ModelTransformationMatrix leftMat, final ModelTransformationMatrix rightMat) {
		final ModelTransformationMatrix result = new ModelTransformationMatrix();
		result.init();

		float m00, m01, m02, m03, m10, m11, m12, m13, m20, m21, m22, m23, m30, m31, m32, m33; // vars
		// for
		// temp
		// result
		// matrix

		m00 = leftMat.values[0] * rightMat.values[0] + leftMat.values[4] * rightMat.values[1] + leftMat.values[8] * rightMat.values[2]
				+ leftMat.values[12] * rightMat.values[3];
		m01 = leftMat.values[0] * rightMat.values[4] + leftMat.values[4] * rightMat.values[5] + leftMat.values[8] * rightMat.values[6]
				+ leftMat.values[12] * rightMat.values[7];
		m02 = leftMat.values[0] * rightMat.values[8] + leftMat.values[4] * rightMat.values[9] + leftMat.values[8] * rightMat.values[10]
				+ leftMat.values[12] * rightMat.values[11];
		m03 = leftMat.values[0] * rightMat.values[12] + leftMat.values[4] * rightMat.values[13] + leftMat.values[8] * rightMat.values[14]
				+ leftMat.values[12] * rightMat.values[15];

		m10 = leftMat.values[1] * rightMat.values[0] + leftMat.values[5] * rightMat.values[1] + leftMat.values[9] * rightMat.values[2]
				+ leftMat.values[13] * rightMat.values[3];
		m11 = leftMat.values[1] * rightMat.values[4] + leftMat.values[5] * rightMat.values[5] + leftMat.values[9] * rightMat.values[6]
				+ leftMat.values[13] * rightMat.values[7];
		m12 = leftMat.values[1] * rightMat.values[8] + leftMat.values[5] * rightMat.values[9] + leftMat.values[9] * rightMat.values[10]
				+ leftMat.values[13] * rightMat.values[11];
		m13 = leftMat.values[1] * rightMat.values[12] + leftMat.values[5] * rightMat.values[13] + leftMat.values[9] * rightMat.values[14]
				+ leftMat.values[13] * rightMat.values[15];

		m20 = leftMat.values[2] * rightMat.values[0] + leftMat.values[6] * rightMat.values[1] + leftMat.values[10] * rightMat.values[2]
				+ leftMat.values[14] * rightMat.values[3];
		m21 = leftMat.values[2] * rightMat.values[4] + leftMat.values[6] * rightMat.values[5] + leftMat.values[10] * rightMat.values[6]
				+ leftMat.values[14] * rightMat.values[7];
		m22 = leftMat.values[2] * rightMat.values[8] + leftMat.values[6] * rightMat.values[9] + leftMat.values[10] * rightMat.values[10]
				+ leftMat.values[14] * rightMat.values[11];
		m23 = leftMat.values[2] * rightMat.values[12] + leftMat.values[6] * rightMat.values[13] + leftMat.values[10] * rightMat.values[14]
				+ leftMat.values[14] * rightMat.values[15];

		m30 = leftMat.values[3] * rightMat.values[0] + leftMat.values[7] * rightMat.values[1] + leftMat.values[11] * rightMat.values[2]
				+ leftMat.values[15] * rightMat.values[3];
		m31 = leftMat.values[3] * rightMat.values[4] + leftMat.values[7] * rightMat.values[5] + leftMat.values[11] * rightMat.values[6]
				+ leftMat.values[15] * rightMat.values[7];
		m32 = leftMat.values[3] * rightMat.values[8] + leftMat.values[7] * rightMat.values[9] + leftMat.values[11] * rightMat.values[10]
				+ leftMat.values[15] * rightMat.values[11];
		m33 = leftMat.values[3] * rightMat.values[12] + leftMat.values[7] * rightMat.values[13] + leftMat.values[11] * rightMat.values[14]
				+ leftMat.values[15] * rightMat.values[15];

		result.values[0] = m00;
		result.values[4] = m01;
		result.values[8] = m02;
		result.values[12] = m03;
		result.values[1] = m10;
		result.values[5] = m11;
		result.values[9] = m12;
		result.values[13] = m13;
		result.values[2] = m20;
		result.values[6] = m21;
		result.values[10] = m22;
		result.values[14] = m23;
		result.values[3] = m30;
		result.values[7] = m31;
		result.values[11] = m32;
		result.values[15] = m33;

		return result;
	}

	// create the transformation matrix using a full symmetry operation string.
	public static final Pattern spaces = Pattern.compile("\\s+");

	public static final Pattern commaSpaces = Pattern.compile("\\s*,\\s*");

	public static final Pattern slash = Pattern.compile("/");

	public void setFullSymmetryOperation(final String fullSymmetryOperation_) {
		// if(this.symmetryShorthand.equals("7_555")) {
		// fullSymmetryOperation_ = "y,x,1/2-z";
		// }

		if (fullSymmetryOperation_ == null) {
			return;
		}

		final String fullSymmetryOperation = ModelTransformationMatrix.spaces.matcher(fullSymmetryOperation_)
				.replaceAll("");
		final String[] xyzRawArray = ModelTransformationMatrix.commaSpaces.split(fullSymmetryOperation);

		if (xyzRawArray == null || xyzRawArray.length != 3) {
			return;
		}

		this.init();

		// swap and/or invert the axes...
		// this sets up the values in the upper left 3x3 matrix.
		for (int i = 0; i < xyzRawArray.length; i++) {
			final String xyzRaw = xyzRawArray[i];

			// zero out the values first.
			this.values[i] = 0f;
			this.values[i + 4] = 0f;
			this.values[i + 8] = 0f;

			int coordIndexX = xyzRaw.indexOf('x');
			if (coordIndexX >= 0) {
				if (coordIndexX != 0 && xyzRaw.charAt(coordIndexX - 1) == '-') {
					this.values[i] = -1f;
				} else {
					this.values[i] = 1f;
				}
			}

			int coordIndexY = xyzRaw.indexOf('y');

			if (coordIndexY >= 0) {
				if (coordIndexY != 0 && xyzRaw.charAt(coordIndexY - 1) == '-') {
					this.values[i + 4] = -1f;
				} else {
					this.values[i + 4] = 1f;
				}
			}

			int coordIndexZ = xyzRaw.indexOf('z');

			if (coordIndexZ >= 0) {
				if (coordIndexZ != 0
						&& xyzRaw.charAt(coordIndexZ - 1) == '-') {
					this.values[i + 8] = -1f;
				} else {
					this.values[i + 8] = 1f;
				}
			}

			int coordIndex = -1;
			if(coordIndexX >=0) {
				coordIndex = coordIndexX;
			}
			if(coordIndex >= 0) {
				if(coordIndexY >= 0) {
					coordIndex = Math.min(coordIndex, coordIndexY);
				}
			} else {
				coordIndex = coordIndexY;
			}
			if(coordIndex >= 0) {
				if(coordIndexZ >= 0) {
					coordIndex = Math.min(coordIndex, coordIndexZ);
				}
			} else {
				coordIndex = coordIndexZ;
			}
			if(coordIndex < 0) {
				// error
				this.values = null;
				return;
			}


			// set the translation coordinate in case the translation is zero.
			final int translationIndex = 12 + i;
			// if(isX) {
			// translationIndex = 12;
			// } else if(isY) {
			// translationIndex = 13;
			// } else { // isZ
			// translationIndex = 14;
			// }
			this.values[translationIndex] = 0f;

			if (coordIndex != 0) {
				// find out where the fraction (if any) ends...
				final char tmp = xyzRaw.charAt(coordIndex - 1);
				if (tmp == '-' || tmp == '+') {
					coordIndex--;
				}

				if (coordIndex != 0) {
					int fractionStartIndex = 0;
					boolean isNegated = false;
					if (xyzRaw.charAt(0) == '-') {
						isNegated = true;
						fractionStartIndex = 1;
					}

					// get the fraction and convert it to a real
					final String fractionSt = xyzRaw.substring(fractionStartIndex,
							coordIndex);
					final String[] fractionParts = ModelTransformationMatrix.slash.split(fractionSt);
					if (fractionParts == null || fractionParts.length == 0
							|| fractionParts.length > 2) { // error
						return;
					}
					float fraction = Float.parseFloat(fractionParts[0]);
					if (fractionParts.length == 2) {
						fraction /= Float.parseFloat(fractionParts[1]);
					}
					if (isNegated) {
						fraction = -fraction;
					}
					// fraction *= 193.800f;

					// do the translation...
					// this sets up the bottom left 3 element horizontal vector.
					this.values[translationIndex] = fraction;
				}
			}
		}

		// finish by setting up the right column of the matrix.
		this.values[3] = 0f;
		this.values[7] = 0f;
		this.values[11] = 0f;
		this.values[15] = 1f;

		this.printMatrix(fullSymmetryOperation);
	}

	public void printMatrix(final String fullSymmetryOperation)
	{
		/* ** DEBUGGING - printMatrix
		System.err.println("Generated transformation matrix from full symmetry "
						+ fullSymmetryOperation + ", chain id "
						+ this.ndbChainId + ": ");

		for (int row = 0; row < 4; row++)
		{
			for (int column = 0; column < 4; column++)
				System.err.print(this.modelTransform[column * 4 + row] + "\t");

			System.err.println();
		}
		System.err.println();
		 * **/
	}

	public void init() {
		if (this.values == null)
			values = new float[16];
		setIdentity();
	}



	public Matrix getMatrix(){



		Matrix m = new Matrix(3,3);
		m.set(0,0,values[0]);
		m.set(1,0,values[1]);
		m.set(2,0,values[2]);
		m.set(0,1,values[4]);
		m.set(1,1,values[5]);
		m.set(2,1,values[6]);
		m.set(0,2,values[8]);
		m.set(1,2,values[9]);
		m.set(2,2,values[10]);

		return m;
	}

	public void setMatrix(Matrix m){
		values[0] = (float)m.get(0,0);
		values[1] = (float)m.get(0,1);
		values[2] = (float)m.get(0,2);
		values[4] = (float)m.get(1,0);
		values[5] = (float)m.get(1,1);
		values[6] = (float)m.get(1,2);
		values[8] = (float)m.get(2,0);
		values[9] = (float)m.get(2,1);
		values[10] = (float)m.get(2,2);
	}
	public double[] getVector(){
		double[] v = new double[3];
		v[0] = values[12];
		v[1] = values[13];
		v[2] = values[14];

		return v;

	}
	public void setVector(double[] v){

		values[12] = (float) v[0];
		values[13] = (float) v[1];
		values[14] = (float) v[2];

	}

	@Override
	public String toString() {

		Matrix m = getMatrix();
		double[] v = getVector();
		return "ModelTransformationMatrix [id=" + id + ", ndbChainId="
		+ ndbChainId + ", symmetryShorthand=" + symmetryShorthand
		+ ", code=" + code + ", values=" + m.toString() + " " +  Arrays.toString(v)
		+ "]";
	}

	
	

	public String toXML() throws IOException{

		StringWriter sw = new StringWriter();
		PrintWriter writer = new PrintWriter(sw);

		PrettyXMLWriter xml = new PrettyXMLWriter(new PrintWriter(writer));
		
		toXML(xml);
		
		xml.close();
		writer.close();
		sw.close();
		return sw.toString();
	}
	
	public void toXML(PrettyXMLWriter xml) throws IOException{
		xml.openTag("transformation");
		xml.attribute("index",id);
		//xml.attribute("chainId", ndbChainId);

		String shorthand = symmetryShorthand;
		if ( shorthand != null)
			xml.attribute("symmetryShorthand", shorthand);


		if ( code != null) 
			xml.attribute("code",code);

		xml.openTag("matrix");
		Matrix m = getMatrix();
		
			for ( int i = 0 ; i<3 ; i++){
				for ( int j = 0 ; j<3 ;j++){	
				xml.attribute("m" +  (i+1) + (j+1), String.format("%.8f",m.get(i,j)));
			}
		}
		xml.closeTag("matrix");

		double[] shift = getVector();
		xml.openTag("shift");
		for ( int i = 0 ; i<3 ; i++) {
			xml.attribute("v"+(i+1),String.format("%.8f", shift[i]));
		}
		xml.closeTag("shift");

		xml.closeTag("transformation");
		
	}

	public static ModelTransformationMatrix fromXML(String xml) 
			throws SAXException, 
			IOException, 
			ParserConfigurationException{

	
		List<ModelTransformationMatrix> transformations = fromMultiXML(xml);
		
		if ( transformations.size() > 0)
			return transformations.get(0);
		
		else
			return null;
	}
	
	
	public static List<ModelTransformationMatrix> fromMultiXML(String xml) throws ParserConfigurationException, SAXException, IOException{
		
	
		List<ModelTransformationMatrix> transformations = new ArrayList<ModelTransformationMatrix>();
		
		// read the XML of a string and returns a ModelTransformationmatrix
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = factory.newDocumentBuilder();

		InputSource inStream = new InputSource();
		inStream.setCharacterStream(new StringReader(xml));
		Document doc = db.parse(inStream);

		// normalize text representation
		doc.getDocumentElement().normalize();


		//Element rootElement = doc.getDocumentElement();

		NodeList listOfTransforms = doc.getElementsByTagName("transformation");

		for(int pos=0; pos<listOfTransforms.getLength() ; pos++) {
			Node rootElement       = listOfTransforms.item(pos);

			ModelTransformationMatrix max = new ModelTransformationMatrix();
			
			max.id = getAttribute(rootElement,"index");				
			max.ndbChainId = getAttribute(rootElement,"chainId");

			max.code = getAttribute(rootElement, "code");
			max.symmetryShorthand = getAttribute(rootElement, "symmetryShorthand");



			NodeList listOfChildren = rootElement.getChildNodes();


			for(int i=0; i<listOfChildren.getLength() ; i++)
			{
				// and now the matrix ...
				Node block       = listOfChildren.item(i);

				// we only look at blocks.
				if ( block.getNodeName().equals("matrix")) 
					max.setMatrix(getMatrixFromXML(block));

				if ( block.getNodeName().equals("shift")) 
					max.setVector(getVectorFromXML(block));

			}

			transformations.add(max);
		}

		return transformations;
	}

	private static double[] getVectorFromXML(Node block) {
		double[] d = new double[3];
		for ( int i = 0 ; i<3 ; i++){
			d[i] = Float.parseFloat(getAttribute(block, "v" + (i+1) ));
		}
		return d;
	}

	private static Matrix getMatrixFromXML(Node block) {
		Matrix m  = new Matrix(3,3);
		for ( int i = 0 ; i<3 ; i++){
			for ( int j = 0 ; j<3 ;j++){
				String val = getAttribute(block, "m" + (j+1)+(i+1));			
				m.set(j,i, Float.parseFloat(val));
			}
		}
		return m;
	}

	private static String getAttribute(Node node, String attr){
		if( ! node.hasAttributes()) 
			return null;

		NamedNodeMap atts = node.getAttributes();

		if ( atts == null)
			return null;

		Node att = atts.getNamedItem(attr);
		if ( att == null)
			return null;

		String value = att.getTextContent();

		return value;

	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getNdbChainId() {
		return ndbChainId;
	}

	public void setNdbChainId(String ndbChainId) {
		this.ndbChainId = ndbChainId;
	}
	

	public String getSymmetryShorthand() {
		return symmetryShorthand;
	}

	public void setSymmetryShorthand(String symmetryShorthand) {
		this.symmetryShorthand = symmetryShorthand;
	}
	
	

	public String getCode() {
		return code;
	}

	public void setCode(String code) {
		this.code = code;
	}

	public ModelTransformationMatrix clone(){
		ModelTransformationMatrix m = new ModelTransformationMatrix();
		
		m.setMatrix(getMatrix());
		m.setVector(getVector());
		m.setId(getId());
		m.setNdbChainId(getNdbChainId());
		m.setSymmetryShorthand(getSymmetryShorthand());
		m.setCode(getCode());
		return m;
		
	}

	
		
	
	
}