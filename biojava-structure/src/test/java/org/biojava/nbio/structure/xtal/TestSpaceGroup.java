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
package org.biojava.nbio.structure.xtal;

import org.biojava.nbio.structure.jama.Matrix;
import org.junit.Assert;
import org.junit.Test;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import java.util.Collection;

/**
 * Testing of space group symop.lib parsing and for crystal operator calculation correctness
 *
 * @author duarte_j
 *
 */
public class TestSpaceGroup {

	private static final double DELTA = 0.000001;
	// use true for printing operators of all space groups (including non-enantiomorphics), or false to print only enantiomorphics
	private static boolean PRINT_OPERATORS_FROM_ALL_SGS = false;

	// print information per operator
	private static boolean VERBOSE = false;

	@Test
	public void testTransfConversion() {
		Collection<SpaceGroup> allSGs = SymoplibParser.getAllSpaceGroups().values();

		int countEn = 0;
		int countNonEn = 0;
		int countSpecial = 0;
		for (SpaceGroup spaceGroup:allSGs) {

			if (spaceGroup.isEnantiomorphic() && spaceGroup.getId()<1000) {
				countEn++;
			}
			if (!spaceGroup.isEnantiomorphic() && spaceGroup.getId()<1000) {
				countNonEn++;
			}
			if (spaceGroup.getId()>1000) {
				countSpecial++;
			}

			if (VERBOSE)
				System.out.println(spaceGroup.getId()+"  "+spaceGroup.getShortSymbol()+" -- "+spaceGroup.getMultiplicity()+" "+spaceGroup.getPrimitiveMultiplicity());


			// cell translations must be the same in each subgroup of operators (applies to I, C, F and H SGs)
			if (spaceGroup.getId()<1000) {
				int fold = spaceGroup.getMultiplicity()/spaceGroup.getPrimitiveMultiplicity();
				for (int n=1;n<fold;n++) {
					for (int j=0;j<spaceGroup.getPrimitiveMultiplicity();j++) {
						Matrix4d t = spaceGroup.getTransformation(n*spaceGroup.getPrimitiveMultiplicity()+j);
						Matrix4d tPrimitive = spaceGroup.getTransformation(j);
						Vector3d cellTransl = new Vector3d(t.m03,t.m13,t.m23);
						Vector3d primitive = new Vector3d(tPrimitive.m03,tPrimitive.m13,tPrimitive.m23);
						cellTransl.sub(spaceGroup.getCellTranslation(n));

						double diffx = primitive.x - cellTransl.x;
						double diffy = primitive.y - cellTransl.y;
						double diffz = primitive.z - cellTransl.z;

						if (!(
								(Math.abs(diffx)<0.000001 &&
								 Math.abs(diffy)<0.000001 &&
								 Math.abs(diffz)<0.000001)
								 )) {

							//System.out.printf("DIFF NOT 0: \n" +
							//		" primitive (%5.2f,%5.2f,%5.2f), multiple (%5.2f,%5.2f,%5.2f), difference (%5.2f,%5.2f,%5.2f)\n",
							//		primitive.x, primitive.y, primitive.z,
							//		cellTransl.x ,cellTransl.y, cellTransl.z,
							//		diffx, diffy, diffz);

							// they are positive and not larger than 1
							Assert.assertFalse(diffx>1 && diffx<0);
							Assert.assertFalse(diffy>1 && diffy<0);
							Assert.assertFalse(diffz>1 && diffz<0);

							// they are integer, i.e. 0 or 1
							Assert.assertTrue(diffx == (int)diffx);
							Assert.assertTrue(diffy == (int)diffy);
							Assert.assertTrue(diffz == (int)diffz);
						}

					}
				}

			}

			for (int i=0;i<spaceGroup.getNumOperators();i++){
				CrystalCell unitCell = spaceGroup.getBravLattice().getExampleUnitCell();
				Matrix4d m = spaceGroup.getTransformation(i);
				Matrix4d mT = unitCell.transfToOrthonormal(m);

				// as stated in PDB documentation for SCALE matrix (our MTranspose matrix) the inverse determinant should match the cell volume
				Assert.assertEquals(unitCell.getVolume(),1.0/unitCell.getMTranspose().determinant(),DELTA);
				// and checking that our method to check scale matrix works as expected
				Matrix4d scaleMat = new Matrix4d(unitCell.getMTranspose(),new Vector3d(0,0,0),1);
				Assert.assertTrue(unitCell.checkScaleMatrixConsistency(scaleMat));

				// traces before and after transformation must coincide
				Assert.assertEquals(m.m00+m.m11+m.m22+m.m33, mT.m00+mT.m11+mT.m22+mT.m33, DELTA);

				Matrix3d rot = new Matrix3d(mT.m00,mT.m01,mT.m02,mT.m10,mT.m11,mT.m12,mT.m20,mT.m21,mT.m22);

				// determinant is either 1 or -1 (for improper rotations i.e. mirrors)
				Assert.assertTrue(Math.abs(rot.determinant()-1)<DELTA || Math.abs(rot.determinant()+1)<DELTA);

				CrystalTransform ct = new CrystalTransform(spaceGroup, i);

				if (spaceGroup.isEnantiomorphic() && spaceGroup.getId()<1000) {

					// determinant must be 1
					Assert.assertEquals(1.0,rot.determinant(), DELTA);

					// at least 1 eigenvalue must be 1 (there's one direction that remains unchanged under rotation)
					double[][] ar = {{mT.m00,mT.m01,mT.m02},{mT.m10,mT.m11,mT.m12},{mT.m20,mT.m21,mT.m22}};
					Matrix mat = new Matrix(ar);
					double[] eigenv = mat.svd().getSingularValues();
					Assert.assertTrue(eq(eigenv[0],1.0) || eq(eigenv[1],1.0) || eq(eigenv[2],1.0));

					// transpose must be equals to inverse
					Matrix3d rotTransp = new Matrix3d();
					Matrix3d rotInv = new Matrix3d();
					rotTransp.transpose(rot);
					rotInv.invert(rot);

					assertMatrixEquals(rotTransp,rotInv);

					int foldType = spaceGroup.getAxisFoldType(i);

					Assert.assertTrue(foldType==1 || foldType==2 || foldType==3 || foldType==4 || foldType==6);

					if (!ct.isPureTranslation()) {
						Assert.assertTrue(ct.isIdentity() || ct.isRotation());
					}
					if (!ct.isRotation()) {
						Assert.assertTrue(ct.isIdentity() || ct.isPureTranslation());
					}
					if (!ct.isPureTranslation() && ct.isFractionalTranslation()) {
						Assert.assertTrue(foldType!=1);
						Assert.assertTrue(ct.isRotation());
					}

				}

				// i=0 is identity
				if (i==0) {
					Assert.assertTrue(ct.isIdentity());
					Assert.assertFalse(ct.isPureCrystalTranslation());
					Assert.assertTrue(ct.getTransformType()==TransformType.AU);
				}

				Assert.assertFalse(ct.isPureCrystalTranslation());
				Assert.assertFalse(ct.getTransformType()==TransformType.XTALTRANSL);


				// rotation axes and type

				int foldType = spaceGroup.getAxisFoldType(i);
				Matrix3d W = new Matrix3d(m.m00,m.m01,m.m02,m.m10,m.m11,m.m12,m.m20,m.m21,m.m22);
				AxisAngle4d axisAngle = spaceGroup.getRotAxisAngle(i);
				if (W.determinant()>0) { // i.e. 1
					switch (foldType) {
					case 1:
						Assert.assertEquals(0, axisAngle.angle, DELTA);
						break;
					case 2:
						Assert.assertEquals(Math.PI, axisAngle.angle, DELTA);
						break;
					case 3:
						Assert.assertEquals(2.0*Math.PI/3.0, axisAngle.angle, DELTA);
						break;
					case 4:
						Assert.assertEquals(Math.PI/2.0, axisAngle.angle, DELTA);
						break;
					case 6:
						Assert.assertEquals(Math.PI/3.0, axisAngle.angle, DELTA);
						break;
					}

				} else { // i.e. determinant -1
					switch (foldType) {
					case -1:
						Assert.assertEquals(0,axisAngle.angle,DELTA);
						// no glide planes
						Assert.assertTrue(ct.getTranslScrewComponent().epsilonEquals(new Vector3d(0,0,0), DELTA));
						break;
					case -2:
						Assert.assertEquals(0,axisAngle.angle,DELTA);
						// glide planes can happen
						break;
					case -3:
						Assert.assertEquals(0,axisAngle.angle,DELTA);
						// no glide planes
						Assert.assertTrue(ct.getTranslScrewComponent().epsilonEquals(new Vector3d(0,0,0), DELTA));
						break;
					case -4:
						Assert.assertEquals(0,axisAngle.angle,DELTA);
						// no glide planes
						Assert.assertTrue(ct.getTranslScrewComponent().epsilonEquals(new Vector3d(0,0,0), DELTA));
						break;
					case -6:
						Assert.assertEquals(0,axisAngle.angle,DELTA);
						// no glide planes
						Assert.assertTrue(ct.getTranslScrewComponent().epsilonEquals(new Vector3d(0,0,0), DELTA));
						break;

					}

				}

				Vector3d axis = new Vector3d(axisAngle.x, axisAngle.y, axisAngle.z);
				Vector3d translScrewComponent = ct.getTranslScrewComponent();

				// if both non-0, then both have to be on the same direction (with perhaps different sense)
				if (!axis.epsilonEquals(new Vector3d(0,0,0), DELTA) &&
					!translScrewComponent.epsilonEquals(new Vector3d(0,0,0), DELTA)) {

					Assert.assertTrue ((eq(axis.angle(translScrewComponent),0.0)) ||
										eq(axis.angle(translScrewComponent),Math.PI) );
				}

				// we only print the enantiomorphic ones, otherwise there's too much output

				if (spaceGroup.isEnantiomorphic() || PRINT_OPERATORS_FROM_ALL_SGS)
					if (VERBOSE)
						System.out.print(
							String.format("%2d %2d",i,1+i-spaceGroup.getPrimitiveMultiplicity()*(i/spaceGroup.getPrimitiveMultiplicity()))+" "+
							String.format("%20s",spaceGroup.getTransfAlgebraic(i))+" "+
							String.format("(%5.2f,%5.2f,%5.2f)",axis.x,axis.y,axis.z)+" "+
							String.format("%2d",foldType)+" " +
							String.format("(%5.2f,%5.2f,%5.2f)",translScrewComponent.x,translScrewComponent.y,translScrewComponent.z));


				if (ct.getTransformType().isScrew()) { // tests for both screw or glide planes (i.e. a non-zero scre transl component)

					if (spaceGroup.isEnantiomorphic()  || PRINT_OPERATORS_FROM_ALL_SGS) {
						if (VERBOSE)
							System.out.print(" -- SCREW AXIS");
						if (foldType==-2) {
							if (VERBOSE)
								System.out.print(" (GLIDE)");
						}
					}


					Assert.assertTrue(
							ct.getTransformType()==TransformType.GLIDE ||
							ct.getTransformType()==TransformType.TWOFOLDSCREW ||
							ct.getTransformType()==TransformType.THREEFOLDSCREW ||
							ct.getTransformType()==TransformType.FOURFOLDSCREW ||
							ct.getTransformType()==TransformType.SIXFOLDSCREW);

					Assert.assertFalse(ct.isPureTranslation());
				}

				if (ct.isPureTranslation()) {
					Assert.assertEquals(1,foldType);

					if (spaceGroup.isEnantiomorphic()  || PRINT_OPERATORS_FROM_ALL_SGS) {
						if (VERBOSE)
							System.out.print(" -- FRACTIONAL TRANSLATION");
					}

					Assert.assertTrue(ct.isFractionalTranslation());
					Assert.assertTrue(ct.getTransformType()==TransformType.CELLTRANSL);
				}
				if (ct.isRotation() && !ct.isFractionalTranslation()) {
					Assert.assertTrue(ct.getTransformType()==TransformType.TWOFOLD ||
							ct.getTransformType()==TransformType.THREEFOLD ||
							ct.getTransformType()==TransformType.FOURFOLD ||
							ct.getTransformType()==TransformType.SIXFOLD);
					Assert.assertTrue(!ct.getTransformType().isScrew());
				}

				if (spaceGroup.isEnantiomorphic() || PRINT_OPERATORS_FROM_ALL_SGS) {
					if (VERBOSE) {
						System.out.print(" -- "+ct.getTransformType().getShortName());

						System.out.println();
					}

				}

			}


		}


		Assert.assertEquals(266, allSGs.size()); // the total count must be 266
		Assert.assertEquals(65, countEn);        // enantiomorphic groups (protein crystallography groups)
		Assert.assertEquals(165, countNonEn);    // i.e. 266-65-36
		Assert.assertEquals(36, countSpecial);   // the rest of the groups present un symop.lib (sometimes used in PDB)
	}

	private static void assertMatrixEquals(Matrix3d m1, Matrix3d m2) {
		for (int i=0;i<3;i++) {
			for (int j=0;j<3;j++) {
				 Assert.assertEquals(m1.getElement(i, j),m2.getElement(i, j),DELTA);
			}
		}
	}

	private static boolean eq(double d1, double d2) {
		return (Math.abs(d1-d2)<DELTA);
	}

	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		TestSpaceGroup test = new TestSpaceGroup();
		test.testTransfConversion();
	}

}
