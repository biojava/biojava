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

package org.biojava.nbio.structure.symmetry.geometry;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;


/**
 *
 * @author Peter
 */
public final class SuperPositionQCP {
    double evecprec = 1d-6;
    double evalprec = 1d-11;
    
    private Point3d[] x = null;
    private Point3d[] y = null;
    
    private double[] weight = null;
    
    private Point3d[] xref = null;
    private Point3d[] yref = null;
    private Point3d xtrans;
    private Point3d ytrans;
    
    private double e0;
    private Matrix3d rotmat = new Matrix3d();
    private Matrix4d transformation = new Matrix4d();
    private double rmsd = 0;
    private double Sxy, Sxz, Syx, Syz, Szx, Szy;
    private double SxxpSyy, Szz, mxEigenV, SyzmSzy,SxzmSzx, SxymSyx;
    private double SxxmSyy, SxypSyx, SxzpSzx;
    private double Syy, Sxx, SyzpSzy;
    private boolean rmsdCalculated = false;
    private boolean transformationCalculated = false;
    private boolean centered = false;
    
    
    public void set(Point3d[] x, Point3d[] y) {
    	this.x = x;
    	this.y = y;
        rmsdCalculated = false;
        transformationCalculated = false;
    }
    
    public void set(Point3d[] x, Point3d[] y, double[] weight) {
    	this.x = x;
    	this.y = y;
    	this.weight = weight;
        rmsdCalculated = false;
        transformationCalculated = false;
    }
    
    public void setCentered(boolean centered) {
    	this.centered = centered;
    }
    
    public double getRmsd() {
    	if (! rmsdCalculated) {
    		calcRmsd(x, y);
    	}
    	return rmsd;
    }
    
    public Matrix4d getTransformationMatrix() {
    	getRotationMatrix();
    	if (! centered) {
    		calcTransformation();
    	} else {
    		transformation.set(rotmat);
    	}
    	return transformation; 	
    }
    
    public Matrix3d getRotationMatrix() {
    	getRmsd();
    	if (! transformationCalculated) {
    		calcRotationMatrix();
    	}
    	return rotmat; 	
    }
    
    public Point3d[] getTransformedCoordinates() {
    	SuperPosition.transform(transformation, x);
    	return x;
    }
    
    /**
     * this requires the coordinates to be precentered
     * @param x
     * @param y
     */
    private void calcRmsd(Point3d[] x, Point3d[] y) {
    	if (centered) {
    		innerProduct(y, x);
    	} else {
    		// translate to origin
    		xref = SuperPosition.clonePoint3dArray(x);
    		xtrans = SuperPosition.centroid(xref);
    		//  	System.out.println("x centroid: " + xtrans);
    		xtrans.negate();
    		SuperPosition.translate(xtrans, xref);

    		yref = SuperPosition.clonePoint3dArray(y);
    		ytrans = SuperPosition.centroid(yref);
    		// 	System.out.println("y centroid: " + ytrans);
    		ytrans.negate();
    		SuperPosition.translate(ytrans, yref); 
    		innerProduct(yref, xref);
    	}
    	calcRmsd(x.length);
    }
  
    /* Superposition coords2 onto coords1 -- in other words, coords2 is rotated, coords1 is held fixed */
    private void calcTransformation() {
     //   transformation.set(rotmat,new Vector3d(0,0,0), 1);
        transformation.set(rotmat);
  //      long t2 = System.nanoTime();
  //      System.out.println("create transformation: " + (t2-t1));
 //       System.out.println("m3d -> m4d");
 //       System.out.println(transformation);

        // combine with x -> origin translation
        Matrix4d trans = new Matrix4d();
        trans.setIdentity();
        trans.setTranslation(new Vector3d(xtrans));
        transformation.mul(transformation, trans);
//        System.out.println("setting xtrans");
//        System.out.println(transformation);
//
//        // combine with origin -> y translation
        ytrans.negate();  
        Matrix4d transInverse = new Matrix4d(); 
        transInverse.setIdentity();     
        transInverse.setTranslation(new Vector3d(ytrans));
        transformation.mul(transInverse, transformation);
//        System.out.println("setting ytrans");
//        System.out.println(transformation);
    }
    
    /** 
     * http://theobald.brandeis.edu/qcp/qcprot.c
     * @param A
     * @param coords1
     * @param coords2
     * @return
     */
    private void innerProduct(Point3d[] coords1, Point3d[] coords2) {
    	double          x1, x2, y1, y2, z1, z2;
        double          g1 = 0.0, g2 = 0.0;
        
        Sxx = 0;
        Sxy = 0;
        Sxz = 0;
        Syx = 0;
        Syy = 0;
        Syz = 0;
        Szx = 0;
        Szy = 0;
        Szz = 0; 

        if (weight != null)
        {
            for (int i = 0; i < coords1.length; i++)
            {
                x1 = weight[i] * coords1[i].x;
                y1 = weight[i] * coords1[i].y;
                z1 = weight[i] * coords1[i].z;

                g1 += x1 * coords1[i].x + y1 * coords1[i].y + z1 * coords1[i].z;

                x2 = coords2[i].x;
                y2 = coords2[i].y;
                z2 = coords2[i].z;

                g2 += weight[i] * (x2 * x2 + y2 * y2 + z2 * z2);

                Sxx +=  (x1 * x2);
                Sxy +=  (x1 * y2);
                Sxz +=  (x1 * z2);

                Syx +=  (y1 * x2);
                Syy +=  (y1 * y2);
                Syz +=  (y1 * z2);

                Szx +=  (z1 * x2);
                Szy +=  (z1 * y2);
                Szz +=  (z1 * z2);   
            }
        }
        else
        {
            for (int i = 0; i < coords1.length; i++)
            { 
                g1 += coords1[i].x * coords1[i].x + coords1[i].y * coords1[i].y + coords1[i].z * coords1[i].z;
                g2 += coords2[i].x * coords2[i].x + coords2[i].y * coords2[i].y + coords2[i].z * coords2[i].z;

                Sxx +=  coords1[i].x * coords2[i].x;
                Sxy +=  coords1[i].x * coords2[i].y;
                Sxz +=  coords1[i].x * coords2[i].z;

                Syx +=  coords1[i].y * coords2[i].x;
                Syy +=  coords1[i].y * coords2[i].y;
                Syz +=  coords1[i].y * coords2[i].z;

                Szx +=  coords1[i].z * coords2[i].x;
                Szy +=  coords1[i].z * coords2[i].y;
                Szz +=  coords1[i].z * coords2[i].z;  
            }
        }

        e0 = (g1 + g2) * 0.5;
    }

    private int calcRmsd(int len)
    {     
        double Sxx2 = Sxx * Sxx;
        double Syy2 = Syy * Syy;
        double Szz2 = Szz * Szz;

        double Sxy2 = Sxy * Sxy;
        double Syz2 = Syz * Syz;
        double Sxz2 = Sxz * Sxz;

        double Syx2 = Syx * Syx;
        double Szy2 = Szy * Szy;
        double Szx2 = Szx * Szx;

        double SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
        double Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

        double c2 = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
        double c1 = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

        SxzpSzx = Sxz + Szx;
        SyzpSzy = Syz + Szy;
        SxypSyx = Sxy + Syx;
        SyzmSzy = Syz - Szy;
        SxzmSzx = Sxz - Szx;
        SxymSyx = Sxy - Syx;
        SxxpSyy = Sxx + Syy;
        SxxmSyy = Sxx - Syy;
        
        double Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

        double c0 = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
             + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
             + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
             + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
             + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

        mxEigenV = e0;      
 
        int i;
        for (i = 0; i < 50; ++i)
        {
            double oldg = mxEigenV;
            double x2 = mxEigenV*mxEigenV;
            double b = (x2 + c2)*mxEigenV;
            double a = b + c1;
            double delta = ((a*mxEigenV + c0)/(2.0*x2*mxEigenV + b + a));
            mxEigenV -= delta;

            if (Math.abs(mxEigenV - oldg) < Math.abs(evalprec*mxEigenV))
                break;
        }
        

        if (i == 50) 
           System.err.println("More than %d iterations needed!" + i);

        /* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
        rmsd = Math.sqrt(Math.abs(2.0 * (e0 - mxEigenV)/len));
       
        return 1;
    }
    
    private int calcRotationMatrix() {
        double a11 = SxxpSyy + Szz-mxEigenV;
        double a12 = SyzmSzy; 
        double a13 = - SxzmSzx; 
        double a14 = SxymSyx;
        double a21 = SyzmSzy; 
        double a22 = SxxmSyy - Szz-mxEigenV; 
        double a23 = SxypSyx; 
        double a24= SxzpSzx;
        double a31 = a13; 
        double a32 = a23; 
        double a33 = Syy-Sxx-Szz - mxEigenV; 
        double a34 = SyzpSzy;
        double a41 = a14; 
        double a42 = a24; 
        double a43 = a34; 
        double a44 = Szz - SxxpSyy - mxEigenV;
        double a3344_4334 = a33 * a44 - a43 * a34; 
        double a3244_4234 = a32 * a44-a42*a34;
        double a3243_4233 = a32 * a43 - a42 * a33; 
        double a3143_4133 = a31 * a43-a41*a33;
        double a3144_4134 = a31 * a44 - a41 * a34; 
        double a3142_4132 = a31 * a42-a41*a32;
        double q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
        double q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
        double q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
        double q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

        double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

    /* The following code tries to calculate another column in the adjoint matrix when the norm of the 
       current column is too small.
       Usually this commented block will never be activated.  To be absolutely safe this should be
       uncommented, but it is most likely unnecessary.  
    */
        if (qsqr < evecprec)
        {
            q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
            q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
            q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
            q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < evecprec)
            {
                double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
                double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
                double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

                q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
                q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
                q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
                q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

                if (qsqr < evecprec)
                {
                    q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                    q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                    q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                    q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                    qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                    
                    if (qsqr < evecprec)
                    {
                        /* if qsqr is still too small, return the identity matrix. */
                        rotmat.setIdentity();

                        return 0;
                    }
                }
            }
        }

        double normq = Math.sqrt(qsqr);
        q1 /= normq;
        q2 /= normq;
        q3 /= normq;
        q4 /= normq;
        
        System.out.println("q: " + q1 + " " + q2 + " " + q3 + " " + q4);

        double a2 = q1 * q1;
        double x2 = q2 * q2;
        double y2 = q3 * q3;
        double z2 = q4 * q4;

        double xy = q2 * q3;
        double az = q1 * q4;
        double zx = q4 * q2;
        double ay = q1 * q3;
        double yz = q3 * q4;
        double ax = q1 * q2;

        rotmat.m00 = a2 + x2 - y2 - z2;
        rotmat.m01 = 2 * (xy + az);
        rotmat.m02 = 2 * (zx - ay);

        rotmat.m10 = 2 * (xy - az);
        rotmat.m11 = a2 - x2 + y2 - z2;
        rotmat.m12 = 2 * (yz + ax);

        rotmat.m20 = 2 * (zx + ay);
        rotmat.m21 = 2 * (yz - ax);
        rotmat.m22 = a2 - x2 - y2 + z2;

        return 1;
    }
	
}
