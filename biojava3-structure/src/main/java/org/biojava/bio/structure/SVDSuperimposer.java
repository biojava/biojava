/*
 *                  BioJava development code
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
 * Created on Dec 4, 2005
 *
 */
package org.biojava.bio.structure;

import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.jama.SingularValueDecomposition;


/** A class that calculates the superimposition between two sets of atoms
 * inspired by the biopython SVDSuperimposer class...
 *
 *
 * example usage:
 * <pre>
 * try{

            // get some arbitrary amino acids from somewhere
            String filename   =  "/Users/ap3/WORK/PDB/5pti.pdb" ;

            PDBFileReader pdbreader = new PDBFileReader();
            Structure struc = pdbreader.getStructure(filename);
            Group g1 = (Group)struc.getChain(0).getGroup(21).clone();
            Group g2 = (Group)struc.getChain(0).getGroup(53).clone();

            if ( g1.getPDBName().equals("GLY")){
                if ( g1 instanceof AminoAcid){
                    Atom cb = Calc.createVirtualCBAtom((AminoAcid)g1);
                    g1.addAtom(cb);
                }
            }

            if ( g2.getPDBName().equals("GLY")){
                if ( g2 instanceof AminoAcid){
                    Atom cb = Calc.createVirtualCBAtom((AminoAcid)g2);
                    g2.addAtom(cb);
                }
            }

            Structure struc2 = new StructureImpl((Group)g2.clone());

            System.out.println(g1);
            System.out.println(g2);


            Atom[] atoms1 = new Atom[3];
            Atom[] atoms2 = new Atom[3];

            atoms1[0] = g1.getAtom("N");
            atoms1[1] = g1.getAtom(" CA ");
            atoms1[2] = g1.getAtom("CB");


            atoms2[0] = g2.getAtom("N");
            atoms2[1] = g2.getAtom(" CA ");
            atoms2[2] = g2.getAtom("CB");


            SVDSuperimposer svds = new SVDSuperimposer(atoms1,atoms2);


            Matrix rotMatrix = svds.getRotation();
            Atom tranMatrix = svds.getTranslation();


            // now we have all the info to perform the rotations ...

            Calc.rotate(struc2,rotMatrix);

            //          shift structure 2 onto structure one ...
            Calc.shift(struc2,tranMatrix);

            //
            // write the whole thing to a file to view in a viewer

            String outputfile = "/Users/ap3/WORK/PDB/rotated.pdb";

            FileOutputStream out= new FileOutputStream(outputfile);
            PrintStream p =  new PrintStream( out );

            Structure newstruc = new StructureImpl();

            Chain c1 = new ChainImpl();
            c1.setName("A");
            c1.addGroup(g1);
            newstruc.addChain(c1);

            Chain c2 = struc2.getChain(0);
            c2.setName("B");
            newstruc.addChain(c2);

            // show where the group was originally ...
            Chain c3 = new ChainImpl();
            c3.setName("C");
            //c3.addGroup(g1);
            c3.addGroup(g2);

            newstruc.addChain(c3);
            p.println(newstruc.toPDB());

            p.close();

            System.out.println("wrote to file " + outputfile);

        } catch (Exception e){
            e.printStackTrace();
        }
        </pre>
 *
 *
 * @author Andreas Prlic
 * @since 1.5
 * @version %I% %G%

 */
public class SVDSuperimposer {

    Matrix rot;
    Matrix tran;

    Matrix centroidA;
    Matrix centroidB;

    /** Create a SVDSuperimposer object and calculate a SVD superimposition of two sets of atoms.
     *
     * @param atomSet1 Atom array 1
     * @param atomSet2 Atom array 2
     * @throws StructureException
     */
    public SVDSuperimposer(Atom[] atomSet1,Atom[]atomSet2)
    throws StructureException{

        if ( atomSet1.length != atomSet2.length ){
            throw new StructureException("The two atom sets are not of same length!");
        }

        Atom cena = Calc.getCentroid(atomSet1);
        Atom cenb = Calc.getCentroid(atomSet2);

        double[][] centAcoords = new double[][]{{cena.getX(),cena.getY(),cena.getZ()}};
        centroidA = new Matrix(centAcoords);

        double[][] centBcoords = new double[][]{{cenb.getX(),cenb.getY(),cenb.getZ()}};
        centroidB = new Matrix(centBcoords);

        //      center at centroid

        Atom[] ats1 = Calc.centerAtoms(atomSet1,cena);
        Atom[] ats2 = Calc.centerAtoms(atomSet2,cenb);

        double[][] coordSet1 = new double[ats1.length][3];
        double[][] coordSet2 = new double[ats2.length][3];

        // copy the atoms into the internal coords;
        for (int i =0 ; i< ats1.length;i++) {
            coordSet1[i] = ats1[i].getCoords();
            coordSet2[i] = ats2[i].getCoords();
        }

        calculate(coordSet1,coordSet2);


    }





    /** Do the actual calculation.
     *
     * @param coordSet1 coordinates for atom array 1
     * @param coordSet2 coordiantes for atom array 2
     */
    private void calculate(double[][] coordSet1, double[][]coordSet2){
        // now this is the bridge to the Jama package:
        Matrix a = new Matrix(coordSet1);
        Matrix b = new Matrix(coordSet2);


//      # correlation matrix

        Matrix b_trans = b.transpose();
        Matrix corr = b_trans.times(a);


        SingularValueDecomposition svd = corr.svd();

        // A = U*S*V'.


        Matrix u = svd.getU();
        // v is alreaady transposed ! difference to numermic python ...
        Matrix vt =svd.getV();

        Matrix vt_orig = (Matrix) vt.clone();
        Matrix u_transp = u.transpose();

        Matrix rot_nottrans = vt.times(u_transp);
        rot = rot_nottrans.transpose();

        // check if we have found a reflection

        //printMatrix(rot);

        double det = rot.det();
        //System.out.println(det);

         if (det<0) {
            vt = vt_orig.transpose();
            vt.set(2,0,(0 - vt.get(2,0)));
            vt.set(2,1,(0 - vt.get(2,1)));
            vt.set(2,2,(0 - vt.get(2,2)));

            Matrix nv_transp = vt.transpose();
            rot_nottrans = nv_transp.times(u_transp);
            rot = rot_nottrans.transpose();

        }

        Matrix cb_tmp = centroidB.times(rot);
        tran = centroidA.minus(cb_tmp);


    }

    /** Calculate the RMS (root mean square) deviation of two sets of atoms.
     *
     * Atom sets must be pre-rotated.
     * 
     * @param atomSet1 atom array 1
     * @param atomSet2 atom array 2
     * @return the RMS of two atom sets
     * @throws StructureException
     */
    public static double getRMS(Atom[] atomSet1, Atom[] atomSet2) throws StructureException {
        if ( atomSet1.length != atomSet2.length ){
            throw new StructureException("The two atom sets are not of same length!");
        }

        double sum = 0.0;
        for ( int i =0 ; i < atomSet1.length;i++){
            double d = Calc.getDistance(atomSet1[i],atomSet2[i]);
            sum += (d*d);

        }

        double avd = ( sum/ atomSet1.length);
        //System.out.println("av dist" + avd);
        return Math.sqrt(avd);


    }

    /**
     * Calculate the TM-Score for the superposition.
     * 
     * Atom sets must be pre-rotated.
     * 
     * <p>Citation:<br/>
     * <i>Zhang Y and Skolnick J (2004). "Scoring function for automated assessment 
     * of protein structure template quality". Proteins 57: 702 - 710.</i>
     * 
     * @param atomSet1 atom array 1
     * @param atomSet2 atom array 2
     * @param len1 The full length of the protein supplying atomSet1
     * @param len2 The full length of the protein supplying atomSet2
     * @return The TM-Score 
     * @throws StructureException
     */
    public static double getTMScore(Atom[] atomSet1, Atom[] atomSet2, int len1, int len2) throws StructureException {
        if ( atomSet1.length != atomSet2.length ){
            throw new StructureException("The two atom sets are not of same length!");
        }
        if ( atomSet1.length > len1 ){
            throw new StructureException("len1 must be greater or equal to the alignment length!");
        }
        if ( atomSet2.length > len2 ){
            throw new StructureException("len2 must be greater or equal to the alignment length!");
        }
        
        int Lmin = Math.min(len1,len2);
        int Laln = atomSet1.length;
        
        double d0 = 1.24 * Math.cbrt(Lmin - 15.) - 1.8;
        double d0sq = d0*d0;
        
        double sum = 0;
        for(int i=0;i<Laln;i++) {
            double d = Calc.getDistance(atomSet1[i],atomSet2[i]);
        	sum+= 1./(1+d*d/d0sq);
        }
        
        return sum/Lmin;
    }

    /**  Get the Rotation matrix that is required to superimpose the two atom sets.
     *
     * @return a rotation matrix.
     */
    public Matrix getRotation(){
        return rot;
    }

    /** Get the shift vector.
     *
     * @return the shift vector
     */
    public Atom getTranslation(){

        Atom a = new AtomImpl();
        a.setX(tran.get(0,0));
        a.setY(tran.get(0,1));
        a.setZ(tran.get(0,2));
        return a;
    }

    /** Simple debug method to print a Matrix object on System.out.
     *
     * @param m a Matrix
     */
    public void printMatrix(Matrix m){
        for (int i = 0 ; i < m.getRowDimension(); i++){
            for (int j = 0 ; j< m.getColumnDimension(); j++){
                System.out.print("\t" + m.get(i,j) + " ");
            }
            System.out.println("");
        }
    }

}
