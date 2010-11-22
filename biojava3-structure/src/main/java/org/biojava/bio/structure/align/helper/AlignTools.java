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
 * Created on May 27, 2006
 *
 */
package org.biojava.bio.structure.align.helper;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.StructureException;
import org.biojava.bio.structure.jama.Matrix;

public class AlignTools {

    /** get a subset of Atoms based by their positions
     *  
     * @param caall
     * @param idx an array where each element is a position of all the Atoms to return
     * @return at Atom[] array
     */
    public static Atom[] getFragmentFromIdxList(Atom[] caall, int[] idx){
        Atom[] subset = new Atom[idx.length];
        
        for (int p = 0 ; p < idx.length;p++){
        
            int pos1 =  idx[p];
            subset[p] =  (Atom) caall[pos1].clone();            
        }
        return subset;
    }
    
    /** get a continue subset of Atoms based by the starting position and the length
     * 
     * @param caall
     * @param pos ... the start position
     * @param fragmentLength .. the length of the subset to extract.
     * @return an Atom[] array
     */
    public static Atom[] getFragment(Atom[] caall, int pos, int fragmentLength){
     
        if ( pos+fragmentLength > caall.length)
            return null;
        
        Atom[] tmp = new Atom[fragmentLength];
        
        for (int i=0;i< fragmentLength;i++){
            tmp[i] = (Atom)caall[i+pos].clone();
        }
        return tmp;
        
    }
    
    
    /** get a continue subset of Atoms based by the starting position and the length
     * does not clone the original atoms.
     * 
     * @param caall
     * @param pos ... the start position
     * @param fragmentLength .. the length of the subset to extract.
     * @return an Atom[] array
     */
    public static Atom[] getFragmentNoClone(Atom[] caall, int pos, int fragmentLength){
     
        if ( pos+fragmentLength > caall.length)
            return null;
        
        Atom[] tmp = new Atom[fragmentLength];
        
        for (int i=0;i< fragmentLength;i++){
            tmp[i] = (Atom)caall[i+pos];
        }
        return tmp;
        
    }
    
    /** get the centroid for the set of atoms starting fromposition pos, length fragmentLenght
     * 
     * @param ca
     * @param pos
     * @param fragmentLength
     * @return an Atom
     */
    public static Atom getCenter(Atom[] ca, int pos, int fragmentLength){
        Atom center = new AtomImpl();
        
        if ( pos+fragmentLength > ca.length) {
            System.err.println("pos ("+pos+"), fragL ("+fragmentLength +") > ca.length"+ca.length);
            return center;
        }
        
        Atom[] tmp = getFragmentNoClone(ca,pos,fragmentLength);
        
        return Calc.getCentroid(tmp);
    }
    
    
    /* Get distances along diagonal k from coordinate array coords.
     * 
     * @param atoms set of atoms to be used
     * @param k number of diagonal to be used
     */
    public static double[] getDiagonalAtK (Atom[] atoms,int k) {
        
        int l = atoms.length;
                
        double[] dk = new double[(l-k)];
        
        for ( int i = 0 ; i< (l-k); i++){
            try {
                double dist = Calc.getDistance(atoms[i],atoms[i+k]);                
                dk[i] = dist; 
            } catch (StructureException e){
                // this should not happen ...
                // atoms should never be null here ...
            }
        }
        
        return dk;
    }

    /**
     * Given distance matrix diagonals dk1, dk2, get the rmsd of a fpair.
     * i,j is the fpair start in mol1 and mol2, l is the length of the fragment
     * 
     * @param dk1 distances of structure 1
     * @param dk2 distance of structure 2
     * @param i position in structure 1
     * @param j position in structure 2
     * @param l length of the fragments
     * @param k diagonal used
     * @return a double
     */
    public static double rms_dk_diag(double[] dk1, double[] dk2, int i, int j, int l, int k) {
//        dk = 0
//        for x in range(l-k):
//            dk += (dk1[x+i]-dk2[x+j])*(dk1[x+i]-dk2[x+j])
//        dk /= (l-k)
//        return math.sqrt(dk)
        
        double dk = 0.0;
        for (int x = 0 ; x < (l-k) ; x++)
            dk += (dk1[x+i]-dk2[x+j])*(dk1[x+i]-dk2[x+j]);
        
        dk /= ( l-k);
        return Math.sqrt(dk);
        
    }
    
    /** matrix of all distances between two sets of 3d coords"
     * 
     * @param ca1
     * @param ca2
     * @return a Matrixd
     */
    public static Matrix getDistanceMatrix(Atom[] ca1, Atom[]ca2){
        
      int r = ca1.length;
      int c = ca2.length;
      
      Matrix out = new Matrix(r,c);
      
      for (int i=0; i<r; i++) {
          Atom a1 = ca1[i];
          for (int j=0;j<c;j++){
              Atom b1 = ca2[j];
            
//              double d = Math.sqrt(
//                      (a1.getX()-b1.getX()) * (a1.getX()-b1.getX()) + 
//                      (a1.getY()-b1.getY()) * (a1.getY()-b1.getY()) +
//                      (a1.getZ()-b1.getZ()) * (a1.getZ()-b1.getZ()));    
              try {
                  double d = Calc.getDistance(a1,b1);
                  out.set(i,j,d);
              } catch (StructureException e) {
                  e.printStackTrace();
                  out.set(i,j,999);
              }
          }
      }
      return out;  
    }   
    

}
