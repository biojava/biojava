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
 * Created on Mar 1, 2006
 *
 */
package org.biojava.bio.structure.align.pairwise;

import org.biojava.bio.structure.align.StrucAligParameters;
import org.biojava.bio.structure.align.helper.AligMatEl;
import org.biojava.bio.structure.jama.Matrix;

public class AligNPE {

    public AligNPE() {
        super();

    }
    
    /** Align w/o penalizing endpags. Return alignment and score
     * 
     * @param sim the similarity matrix   
     * @param params the structure alignment parameters to be used
     * @return an Alignable
     */
    public static Alignable align_NPE(Matrix sim,StrucAligParameters params){
        //System.out.println("align_NPE");
        
        float gapOpen = params.getGapOpen();
        float gapExtension = params.getGapExtension();
        
        int rows = sim.getRowDimension();
        int cols = sim.getColumnDimension();
        
        Alignable al = new StrCompAlignment(rows,cols);
        al.setGapExtCol(gapExtension);
        al.setGapExtRow(gapExtension);
        al.setGapOpenCol(gapOpen);
        al.setGapOpenRow(gapOpen);
        //System.out.println("size of aligmat: " + rows+1 + " " + cols+1);
        //AligMatEl[][] aligmat = new AligMatEl[rows+1][cols+1];
        AligMatEl[][] aligmat = al.getAligMat();
        
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {                
             
                int e=0;
                //if ( ( i < rows) &&
                //        ( j < cols)) {
                    //TODO: the ALIGFACTOR calc should be hidden in Gotoh!!
                
                e = (int)Math.round(Gotoh.ALIGFACTOR * sim.get(i,j));
                //}
                //System.out.println(e);
                AligMatEl am = new AligMatEl();
                am.setValue(e);
                //am.setTrack(new IndexPair((short)-99,(short)-99));
                aligmat[i+1][j+1] = am;
                
            }
        }
        //al.setAligMat(aligmat);
        
         new Gotoh(al);
        
        return al;
    }
    
    

}
