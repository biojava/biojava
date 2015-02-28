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
 * Created on Jan 29, 2006
 *
 */
package org.biojava.nbio.structure.align.pairwise;

import org.biojava.nbio.structure.align.helper.AligMatEl;
import org.biojava.nbio.structure.align.helper.IndexPair;


public class StrCompAlignment 
implements Alignable
{
    
   
    AligMatEl[][] aligmat;
    int rows;
    int cols;
    IndexPair[] path;
    int pathSize;
    float gapOpenCol;
    float gapOpenRow;
    float gapExtCol;
    float gapExtRow;
    float score;
    
    public StrCompAlignment(int rows, int cols){
        //System.out.println("new alignment " + rows + " " + cols);
        this.rows = rows;
        this.cols = cols;
        aligmat = new AligMatEl[rows+1][cols+1];
    
        /*for (int i=0;i<rows+1;i++){
            for(int j=0;j<cols+1;j++){
                aligmat[i][j] = new AligMatEl();
            }
        }*/
        path = new IndexPair[0];
        score = 0;
    }
    
    @Override
	public int getRows(){
        return rows;
    }
    
    @Override
	public int getCols() {
        return cols;
    }
    
    public void setAligMat(int i, int j,AligMatEl el){
        aligmat[i][j] = el;
    }
    public AligMatEl getAligMat(int i,int j){
        return aligmat[i][j];
    }
    
    @Override
	public AligMatEl[][] getAligMat(){
        return aligmat;
    }
    
    @Override
	public void setAligMat(AligMatEl[][] al){
        //System.out.println("setting alig mat: " + al.length + " " + al[0].length);
        rows = al.length -1;
        cols = al[0].length -1;
        aligmat = al;
    }
    
    @Override
	public float getGapExtCol() {
        return gapExtCol;
    }

    @Override
	public void setGapExtCol(float gapExtCol) {
        this.gapExtCol = gapExtCol;
    }

    @Override
	public float getGapExtRow() {
        return gapExtRow;
    }

    @Override
	public void setGapExtRow(float gapExtRow) {
        this.gapExtRow = gapExtRow;
    }

    @Override
	public float getGapOpenCol() {
        return gapOpenCol;
    }

    @Override
	public void setGapOpenCol(float gapOpenCol) {
        this.gapOpenCol = gapOpenCol;
    }

    @Override
	public float getGapOpenRow() {
        return gapOpenRow;
    }

    @Override
	public void setGapOpenRow(float gapOpenRow) {
        this.gapOpenRow = gapOpenRow;
    }

    @Override
	public float getScore() {
      return score;
    }

    @Override
	public void setScore(float score) {
        //System.out.println("StrCompAlig got score " +score);
       this.score = score;
    }

    @Override
	public IndexPair[] getPath() {
        return path;
    }

    @Override
	public void setPath(IndexPair[] path) {
        this.path = path;
    }

    @Override
	public int getPathSize() {
        return pathSize;
    }

    @Override
	public void setPathSize(int pathSize) {
        this.pathSize = pathSize;
    }
    
    
    
}





