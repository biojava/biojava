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
import org.biojava.nbio.structure.align.helper.GapArray;
import org.biojava.nbio.structure.align.helper.IndexPair;

import java.util.ArrayList;
import java.util.List;

/** a class to perform Gotoh algorithm
 *
 * @author Andreas Prlic (Java), Peter Lackner (original C code)
 * @since 10:56:53 AM
 * @version %I% %G%
 */
public class Gotoh {
    public static int ALIGFACTOR = 1000; // constant to shift floats to ints

    Alignable a;

    int k,openPen,elgPen,rowDim,colDim,openVal,elgVal;

    AligMatEl currentCell;
    GapArray currentGap;



    public Gotoh(Alignable alignable) {
        super();
        a = alignable;

        align();

    }



    private void align() {

        rowDim = a.getRows()+1;
        colDim = a.getCols()+1;

        openPen = Math.round(ALIGFACTOR * a.getGapOpenCol());
        elgPen  = Math.round(ALIGFACTOR * a.getGapExtCol());

        GapArray[] gapCol = new GapArray[colDim];
        GapArray[] gapRow = new GapArray[rowDim];

        // init the arrays
        for ( int i = 0 ; i< colDim;i++){
            gapCol[i] = new GapArray();
        }
        for ( int j = 0 ; j< rowDim;j++){
            gapRow[j] = new GapArray();
        }
        currentGap = new GapArray();

        AligMatEl[][]aligmat = a.getAligMat();
        int lastValue = aligmat[rowDim-1][colDim-1].getValue();

//      first cell
        aligmat[0][0].setValue(0);
        gapCol[0].setValue(0);
        gapCol[0].setIndex(0);
        gapRow[0].setValue(0);
        gapRow[0].setIndex(0);

        // set row 0 margin
        for(int j=1;j<colDim;j++) {
            aligmat[0][j].setValue(0);
            gapCol[j].setValue(-(openPen+elgPen));
            gapCol[j].setIndex(0);
        }


        for(int rowCounter=1;rowCounter<rowDim;rowCounter++){

            // set column 0 margin
            aligmat[rowCounter][0].setValue(0);
            gapRow[rowCounter].setValue(-(openPen+elgPen));
            gapRow[rowCounter].setIndex(0);

            for(int colCounter=1;colCounter<colDim;colCounter++) {

                currentCell = aligmat[rowCounter][colCounter];

                // no gap 
                currentCell.setValue( aligmat[rowCounter-1][colCounter-1].getValue() + currentCell.getValue());
                currentCell.setRow((short)(rowCounter-1));
                currentCell.setCol((short)(colCounter-1));

                // scan column j for gap, gap in seqB
                openVal = aligmat[rowCounter-1][colCounter].getValue() - (openPen+elgPen);
                elgVal  = gapCol[colCounter].getValue()-elgPen;

                currentGap = new GapArray();

                if ( openVal >= elgVal){
                    currentGap.setValue(openVal);
                    currentGap.setIndex(rowCounter-1);
                } else {
                    currentGap.setValue(elgVal);
                    currentGap.setIndex(gapCol[colCounter].index);
                }
                gapCol[colCounter] = currentGap;

                if (currentGap.getValue() > currentCell.getValue()){
                    if ( currentGap.getIndex() >= rowDim)
                        System.err.println("col gap at" + rowCounter + " " + colCounter + " to " + currentGap.getIndex());
                    currentCell.setValue( currentGap.getValue());
                    currentCell.setRow((short)currentGap.getIndex());
                    currentCell.setCol((short)colCounter);
                }

//              scan row i for gap, gap in row
                openVal = aligmat[rowCounter][colCounter-1].getValue()-(openPen+elgPen);
                elgVal  = gapRow[rowCounter].getValue() - elgPen;

                currentGap = new GapArray();

                if (openVal >= elgVal){
                    currentGap.setValue(openVal);
                    currentGap.setIndex(colCounter-1);
                } else {
                    currentGap.setValue(elgVal);
                    currentGap.setIndex(gapRow[rowCounter].getIndex());
                }
                gapRow[rowCounter] = currentGap;


                if ( currentGap.getValue() > currentCell.getValue() ) {
                    if ( currentGap.getIndex() >= colDim)
                        System.err.println("row gap at" + rowCounter + " " + colCounter + " to " + currentGap.getIndex());
                    currentCell.setValue(currentGap.getValue());
                    currentCell.setRow((short)rowCounter);
                    currentCell.setCol((short)currentGap.getIndex());
                }

                aligmat[rowCounter][colCounter]=currentCell;
            }

        }


        // last cell in alignment matrix
        // do not penalize end gaps
        int rowCount = rowDim -1;
        int colCount = colDim -1;
        currentCell = aligmat[rowCount][colCount];

        // no gap
        currentCell.setValue(aligmat[rowCount-1][colCount-1].getValue() + lastValue);
        currentCell.setRow((short)(rowCount-1));
        currentCell.setCol((short)(colCount-1));

        // scan last column j for gap, gap in seqB
        // do not penalyze gaps
        for (int k=1;k<=rowCount;k++) {
            int val = aligmat[rowCount-k][colCount].getValue();
            if ( val>currentCell.getValue()){
                currentCell.setValue(val);
                //System.out.println("setting row to " + (rowCount ) );
                currentCell.setRow((short)(rowCount-k ));
                currentCell.setCol((short)(colCount  ));
            }
        }

        // scan row i for gap, gap in SeqA 
        // do not penalyze gaps 
        for (int k=1;k<=colCount;k++) {
            int val = aligmat[rowCount][colCount-k].getValue();
            if ( val > currentCell.getValue()){
                currentCell.setValue(val);
                currentCell.setRow((short) (rowCount ) );
                currentCell.setCol((short)(colCount-k  ));
            }
        }
        a.setScore(aligmat[rowDim-1][colDim-1].getValue() / (float)ALIGFACTOR);
        setPath();

    }

    private void setPath(){

        // dmyersturnbull: TODO fix exception handling and logging

        int  n;
        IndexPair[] backId = new IndexPair[a.getRows()+1+a.getCols()+1];
        List<IndexPair> path = new ArrayList<IndexPair>();

        backId[0] = new IndexPair((short)(a.getRows()),(short)(a.getCols()));

        // backtrace, get backId indices, the indices in diagonal store in path


        int pathsize = 0;

        AligMatEl[][] aligmat = a.getAligMat();


        n=1;
        while ( (backId[n-1].getRow()>=1) &&
                (backId[n-1].getCol()>=1)
                )
        {
            // get backtrace index
            int x = backId[n-1].getRow();
            int y = backId[n-1].getCol();

            try {

                AligMatEl el = null ;
                try {
                    el =aligmat[x][y];
                } catch(Exception e){


                    e.printStackTrace();
                    for (int f=0; f< n;f++){
                        System.out.println(backId[f]);
                    }

                }

                if ( el == null)
                    System.out.println("el = null! x:"+ x + " y " + y);
                backId[n] = el;
            } catch (Exception e){
                e.printStackTrace();
                System.out.println("x " + x);
                System.out.println("y " + y);
                System.out.println(backId[n-2]);
                System.exit(0);
            }
            // get diagonal indeces into path
            if (((backId[n-1].getRow() - backId[n].getRow()) == 1)
                    && (( backId[n-1].getCol() - backId[n].getCol()) == 1)) {
                path.add(backId[n-1]);
                pathsize++;

            }
            n++;
        }

        // path is reverse order..
        // switch order
        IndexPair[] newpath = new IndexPair[pathsize];
        for (int i = 0 ; i < pathsize; i++){
            IndexPair o = path.get(pathsize-1-i);
            IndexPair np = new IndexPair((short)(o.getRow()-1),(short)(o.getCol()-1));
            newpath[i] = np;
        }

        a.setPath(newpath);
        a.setPathSize(pathsize);


    }
}
