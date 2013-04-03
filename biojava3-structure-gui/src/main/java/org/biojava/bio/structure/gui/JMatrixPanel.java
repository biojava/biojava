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
 * Created on Aug 3, 2007
 * 
 */

package org.biojava.bio.structure.gui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.lang.reflect.Method;

import javax.swing.JPanel;

import org.biojava.bio.structure.align.StrucAligParameters;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.align.pairwise.FragmentPair;
import org.biojava.bio.structure.gui.util.color.ContinuousColorMapper;
import org.biojava.bio.structure.gui.util.color.DefaultMatrixMapper;
import org.biojava.bio.structure.jama.Matrix;


/** a JPanel that can display a difference of distance matrix and paths that have been
 * taken for the alignment
 * 
 * <p>Note: This panel displays the transpose of its underlying matrix.
 * Thus its width will be the same as {@link Matrix#getRowDimension()} and its
 * height the same as {@link Matrix#getColumnDimension()}. This stems from the
 * unfortunate ordering of {@link Matrix#get(int, int)} parameters as (row, col),
 * which is opposite from the normal (x,y) order used when displaying graphics.
 *  
 * @author Andreas Prlic
 *
 */
public class JMatrixPanel extends JPanel{

	/**
	 * 
	 */
	private static final long serialVersionUID = -1720879395453257846L;
	BufferedImage _bufImage;
	Matrix matrix;
	ContinuousColorMapper cellColor; //Maps matrix elements to a color
	float scale;

	FragmentPair[] fragmentPairs;
	AlternativeAlignment[] aligs;
	int selectedAlignmentPos;
	
	final static BasicStroke stroke = new BasicStroke(2.0f);
	StrucAligParameters params;
	
	public JMatrixPanel(){
		scale = 1;
		cellColor = new DefaultMatrixMapper(10, 0.9f);
		//saturation = 0.9f;
		//scalevalue = 10;
		selectedAlignmentPos = -1;
		matrix = new Matrix(0,0);
		params = new StrucAligParameters();
	}

	public int getSelectedAlignmentPos() {
		return selectedAlignmentPos;
	}

	public void setSelectedAlignmentPos(int selectedAlignmentPos) {
		this.selectedAlignmentPos = selectedAlignmentPos;
	}

	public AlternativeAlignment[] getAlternativeAligs() {
		return aligs;
	}

	public void setAlternativeAligs(AlternativeAlignment[] aligs) {
		this.aligs = aligs;
	}



	public FragmentPair[] getFragmentPairs() {
		return fragmentPairs;
	}

	public void setFragmentPairs(FragmentPair[] fragmentPairs) {
		this.fragmentPairs = fragmentPairs;
	}

	public float getScale() {
		return scale;
	}

	public void setPreferredSize(){

		int prefW = Math.round(matrix.getRowDimension() * scale);
		int prefH = Math.round(matrix.getColumnDimension() * scale);

		this.setPreferredSize(new Dimension(prefW,prefH));

	}

	public void setScale(float scale) {

		if ( scale == this.scale)
			return;
		//System.out.println("setting scale " + scale + "current width " + getWidth() + " " + getHeight());

		this.scale = scale;

		setPreferredSize();

		this.repaint();

	}

	public Matrix getMatrix() {
		return matrix;
	}

	/** sets the distance matrix to be displayed
	 * 
	 * @param matrix
	 */
	public void setMatrix(Matrix matrix) {
		this.matrix = matrix;
		setPreferredSize();
	}

	public void paintComponent(Graphics g){

		//super.paintComponent(g);

		Graphics2D g2 = (Graphics2D)g;
		if ( _bufImage == null){						

			int w = getWidth();
			int h = getHeight();
			_bufImage = (BufferedImage) createImage(w,h);
			//Graphics gc = _bufImage.createGraphics();
			//gc.setColor(Color.blue);
			//gc.fillRect(0,0,w,h);

		}


		g2.drawImage(_bufImage,null,0,0);
		drawDistances(g);

		drawPairs(g);

		if ( scale >= 4) {
			drawBoxes(g);
		}
	}

	/** draw alternative alignments
	 * 
	 * @param g
	 */
	public void drawPairs(Graphics g){
		
		if ( aligs == null)
			return;
		
		int nr = aligs.length;
		
		Graphics2D g2D = (Graphics2D)g;
		Stroke oldStroke = g2D.getStroke();
		g2D.setStroke(stroke);
		
		Color color;
		float hue;
		
		int width = Math.round(scale);
		int w2 = width / 2 ;
		
		for (int i = 0; i < aligs.length; i++) {
			AlternativeAlignment a = aligs[i];
			int[] idx1 = a.getIdx1();
			int[] idx2 = a.getIdx2();
			int xold = -1;
			int yold = -1;
			boolean start = true;
			
			if ( (selectedAlignmentPos != -1 ) && 
					( selectedAlignmentPos == i)){
				color = Color.white;
			} else {
			
				hue = i * (1/ (float)nr);
				color = Color.getHSBColor(hue,1.0f,1.0f);
			}
			g.setColor(color);
			
			for (int j = 0; j < idx1.length; j++) {
				int x1 = Math.round(idx1[j]*scale) ;
				int y1 = Math.round(idx2[j]*scale) ;
				if ( ! start){
					//g.drawLine(xold+1,yold,x1+1,y1);
					
					g2D.draw(new Line2D.Double(xold,yold,x1,y1));
					g.fillRect(xold,yold,2,2);
				} else {
					g.fillRect(x1,y1, w2, w2);
					start =false;
				}
				xold = x1;
				yold = y1;
			}
						
			if ( ! start)
				g.fillRect(xold,yold,w2,w2);
			
			
		}
		
		g2D.setStroke(oldStroke);
	}
	

	/** draw high scoring fragments that are used for the initial alignment seed 
	 * selection
	 * 
	 * @param g
	 */
	public void drawBoxes(Graphics g){
		if ( fragmentPairs == null )
			return;
		
		g.setColor(Color.yellow);
		
		
		for (int i = 0; i < fragmentPairs.length; i++) {
			FragmentPair fp =fragmentPairs[i];
			int xp = fp.getPos1();
			int yp = fp.getPos2();

			int width = Math.round(scale);

			g.drawRect(Math.round(xp*scale),Math.round(yp*scale),width, width);

		}
	}


	/**
	 * For each element in matrix, draw it as a colored square or pixel.
	 * 
	 * The color of a matrix element with value x is specified as
	 *   - H: 1-x/scalevalue
	 *   - S: saturation
	 *   - B: 1-x/scalevalue
	 * @param g1
	 */
	public void drawDistances(Graphics g1){
		Graphics2D g = (Graphics2D)g1;

		int c = matrix.getRowDimension();
		int d = matrix.getColumnDimension();

		float scale = getScale();
		int width = Math.round(scale);

		for (int i = 0; i < c; i++) {
			int ipaint = Math.round(i*scale);

			for (int j = 0; j < d; j++) {
				double val = matrix.get(i,j);

				int jpaint = Math.round(j*scale);

				Color color = cellColor.getColor(val);
				g.setColor(color);

				g.fillRect(ipaint,jpaint,width,width);
			}

		}

	}

	
	@Deprecated
	public float getSaturation() {
		try {
			Method getSaturation = cellColor.getClass().getMethod("getSaturation");
			Float saturation = (Float)getSaturation.invoke(cellColor);
			return saturation;
		} catch (Exception e) {
			throw new IllegalStateException("Error calling getSaturation() for "+cellColor.getClass());
		}
	}
	
	@Deprecated
	public void setSaturation(float saturation) {
		try {
			Method setSaturation = cellColor.getClass().getMethod("setSaturation",Float.TYPE);
			setSaturation.invoke(cellColor);
		} catch (Exception e) {
			throw new IllegalStateException("Error calling setSaturation(float) for "+cellColor.getClass());
		}
	}

	@Deprecated
	public float getScalevalue() {
		try {
			Method getScalevalue = cellColor.getClass().getMethod("getScalevalue");
			Float scalevalue = (Float)getScalevalue.invoke(cellColor);
			return scalevalue;
		} catch (Exception e) {
			throw new IllegalStateException("Error calling getScalevalue() for "+cellColor.getClass());
		}
	}

	@Deprecated
	public void setScalevalue(float scalevalue) {
		try{
			Method setScalevalue = cellColor.getClass().getMethod("setScalevalue",Float.TYPE);
			setScalevalue.invoke(cellColor);
		} catch (Exception e) {
			throw new IllegalStateException("Error calling setScalevalue(float) for "+cellColor.getClass());
		}
	}

	/**
	 * @return the color mapping of the JMatrixPanel
	 */
	public ContinuousColorMapper getCellColor() {
		return cellColor;
	}

	/**
	 * @param cellColor the color mapping of the JMatrixPanel to set
	 */
	public void setCellColor(ContinuousColorMapper cellColor) {
		this.cellColor = cellColor;
	}
	


}
