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
 * Created on Nov 8, 2005
 *
 */
package org.biojava.bio.structure.gui.util;

import java.awt.AlphaComposite;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.util.logging.*;
import javax.swing.JPanel;
import org.biojava.bio.structure.*;


import java.awt.Color;
import java.util.*;


/** A class that draws a Sequence as a rectangle, a scale display over it.
 * 
 * @author Andreas Prlic
 * @since 1.7
 */
public class SequenceScalePanel
extends JPanel{

	static final long serialVersionUID = 7893248902423l;

	Logger logger = Logger.getLogger("org.biojava");

	public static final int    DEFAULT_X_START          = 10  ;
	public static final int    DEFAULT_X_RIGHT_BORDER   = 40 ;
	public static final int    DEFAULT_Y_START          = 0 ;
	public static final int    DEFAULT_Y_STEP           = 10 ;
	public static final int    DEFAULT_Y_HEIGHT         = 8 ;// the size of the box
	public static final int    DEFAULT_Y_BOTTOM         = 16 ;
	public static final int    LINE_HEIGHT              = 10 ;    
	public static final int    MINIMUM_HEIGHT           = 20;
	public static final Color  SEQUENCE_COLOR           = Color.LIGHT_GRAY;
	public static final Color  SCALE_COLOR              = Color.black;
	public static final Color  TEXT_SCALE_COLOR         = Color.GRAY;
	public static final Color  IDX_COLOR         		= Color.yellow;
	public static final Color  GAP_COLOR         		= Color.white;
	public static final Color  BACKGROUND_COLOR;
	public static final Font   seqFont ;

	// the scale value after which to show the sequence as text    
	private static final int   SEQUENCE_SHOW = 9;

	// the height of the panel
	public static final int SIZE = 20;

	Chain chain;
	int chainLength;
	float scale;
	Character[] seqArr;    

	CoordManager coordManager;

	

	int position;
	List<AlignedPosition> apos;

	static {

		String fontName = "Helvetica";

		int fsize = 12;
		seqFont = new Font(fontName,Font.PLAIN,fsize);

		String col1 = "#FFFFFF";
		BACKGROUND_COLOR = Color.decode(col1);

	}


	public SequenceScalePanel(int position) {
		super();
		this.position = position;
		this.setBackground(BACKGROUND_COLOR);

		chain = new ChainImpl();
		setDoubleBuffered(true);

		seqArr = new Character[0];       
		chainLength = 0;
		scale = 1.0f;

		setPrefSize();
		coordManager = new CoordManager();

		apos = new ArrayList<AlignedPosition>();
		
	}

	


	private void setPrefSize() {

		int length = chainLength  ; 
		int l = Math.round(length*scale) + DEFAULT_X_START + DEFAULT_X_RIGHT_BORDER ;
		if ( l  < 60){
			l = 60;
		}
		this.setPreferredSize(new Dimension(l,SIZE));

	}

	public void setAligMap(List<AlignedPosition> apos){
		this.apos = apos;

		if ( apos.size() == 0)
			return;

		AlignedPosition last = apos.get(apos.size()-1);
		//System.out.println("got last aligned:" +last);
		if ( last.getPos(position) != -1){
			// add the end of the chain...
			int idxlast = last.getPos(position);


			for (;idxlast < chainLength;idxlast++){
				AlignedPosition m = new AlignedPosition();
				m.setPos(position,idxlast);

				apos.add(m);
			}
		}

	}

	public synchronized void setChain(Chain c){

		List<Group> a = c.getAtomGroups("amino");

		seqArr = new Character[a.size()];

		chain = new ChainImpl();

		Iterator<Group> iter = a.iterator();
		int i = 0;
		while (iter.hasNext()){
			AminoAcid aa = (AminoAcid) iter.next();

			// preserver original hierarchy ... for highlighting in Jmol
			Chain old = aa.getChain();
			chain.addGroup(aa);
			aa.setChain(old);
			seqArr[i] = aa.getAminoType();
			i++;
		}

		chainLength = i;
		coordManager.setLength(chainLength);
		setPrefSize();

		this.repaint();  
	}

	public Chain getChain(){
		return chain;
	}

	public synchronized float getScale(){
		return scale;
	}


	public void setScale(float scale) {

		this.scale=scale;
		coordManager.setScale(scale);
		setPrefSize();

		this.repaint();
		this.revalidate();
	}

	/** set some default rendering hints, like text antialiasing on
	 * 
	 * @param g2D the graphics object to set the defaults on
	 */
	protected void setPaintDefaults(Graphics2D g2D){
		g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
				RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);

		g2D.setFont(seqFont);
	}

	public void paintComponent(Graphics g){

		g.setColor(BACKGROUND_COLOR);

		Rectangle drawHere = g.getClipBounds();        
		g.fillRect(drawHere.x,drawHere.y, drawHere.width, drawHere.height);


		Graphics2D g2D =(Graphics2D) g;

		setPaintDefaults(g2D);

		int y = 1;


		// draw the scale        

		y = drawScale(g2D,1);

		// draw the background for identical residues
		drawIdx(g2D,y);

		//  sequence
		y = drawSequence(g2D,y);


	}



	/** draw the Scale
	 * 
	 * @param g2D
	 * @param y the height on which to draw the scale
	 * @return the new y position
	 */
	protected int drawScale(Graphics2D g2D, int y){

		// only draw within the ranges of the Clip
		Rectangle drawHere = g2D.getClipBounds();        

		g2D.setColor(SCALE_COLOR);

		int aminosize = Math.round(1*scale);
		if ( aminosize < 1)
			aminosize = 1;



		int startpos = coordManager.getSeqPos(drawHere.x);       
		int endpos   = coordManager.getSeqPos(drawHere.x+drawHere.width);

		if ( endpos > apos.size())
			endpos = apos.size();

		int l = endpos - startpos + 1 ;     

		int drawStart = coordManager.getPanelPos(startpos);
		int drawEnd   = coordManager.getPanelPos(l) - DEFAULT_X_START + aminosize;

		/*System.out.println("SeqScalePanel drawing scale s:" + startpos + " e: " + endpos + 
             " ps: " + drawStart + " pe:" + drawEnd  + " draw.x " + drawHere.x + " draw.w " + drawHere.width +
             " scale " + scale);
		 */

//		the frame around the sequence box
		if ( scale < SEQUENCE_SHOW){
			g2D.setColor(SEQUENCE_COLOR);
			//g2D.setColor(Color.blue);
			Rectangle seqline = new Rectangle(drawStart, y, drawEnd, LINE_HEIGHT);

			//g2D=  (Graphics2D)g;
			g2D.fill(seqline);   
			//g2D.setColor(Color.blue);
			//g2D.draw(seqline);
		}

		// the top line for the scale
		g2D.setColor(SCALE_COLOR);
		Rectangle baseline = new Rectangle(drawStart, y, drawEnd, 2);        
		g2D.fill(baseline);


		// draw the vertical ticks


		int lineH = 11;
		if ( scale <= 3)
			lineH = 8;

		for (int gap =startpos ; ((gap<= endpos) && ( gap < apos.size())); gap++){
			int xpos = coordManager.getPanelPos(gap) ;

			AlignedPosition m = apos.get(gap);
			if ( m.getPos(position) == -1 ){
				// a gap position
				g2D.setColor(GAP_COLOR);
				g2D.fillRect(xpos, y+2, aminosize+1, y+lineH);
				g2D.setColor(GAP_COLOR);
				continue;
			}

			int i = m.getPos(position);

			if ( ((i+1)%100) == 0 ) {

				if ( scale> 0.1) {
					g2D.setColor(TEXT_SCALE_COLOR);
					g2D.fillRect(xpos, y+2, aminosize, y+lineH);
					g2D.setColor(SCALE_COLOR);
					if ( scale < SEQUENCE_SHOW)
						g2D.drawString(""+(i+1),xpos,y+DEFAULT_Y_STEP);
				}

			}else if  ( ((i+1)%50) == 0 ) {
				if ( scale>1.4) {                    
					g2D.setColor(TEXT_SCALE_COLOR);
					g2D.fillRect(xpos,y+2, aminosize, y+lineH);  
					g2D.setColor(SCALE_COLOR);
					if ( scale < SEQUENCE_SHOW)
						g2D.drawString(""+(i+1),xpos,y+DEFAULT_Y_STEP);

				}

			} else if  ( ((i+1)%10) == 0 ) {                
				if ( scale> 3) {
					g2D.setColor(TEXT_SCALE_COLOR);
					g2D.fillRect(xpos, y+2, aminosize, y+lineH);
					g2D.setColor(SCALE_COLOR);
					if ( scale < SEQUENCE_SHOW)
						g2D.drawString(""+(i+1),xpos,y+DEFAULT_Y_STEP);

				}
			} 
		}


		int length = chainLength;       
		if ( endpos >= length-1) {

			int endPanel = coordManager.getPanelPos(endpos);
			g2D.drawString(""+length,endPanel+10,y+DEFAULT_Y_STEP);
		}

		return y ;

	}

	protected void drawIdx(Graphics2D g2D, int y){

		int aminosize = Math.round(1*scale);
		if ( aminosize < 1)
			aminosize = 1;

		// only draw within the ranges of the Clip
		Rectangle drawHere = g2D.getClipBounds();        
		int startpos = coordManager.getSeqPos(drawHere.x);       
		//int endpos   = coordManager.getSeqPos(drawHere.x+drawHere.width-2);


		Composite oldComp = g2D.getComposite();
		g2D.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,0.8f));  
		//logger.info("paint l " + l + " length " + length );

		if ( startpos < 0)
			startpos = 999;

		g2D.setColor(IDX_COLOR);
		int lineH = 11;
		if ( scale <= 3)
			lineH = 8;

		int i = startpos;


		// display the actual sequence!;
		for ( int gap = startpos ;  gap < apos.size() ;gap++){
			int xpos = coordManager.getPanelPos(gap) ;

			AlignedPosition m = apos.get(gap);
			if ( m.getEquivalent() == AlignedPosition.NOT_ALIGNED){
				// a gap position				
				continue;
			}

			i = m.getPos(position);




			for (AlignedPosition xi : apos ) {
				if (xi.getPos(position)!= -1) 
					if ( i == xi.getPos(position)){
						g2D.fillRect(xpos, y+2, aminosize, y+lineH);
						break;
					}
			}
			// TODO:
			// color amino acids by hydrophobicity


		}     

		g2D.setComposite(oldComp);

	}

	/** draw the Amino acid sequence
	 * 
	 * @param g2D
	 * @param y .. height of line to draw the sequence onto
	 * @return the new y value
	 */
	protected int drawSequence(Graphics2D g2D,  int y){
		//g2D.drawString(panelName,10,10);

		g2D.setColor(SEQUENCE_COLOR);
		int aminosize = Math.round(1*scale);
		if ( aminosize < 1)
			aminosize = 1;

		// only draw within the ranges of the Clip
		Rectangle drawHere = g2D.getClipBounds();        
		int startpos = coordManager.getSeqPos(drawHere.x);       
		//int endpos   = coordManager.getSeqPos(drawHere.x+drawHere.width-2);


		Composite oldComp = g2D.getComposite();
		g2D.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER,0.8f));  
		//logger.info("paint l " + l + " length " + length );

		if ( startpos < 0)
			startpos = 999;

		if ( scale > SEQUENCE_SHOW){
			g2D.setColor(Color.black);


			//g2D.setColor(SCALE_COLOR);

			int i = startpos;

			// display the actual sequence!;
			for ( int gap = startpos ;   gap < apos.size() ;gap++){
				int xpos = coordManager.getPanelPos(gap) ;

				AlignedPosition m = apos.get(gap);
				if (m.getPos(position) == -1){
					// a gap position	
					g2D.drawString("-",xpos+1,y+2+DEFAULT_Y_STEP);
					continue;
				}

				i = m.getPos(position);

				// TODO:
				// color amino acids by hydrophobicity

				g2D.drawString(seqArr[i].toString(),xpos+1,y+2+DEFAULT_Y_STEP);
			}     

//			in full sequence mode we need abit more space to look nice

			y+=2;  
		}
		g2D.setComposite(oldComp);
		y+= DEFAULT_Y_STEP + 2;
		return y;
	}


}
