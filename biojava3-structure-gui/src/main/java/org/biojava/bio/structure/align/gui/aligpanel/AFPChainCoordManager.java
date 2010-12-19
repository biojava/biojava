package org.biojava.bio.structure.align.gui.aligpanel;

import java.awt.Point;

import org.biojava.bio.structure.align.model.AFPChain;

public class AFPChainCoordManager {

	
	AFPChain afpChain ;

	/** Space on the right side between sequence and legend.
	 * 
	 */
	public static final int DEFAULT_RIGHT_SPACER = 10;
	
	/** number of chars per line
	 *  
	 */
	public static final int DEFAULT_LINE_LENGTH = 70;
	
	/** size of space between rows
	 * 
	 */
	public static final int DEFAULT_Y_STEP = 60;
	
	/** size per character
	 * 
	 */
	public static final int DEFAULT_CHAR_SIZE = 12;
	
	
	/** separation of line 1 and 2 in alignment
	 * 
	 */
	public static final int DEFAULT_LINE_SEPARATION = 20;
	
	
	/** left boundary
	 * 
	 */
	public static final int DEFAULT_X_SPACE = 20;

	/** top boundary
	 * 
	 */
	public static final int DEFAULT_Y_SPACE = 40;

	/** Position at which the alignment summary is printed
	 * 
	 */
	public static final int SUMMARY_POS = 20;
	
	private static final int DEFAULT_LEGEND_SIZE = 50;
	
	
	public int getSummaryPos(){
		return SUMMARY_POS;
	}
	
	/** X coordinate size
	 * 
	 * @return the preferred width
	 */
	public int getPreferredWidth(){
		return 2* DEFAULT_X_SPACE + DEFAULT_LINE_LENGTH * DEFAULT_CHAR_SIZE + DEFAULT_LEGEND_SIZE +DEFAULT_RIGHT_SPACER + DEFAULT_LEGEND_SIZE;
	}
	
	/** Y coordinate size
	 * 
	 * @return the preferred height
	 */
	public int getPreferredHeight(){
		return 2* DEFAULT_Y_SPACE + (afpChain.getAlnLength() / DEFAULT_LINE_LENGTH) * DEFAULT_Y_STEP + DEFAULT_LINE_SEPARATION;
	}
	
	/** Convert from a X position in the JPanel to alignment position
	 * 
	 * @param aligSeq sequence 0 or 1 
	 * @param p point on panel
	 * @return the sequence position for a point on the Panel
	 */
	public int getSeqPos(int aligSeq, Point p) {
		
		int x = p.x - DEFAULT_X_SPACE - DEFAULT_LEGEND_SIZE;
		int y = p.y - DEFAULT_Y_SPACE ;
					
		y -=  (DEFAULT_LINE_SEPARATION * aligSeq) - DEFAULT_CHAR_SIZE  ;
			
		int lineNr = y / DEFAULT_Y_STEP;
		
		//System.out.println("line : " + lineNr);
		
		int linePos = x / DEFAULT_CHAR_SIZE;
		
		//System.out.println("line : " + lineNr + " pos in line: " + linePos);
		return lineNr * DEFAULT_LINE_LENGTH + linePos ;
		
	}

	/** get the position of the sequence position on the Panel
	 * 
	 * @param aligSeq  0 or 1 for which of the two sequences to ask for.
	 * @param i sequence position
	 * @return the point on a panel for a sequence position
	 */
	public Point getPanelPos(int aligSeq, int i) {
		Point p = new Point();
		
		// get line
		// we do integer division since we ignore remainders
		int lineNr = i / DEFAULT_LINE_LENGTH;
		
		// but we want to have the reminder for the line position.
		int linePos = i % DEFAULT_LINE_LENGTH;
		
		int x = linePos * DEFAULT_CHAR_SIZE + DEFAULT_X_SPACE + DEFAULT_LEGEND_SIZE;
		
		int y = lineNr * DEFAULT_Y_STEP + DEFAULT_Y_SPACE;
		
		
		y += DEFAULT_LINE_SEPARATION * aligSeq;
		
		p.setLocation(x, y);
		return p;
	}

	public void setAFPChain(AFPChain afpChain) {
		this.afpChain = afpChain;
		
	}

	/** returns the AligSeq (0 or 1) for a point
	 * returns -1 if not over an alig seq.
	 * @param point
	 * @return which of the two sequences a point on the panel corresponds to
	 */
	public int getAligSeq(Point point) {
		
		
		int i1 = getSeqPos(0, point);
		
		
		Point t1 = getPanelPos(0,i1);
		
		if ( Math.abs(t1.x - point.x) <= DEFAULT_CHAR_SIZE &&
				Math.abs(t1.y-point.y) < DEFAULT_CHAR_SIZE ) {
			return 0;
		}
		
		int i2   = getSeqPos(1,point);
		Point t2 = getPanelPos(1,i2);
		
		if ( Math.abs(t2.x - point.x) < DEFAULT_CHAR_SIZE &&
				Math.abs(t2.y-point.y) < DEFAULT_CHAR_SIZE ) {
			return 1;
		}
		
		//System.out.println(" i1: " + i1 +" t1 : " + Math.abs(t1.x - point.x) + " " + Math.abs(t1.y-point.y));
		//System.out.println(i2);
		
		
		return -1;
	}

	/** provide the coordinates for where to draw the legend for line X and if it is chain 1 or 2
	 * 
	 * @param lineNr which line is this for
	 * @param chainNr is it chain 0 or 1
	 * @return get the point where to draw the legend
	 */
	public Point getLegendPosition(int lineNr, int chainNr) {
		int x = DEFAULT_X_SPACE ;
		
		int y = lineNr * DEFAULT_Y_STEP + DEFAULT_Y_SPACE;
		 
		y += chainNr * DEFAULT_LINE_SEPARATION;
		
		Point p = new Point(x,y);
		return p;
	}
	
	public Point getEndLegendPosition(int lineNr, int chainNr) {
		
		int x = DEFAULT_LINE_LENGTH * DEFAULT_CHAR_SIZE + DEFAULT_X_SPACE + DEFAULT_LEGEND_SIZE + DEFAULT_RIGHT_SPACER;
		
		int y = lineNr * DEFAULT_Y_STEP + DEFAULT_Y_SPACE;
		 
		y += chainNr * DEFAULT_LINE_SEPARATION;
		
		Point p = new Point(x,y);
		return p;
	}

}
