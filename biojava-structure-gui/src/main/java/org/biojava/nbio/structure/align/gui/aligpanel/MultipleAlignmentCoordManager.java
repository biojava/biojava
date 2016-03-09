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
package org.biojava.nbio.structure.align.gui.aligpanel;

import java.awt.*;

/**
 * Generalization of the Coodinate Manager to include an arbitrary number of
 * sequences (lines) for MultipleAlignment visualization.
 *
 * @author Aleix Lafita
 *
 */
public class MultipleAlignmentCoordManager {

	private int alignmentLength;     	//number of aligned residues
	private int alignmentSize;			//number of strucures aligned

	/**
	 * Constructor.
	 * @param size number of structures/sequences aligned (rows).
	 * @param length number of aligned residues (columns)
	 */
	public MultipleAlignmentCoordManager(int size, int length){
		alignmentLength = length;
		alignmentSize = size;
		DEFAULT_Y_STEP = 30*size;
	}

	/**
	 * Space on the right side between sequence and legend.
	 */
	public static final int DEFAULT_RIGHT_SPACER = 10;

	/**
	 * Number of chars per line
	 */
	public static final int DEFAULT_LINE_LENGTH = 70;

	/**
	 * Size of space between rows.
	 * Depends on the number of structures aligned.
	 */
	public final int DEFAULT_Y_STEP;

	/**
	 * Size per character
	 */
	public static final int DEFAULT_CHAR_SIZE = 12;

	/**
	 * Separation between sequences in the alignment
	 */
	public static final int DEFAULT_LINE_SEPARATION = 20;

	/**
	 * Left boundary
	 */
	public static final int DEFAULT_X_SPACE = 20;

	/**
	 * Top boundary
	 */
	public static final int DEFAULT_Y_SPACE = 40;

	/**
	 * Position at which the alignment summary is printed
	 */
	public static final int SUMMARY_POS = 20;

	private static final int DEFAULT_LEGEND_SIZE = 50;

	public int getSummaryPos(){
		return SUMMARY_POS;
	}

	/**
	 * X coordinate size
	 *
	 * @return the preferred width
	 */
	public int getPreferredWidth(){
		return alignmentSize* DEFAULT_X_SPACE + DEFAULT_LINE_LENGTH *
				DEFAULT_CHAR_SIZE + DEFAULT_LEGEND_SIZE +
				DEFAULT_RIGHT_SPACER + DEFAULT_LEGEND_SIZE;
	}

	/**
	 * Y coordinate size
	 *
	 * @return the preferred height
	 */
	public int getPreferredHeight(){
		return alignmentSize* DEFAULT_Y_SPACE +
				(alignmentLength / DEFAULT_LINE_LENGTH) *
				DEFAULT_Y_STEP + DEFAULT_LINE_SEPARATION;
	}

	/**
	 * Convert from an X position in the JPanel to the position
	 * in the sequence alignment.
	 *
	 * @param aligSeq sequence number
	 * @param p point on panel
	 * @return the sequence position for a point on the Panel
	 */
	public int getSeqPos(int aligSeq, Point p) {

		int x = p.x - DEFAULT_X_SPACE - DEFAULT_LEGEND_SIZE;
		int y = p.y - DEFAULT_Y_SPACE ;

		y -=  (DEFAULT_LINE_SEPARATION * aligSeq) - DEFAULT_CHAR_SIZE;

		int lineNr = y / DEFAULT_Y_STEP;
		int linePos = x / DEFAULT_CHAR_SIZE;
		return lineNr * DEFAULT_LINE_LENGTH + linePos ;

	}

	/**
	 * Get the X position on the Panel of a particular sequence position.
	 *
	 * @param structure index of the structure for the sequence position.
	 * @param pos sequence position, the aligned position index
	 * @return the point on a panel for a sequence position
	 */
	public Point getPanelPos(int structure, int pos) {
		Point p = new Point();

		int lineNr = pos / DEFAULT_LINE_LENGTH;
		int linePos = pos % DEFAULT_LINE_LENGTH;

		int x = linePos * DEFAULT_CHAR_SIZE + DEFAULT_X_SPACE +
				DEFAULT_LEGEND_SIZE;
		int y = lineNr * DEFAULT_Y_STEP + DEFAULT_Y_SPACE;

		y += DEFAULT_LINE_SEPARATION * structure;

		p.setLocation(x, y);
		return p;
	}

	/**
	 * Returns the index of the structure, for a given point in the Panel.
	 * Returns -1 if not over a position in the sequence alignment.
	 * @param point x and y coordinates in the panel
	 * @return which structure a point on the panel corresponds to
	 */
	public int getAligSeq(Point point) {

		for (int pos=0; pos<alignmentSize; pos++){
			int i = getSeqPos(pos, point);
			Point t = getPanelPos(pos,i);

			if ( Math.abs(t.x - point.x) <= DEFAULT_CHAR_SIZE &&
					Math.abs(t.y-point.y) < DEFAULT_CHAR_SIZE ) return pos;
		}
		return -1;
	}

	/**
	 * Provide the coordinates for where to draw the legend for
	 * line X given the structure index.
	 *
	 * @param lineNr line of the Panel
	 * @param structure the structure index
	 * @return get the point where to draw the legend
	 */
	public Point getLegendPosition(int lineNr, int structure) {

		int x = DEFAULT_X_SPACE ;
		int y = lineNr * DEFAULT_Y_STEP + DEFAULT_Y_SPACE;
		y += structure * DEFAULT_LINE_SEPARATION;

		Point p = new Point(x,y);
		return p;
	}

	public Point getEndLegendPosition(int lineNr, int structure) {

		int x = DEFAULT_LINE_LENGTH * DEFAULT_CHAR_SIZE + DEFAULT_X_SPACE + DEFAULT_LEGEND_SIZE + DEFAULT_RIGHT_SPACER;
		int y = lineNr * DEFAULT_Y_STEP + DEFAULT_Y_SPACE;
		y += structure * DEFAULT_LINE_SEPARATION;

		Point p = new Point(x,y);
		return p;
	}

}
