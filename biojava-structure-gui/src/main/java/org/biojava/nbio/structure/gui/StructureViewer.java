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
package org.biojava.nbio.structure.gui;

import org.biojava.nbio.structure.Structure;

import java.awt.*;

public interface StructureViewer {

	public void setStructure(Structure structure);
	public void repaint();
	public void setSelection(Selection selection);
	public Selection getSelection();

	/** Apply this color to the current Selection
	 *
	 * @param red
	 */
	public void setColor(Color red);
	public Color getColor();


	/** Apply this style to the current selection
	 *
	 * @param wireframe renderstyle
	 */
	public void setStyle(RenderStyle wireframe);

	/** Clear the current display
	 *
	 *
	 */
	public void clear();

	/** Set the Zoom level
	 *
	 * @param i
	 */
	public void setZoom(int i);

}
