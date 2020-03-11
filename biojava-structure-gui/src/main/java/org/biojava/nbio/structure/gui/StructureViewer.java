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

	void setStructure(Structure structure);
	void repaint();
	void setSelection(Selection selection);
	Selection getSelection();

	/** Apply this color to the current Selection
	 *
	 * @param red
	 */
    void setColor(Color red);
	Color getColor();


	/** Apply this style to the current selection
	 *
	 * @param wireframe renderstyle
	 */
    void setStyle(RenderStyle wireframe);

	/** Clear the current display
	 *
	 *
	 */
    void clear();

	/** Set the Zoom level
	 *
	 * @param i
	 */
    void setZoom(int i);

}
