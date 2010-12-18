package org.biojava3.structure.gui;

import java.awt.Color;

import org.biojava.bio.structure.Structure;

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
