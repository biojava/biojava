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
 * Created on Oct 6, 2009
 * Author: Andreas Prlic 
 *
 */

package org.biojava.bio.structure.align.gui.jmol;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;




import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.align.gui.JPrintPanel;
import org.biojava.bio.structure.align.gui.jmol.MyJmolStatusListener;
import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;

import org.jmol.api.JmolStatusListener;
import org.jmol.api.JmolViewer;


public class JmolPanel 
extends JPrintPanel
{
	private static final long serialVersionUID = -3661941083797644242L;

	private JmolViewer viewer;
	private JmolAdapter adapter;
	JmolStatusListener statusListener;
	final Dimension currentSize = new Dimension();
	final Rectangle rectClip = new Rectangle();

	public JmolPanel() {
		super();
		 statusListener = new MyJmolStatusListener();
		adapter = new SmarterJmolAdapter();
		viewer = JmolViewer.allocateViewer(this,
                adapter,
                null,null,null,null,
                 
                 statusListener);
    
	}

	public void paint(Graphics g) {
		
		getSize(currentSize);
		g.getClipBounds(rectClip);
		viewer.renderScreenImage(g, currentSize, rectClip);
	}

	public void evalString(String rasmolScript){

		viewer.evalString(rasmolScript);

	}

	public void openStringInline(String pdbFile){
		viewer.openStringInline(pdbFile);

	}
	public JmolViewer getViewer() {
		return viewer;
	}

	public JmolAdapter getAdapter(){
		return adapter;
	}
	
	public JmolStatusListener getStatusListener(){
		return statusListener;
	}
	public void executeCmd(String rasmolScript) {
		viewer.evalString(rasmolScript);
	}

	public void setStructure(Structure s)
	{
		String pdb = s.toPDB();
		viewer.openStringInline(pdb);
	}

}
