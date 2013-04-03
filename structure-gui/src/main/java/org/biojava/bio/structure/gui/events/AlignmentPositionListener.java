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
 * created at May 28, 2008
 */
package org.biojava.bio.structure.gui.events;

import org.biojava.bio.structure.gui.util.AlignedPosition;

public interface AlignmentPositionListener {

	public void mouseOverPosition(AlignedPosition p);
	public void positionSelected(AlignedPosition p);
	public void toggleSelection(AlignedPosition p);
	public void rangeSelected(AlignedPosition start , AlignedPosition end);
	public void selectionLocked();
	public void selectionUnlocked();
	
}
