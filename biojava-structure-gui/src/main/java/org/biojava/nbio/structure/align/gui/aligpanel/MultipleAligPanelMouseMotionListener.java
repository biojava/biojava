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

import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;

/**
 * Mouse Motion Listener for the {@link MultipleAligPanel},
 * which provides methods to obtain positions of the mouse
 * and connect them to the sequence alignment positions using
 * the information in {@link MultipleAlignmentCoordManager}.
 *
 * @author Aleix Lafita
 *
 */
public class MultipleAligPanelMouseMotionListener
implements MouseMotionListener, MouseListener {

	private MultipleAligPanel parent;
	private List<AlignmentPositionListener> aligPosListeners;
	private int prevPos;

	private boolean isDragging ;
	private AlignedPosition selectionStart;
	private AlignedPosition selectionEnd;
	private boolean selectionLocked;

	public MultipleAligPanelMouseMotionListener(MultipleAligPanel parent){
		this.parent = parent;
		aligPosListeners = new ArrayList<AlignmentPositionListener>();
		prevPos = -1;
		isDragging = false;
		selectionStart = null;
		selectionEnd = null;
		selectionLocked = false;
	}

	public void addAligPosListener(AlignmentPositionListener li){
		aligPosListeners.add(li);
	}

	@Override
	public void mouseDragged(MouseEvent e) {

		AlignedPosition pos = getCurrentAlignedPosition(e);

		if (pos == null) return;
		if (prevPos == pos.getPos1() && isDragging) return;

		if (!isDragging) {
			isDragging = true;
			setSelectionLock(true);
		}

		if (selectionStart == null)	selectionStart = pos;
		if (selectionEnd == null) selectionEnd = pos;

		if (pos.getPos1() <= selectionStart.getPos1()) selectionStart = pos;
		else selectionEnd = pos;

		if (!keyPressed(e)) {
			triggerRangeSelected(selectionStart, selectionEnd);
		} else triggerRangeSelected(selectionStart, selectionEnd);

		prevPos = pos.getPos1();
	}

	private boolean keyPressed(MouseEvent e) {
		if (e.isShiftDown() || e.isControlDown() || e.isAltDown())
			return true;
		return false;
	}

	private void triggerRangeSelected(
			AlignedPosition start, AlignedPosition end) {
		for (AlignmentPositionListener li : aligPosListeners){
			li.rangeSelected(start, end);
		}
	}

	public void triggerSelectionLocked(boolean b) {
		selectionLocked = b;
		for (AlignmentPositionListener li : aligPosListeners){
			if (b) li.selectionLocked();
			else li.selectionUnlocked();
		}
	}

	@Override
	public void mouseMoved(MouseEvent e) {
		if ( selectionLocked) return;
		AlignedPosition pos = getCurrentAlignedPosition(e);
		if ( pos == null) return;

		triggerMouseOverPosition(pos);
	}

	private void triggerMouseOverPosition(AlignedPosition pos) {
		for (AlignmentPositionListener li : aligPosListeners)
			li.mouseOverPosition(pos);
	}

	private AlignedPosition getCurrentAlignedPosition(MouseEvent e) {

		MultipleAlignmentCoordManager coordManager = parent.getCoordManager();
		int aligSeq = coordManager.getAligSeq(e.getPoint());

		//We are not over a position in the sequences
		if (aligSeq == -1) return null;

		//Get sequence positions
		int seqPos = coordManager.getSeqPos(aligSeq, e.getPoint());
		if (seqPos < 0) return null;
		if (seqPos >= parent.length) return null;

		AlignedPosition pos = new AlignedPosition();
		pos.setPos1(seqPos);

		if (parent.getMapSeqToStruct().get(seqPos)!=-1)
			pos.setEquivalent(AlignedPosition.EQUIVALENT);

		return pos;
	}

	public void destroy() {
		aligPosListeners.clear();
		parent = null;
	}

	@Override
	public void mouseClicked(MouseEvent e) {}

	private void triggerToggleSelection(AlignedPosition pos) {
		for (AlignmentPositionListener li : aligPosListeners)
			li.toggleSelection(pos);
	}

	@Override
	public void mouseEntered(MouseEvent e) {}

	@Override
	public void mouseExited(MouseEvent e) {}

	@Override
	public void mousePressed(MouseEvent e) {

		selectionStart = null;
		selectionEnd = null;

		if (!keyPressed(e)) {

			setSelectionLock(false);
			triggerSelectionLocked(false);

			AlignedPosition pos = getCurrentAlignedPosition(e);
			if (pos != null) prevPos = pos.getPos1();
		}
	}

	private void setSelectionLock(boolean flag){
		selectionLocked = flag;
		triggerSelectionLocked(flag);
	}

	@Override
	public void mouseReleased(MouseEvent e) {

		isDragging = false;

		if (keyPressed(e)) {
			boolean keepOn = false;
			if (!selectionLocked) keepOn = true;
			setSelectionLock(true);

			// add to selection
			AlignedPosition pos = getCurrentAlignedPosition(e);
			if (pos == null) return;

			if (keepOn) triggerMouseOverPosition(pos);
			else triggerToggleSelection(pos);
			prevPos = pos.getPos1();
		}
	}
}
