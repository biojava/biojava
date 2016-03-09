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

import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;

public class AligPanelMouseMotionListener implements MouseMotionListener, MouseListener {

	AligPanel parent;

	List<AlignmentPositionListener> aligPosListeners;
	int prevPos;

	boolean isDragging ;
	AlignedPosition selectionStart ;
	AlignedPosition selectionEnd;
	boolean selectionLocked;

	public AligPanelMouseMotionListener(AligPanel parent){
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

		if ( pos == null)
			return;


		int p = pos.getPos1();

		if ( prevPos == p && isDragging) {

			return;
		}


		if ( ! isDragging) {
			isDragging = true;

			setSelectionLock(true);

		}


		if ( selectionStart == null)
			selectionStart = pos;
		if ( selectionEnd == null)
			selectionEnd = pos;

		if ( p <= selectionStart.getPos1()) {
			//selectionEnd = selectionStart;
			selectionStart = pos;

		} else {
			selectionEnd = pos;
		}

		//System.out.println("sel start : " + selectionStart + " - " + selectionEnd);

		if ( ! keyPressed(e)) {
			triggerRangeSelected(selectionStart, selectionEnd);
		} else {
			triggerRangeSelected(selectionStart, selectionEnd);
			//triggerToggleRange(selectionStart, selectionEnd);
		}
		prevPos = p;
	}


	private boolean keyPressed(MouseEvent e) {
		if ( e.isShiftDown() || e.isControlDown() || e.isAltDown())
			return true;
		return false;
	}

	private void triggerRangeSelected(AlignedPosition start,
			AlignedPosition end) {
		for (AlignmentPositionListener li : aligPosListeners){
			li.rangeSelected(start, end);
		}
	}
	public void triggerSelectionLocked(boolean b) {
		selectionLocked = b;
		for (AlignmentPositionListener li : aligPosListeners){
			if ( b)
				li.selectionLocked();
			else
				li.selectionUnlocked();
		}

	}
	@Override
	public void mouseMoved(MouseEvent e) {
		if ( selectionLocked)
			return;
		AlignedPosition pos = getCurrentAlignedPosition(e);
		if ( pos == null)
			return;

		triggerMouseOverPosition(pos);




	}

	private void triggerMouseOverPosition(AlignedPosition pos) {
		for (AlignmentPositionListener li : aligPosListeners){
			li.mouseOverPosition(pos);
		}


	}

	private AlignedPosition getCurrentAlignedPosition(MouseEvent e){
		AFPChainCoordManager coordManager = parent.getCoordManager();

		int aligSeq = coordManager.getAligSeq(e.getPoint());

		// we are over a position in the sequences
		if ( aligSeq == -1) {
			return null;
		}

		//get sequence positions
		int seqPos = coordManager.getSeqPos(aligSeq, e.getPoint());

		//if ( prevPos == seqPos)
		//	return null;


		//prevPos = seqPos;

		if ( seqPos < 0)
			return null;



		AFPChain afpChain = parent.getAFPChain();
		char[] aligs1  = afpChain.getAlnseq1();
		char[] aligs2  = afpChain.getAlnseq2();

		if ( seqPos >= afpChain.getAlnLength()) {
			//System.err.println("seqpos " + seqPos +" >= " + afpChain.getAlnLength());
			return null;
		}

		//System.out.println("alignment " + aligSeq + " " + seqPos + " : ");
		AlignedPosition pos = new AlignedPosition();
		pos.setPos1(seqPos);
		pos.setPos2(seqPos);

		if ( aligs1[seqPos] != '-' && aligs2[seqPos] != '-'){
			pos.setEquivalent(AlignedPosition.EQUIVALENT);
		}

		return pos;
	}

	public void destroy() {
		aligPosListeners.clear();
		parent = null;

	}

	@Override
	public void mouseClicked(MouseEvent e) {

	}

	private void triggerToggleSelection(AlignedPosition pos) {
		for (AlignmentPositionListener li : aligPosListeners){
			li.toggleSelection(pos);
		}

	}


//	private void triggerToggleRange(AlignedPosition start,
//			AlignedPosition end) {
//		for (AlignmentPositionListener li : aligPosListeners){
//			for ( int i = start.getPos1() ; i < end.getPos1() ; i++){
//				AlignedPosition pos = new AlignedPosition();
//				pos.setPos1(i);
//				pos.setPos2(i);
//				li.toggleSelection(pos);
//			}
//		}
//
//	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseExited(MouseEvent e) {


	}

	@Override
	public void mousePressed(MouseEvent e) {

		selectionStart = null;
		selectionEnd = null;

		if ( ! keyPressed(e) ) {
			//System.out.println("mouse pressed " + e.isShiftDown() + " selection locked: " + selectionLocked);

			setSelectionLock(false);
			//System.out.println("selection unlocked by mousePressed");
			triggerSelectionLocked(false);

			AlignedPosition pos = getCurrentAlignedPosition(e);
			if ( pos != null) {
				prevPos = pos.getPos1();
			}


		}


	}


	private void setSelectionLock(boolean flag){
		selectionLocked = flag;
		triggerSelectionLocked(flag);
	}

	@Override
	public void mouseReleased(MouseEvent e) {

		isDragging = false;
		//System.out.println("mouse released... " + e.isShiftDown() + " selection locked:" + selectionLocked);
		if ( keyPressed(e)) {
			boolean keepOn = false;
			if ( ! selectionLocked)
				keepOn = true;
			setSelectionLock(true);


			// add to selection
			AlignedPosition pos = getCurrentAlignedPosition(e);
			if ( pos == null)
				return;

			if ( keepOn)
				triggerMouseOverPosition(pos);
			else
				triggerToggleSelection(pos);
			prevPos = pos.getPos1() ;

		}




	}

}
