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
 * Created on Nov 29, 2005
 *
 */
package org.biojava.nbio.structure.gui.util;

//import java.awt.Cursor;

import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.gui.SequenceDisplay;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;


/** a mouse listener for the AbstractChainRenderer class
 * it listens to all mouse events and triggeres appropriate
 * SequenceListener and FeatureListener events
 *
 * @author Andreas Prlic
 *
 */
public class SequenceMouseListener implements
MouseListener,
MouseMotionListener
{


	boolean selectionLocked;
	boolean dragging;

	int selectionStart;
	int selectionEnd;
	int draggingStart;
	int oldSelectionStart;

	int chainLength;


	CoordManager coordManager;

	//private static final Logger logger = LoggerFactory.getLogger(SequenceMouseListener.class);

	SequenceDisplay parent;

	List<AlignmentPositionListener> alignmentPositionListeners;

	public SequenceMouseListener(SequenceDisplay parent) {
		super();

		this.parent = parent;

		selectionLocked = false;
		dragging = false;

		selectionStart = -1 ;
		selectionEnd = -1;
		oldSelectionStart = -1;
		draggingStart = -1;
		chainLength = 0;
		//clearSequenceListeners();

		coordManager = new CoordManager();

		alignmentPositionListeners = new ArrayList<AlignmentPositionListener>();
		//renderer.getLayeredPane().addMouseListener(popupFrame);

	}

	public void clearListeners(){

		alignmentPositionListeners.clear();
	}
	public void addAlignmentPositionListener(AlignmentPositionListener li){
		alignmentPositionListeners.add(li);
	}

	@Override
	public void mousePressed(MouseEvent event) {

		int pos  = getSeqPos(event);

		draggingStart=pos;
		selectionStart = pos ;
		//selectionEnd   = pos ;
		//triggerClearSelection();
		// triggerSelectionLocked(false);
		triggerMouseOverPosition(pos,event.getY());




	}

	// mouse motion part

	/** get the sequence position of the current mouse event
	 *
	 */
	private int getSeqPos(MouseEvent e) {

		int x = e.getX();
		//int y = e.getY();
		//float scale = seqScale.getScale();
		//int DEFAULT_X_START = SequenceScalePanel.DEFAULT_X_START;

		float scale = parent.getScale();

		coordManager.setScale(scale);
		int seqpos = coordManager.getSeqPos(x-2);


		return seqpos  ;
	}



	public void setChain(Chain c){

		chainLength = c.getAtomLength();
		coordManager.setLength(chainLength);
	}

	private void setSelectionStart(int start){
		if ( start < 0 )
			start = 0;
		if ( start > chainLength)
			start = chainLength;
		selectionStart = start;
	}


	private void setSelectionEnd(int end){
		//if ( end < 0 )
			//end = 0;
		if ( end > chainLength)
			end = chainLength;
		selectionEnd = end;
	}


	@Override
	public void mouseDragged(MouseEvent e) {
		dragging = true;

		int pos = getSeqPos(e) ;


		if (( pos < 0)|| ( pos> chainLength)){
			return;
		}

		if (pos == oldSelectionStart){
			return;
		}
		oldSelectionStart = pos;

		if ( pos > draggingStart){
			selectionStart = draggingStart;
			selectionEnd = pos ;
		} else {
			selectionStart = pos;
			selectionEnd = draggingStart;
		}
		//triggerNewSequenceRange(selectionStart,selectionEnd);

	}




	@Override
	public void mouseMoved(MouseEvent e) {


		if ( selectionLocked )
			return;



		int pos = getSeqPos(e) ;
		//System.out.println("mouse moved " + pos);



		oldSelectionStart = pos;


		this.setSelectionStart(pos);
		this.setSelectionEnd(pos);

		triggerMouseOverPosition(pos,e.getY());



	}



	@Override
	public void mouseClicked(MouseEvent arg0) {


	}

	@Override
	public void mouseEntered(MouseEvent arg0) {


	}

	@Override
	public void mouseExited(MouseEvent arg0) {

	}


	@Override
	public void mouseReleased(MouseEvent event) {
		// logger.info("mouse released");



		draggingStart = -1;


		if ( dragging) {
			if  ( ! selectionLocked) {
				// triggerSelectionLocked(true);
			}
		} else {}

		dragging = false ;
	}



	protected void triggerMouseOverPosition(int pos, int mouseY){
		if ( selectionLocked )
			return;
		List<AlignedPosition> apos = parent.getAligMap();
		if (pos > apos.size()-1)
			return;
		//System.out.println(parent.getAligMap().get(pos));

		for (AlignmentPositionListener li : alignmentPositionListeners) {
			li.mouseOverPosition(apos.get(pos));
		}



	}

}




