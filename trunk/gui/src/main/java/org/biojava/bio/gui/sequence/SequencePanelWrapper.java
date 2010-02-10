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
package org.biojava.bio.gui.sequence;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JPanel;

import org.biojava.bio.gui.sequence.tracklayout.SimpleTrackLayout;
import org.biojava.bio.gui.sequence.tracklayout.TrackLayout;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;

/**
 * Handles multiple SequencePanels and Ranges so that a Sequence can be wrapped
 * over more than one line on screen. This is particularly useful for viewing
 * Protein sequences that would be viewed at a single residue resolution.
 * 
 * The interface is very similar to that of the SequencePanels that it wraps
 * 
 * @author Mark Southern
 * @see SequencePanel
 * @since 1.5
 */
public class SequencePanelWrapper extends JPanel {

	/**
	 * Generated Serial Version UID
	 */
	private static final long serialVersionUID = 8749249181471157230L;
	protected SequencePanel[] seqPanels = new SequencePanel[0];
	private RangeLocation range;
	private Sequence sequence;
	private SequenceRenderer renderer;
	private double scale = 14.0;
	private java.awt.RenderingHints hints;
	private int direction = SequencePanel.HORIZONTAL;
	private TrackLayout trackLayout = new SimpleTrackLayout();
	private List<SequenceViewerListener> viewerListeners = new ArrayList<SequenceViewerListener>();
	private List<SequenceViewerMotionListener> motionListeners = new ArrayList<SequenceViewerMotionListener>();

	/**
	 * Creates a new instance of WrappedSequencePanel
	 */
	public SequencePanelWrapper() {
		initComponents();
	}

	protected void initComponents() {
		setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		setBackground(Color.WHITE);
	}

	/*
	 * Sets the track (line) layout strategy. Layouts include a simple wrapping
	 * at a given number of residues and user-defined layouts that can contain
	 * arbitrary length rows.
	 */
	public void setTrackLayout(TrackLayout tl) {
		this.trackLayout = tl;
		trackLayout.setSequence(getSequence());
		trackLayout.setRange(getRange());

		if (isActive()) {
			refreshSequencePanels();
		}
	}

	public TrackLayout getTrackLayout() {
		return this.trackLayout;
	}

	protected boolean isActive() {
		return (sequence != null) && (renderer != null) && (range != null);
	}

	public void setScale(double scale) {
		this.scale = scale;

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].setScale(scale);
		}
	}

	public double getScale() {
		return scale;
	}

	public void setDirection(int direction) {
		this.direction = direction;

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].setDirection(direction);
		}

		if (isActive()) {
			refreshSequencePanels();
		}
	}

	public int getDirection() {
		return direction;
	}

	/*
	 * a convenience method. The wrap is passed through to the TrackLayout
	 * implementation.
	 */
	public synchronized void setWrap(int w) {
		trackLayout.setWrap(w);

		if (!isActive()) {
			return;
		}

		refreshSequencePanels();
	}

	public int getWrap() {
		return trackLayout.getWrap();
	}

	public void setRenderingHints(java.awt.RenderingHints hints) {
		this.hints = hints;

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].setRenderingHints(hints);
		}
	}

	public java.awt.RenderingHints getRenderingHints() {
		return hints;
	}

	public void setRenderer(SequenceRenderer renderer) {
		this.renderer = renderer;

		for (int i = 0; i < seqPanels.length; i++) {
			try {
				seqPanels[i].setRenderer(renderer);
			} catch (ChangeVetoException e) {
				// should never get here
				e.printStackTrace();
			}
		}
	}

	public SequenceRenderer getRenderer() {
		return renderer;
	}

	public void setSequence(Sequence seq) {
		if (seq == null) {
			removeSeqPanels();

			return;
		}

		trackLayout.setSequence(sequence = seq);
		trackLayout.setRange(range = new RangeLocation(1, seq.length()));
		refreshSequencePanels();
	}

	public Sequence getSequence() {
		return sequence;
	}

	public void setRange(RangeLocation loc) {
		trackLayout.setRange(range = loc);
		refreshSequencePanels();
	}

	public RangeLocation getRange() {
		return range;
	}

	private void removeSeqPanels() {
		for (int i = 0; i < seqPanels.length; i++) {
			remove(seqPanels[i]);
		}

		setSize(0, 0); // make sure view is refreshed properly
		setLayout(null);
		revalidate();
	}

	public void resizeAndValidate() {
		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].resizeAndValidate();
		}
	}

	protected void refreshSequencePanels() {
		removeSeqPanels();

		RangeLocation[] ranges = trackLayout.getRanges();
		seqPanels = new SequencePanel[ranges.length];

		if (getDirection() == SequencePanel.HORIZONTAL) {
			setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		} else {
			setLayout(new BoxLayout(this, BoxLayout.X_AXIS));
		}

		for (int i = 0; i < ranges.length; i++) {
			// logger.debug("Setting sequence panel " + i);
			seqPanels[i] = new SequencePanel();
			seqPanels[i].setFont(getFont());
			seqPanels[i].setAlignmentX(LEFT_ALIGNMENT);
			seqPanels[i].setAlignmentY(TOP_ALIGNMENT);
			seqPanels[i].setSequence(getSequence());

			// seqPanels[i].setRange( ranges[i] );
			// bug in biojava-1.4pre1 - add +1 to max range
			seqPanels[i].setRange(new RangeLocation(ranges[i].getMin(),
					ranges[i].getMax() + 1));
			seqPanels[i].setDirection(getDirection());
			seqPanels[i].setScale(getScale());

			for (Iterator<SequenceViewerListener> it = viewerListeners
					.iterator(); it.hasNext();) {
				seqPanels[i].addSequenceViewerListener(it.next());
			}

			for (Iterator<SequenceViewerMotionListener> jt = motionListeners
					.iterator(); jt.hasNext();) {
				seqPanels[i].addSequenceViewerMotionListener(jt.next());
			}

			add(seqPanels[i]);
		}

		setRenderer(getRenderer());
	}

	public void addSequenceViewerListener(SequenceViewerListener l) {
		viewerListeners.add(l);

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].addSequenceViewerListener(l);
		}
	}

	public void removeSequenceViewerListener(SequenceViewerListener l) {
		viewerListeners.remove(l);

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].removeSequenceViewerListener(l);
		}
	}

	public void addSequenceViewerMotionListener(SequenceViewerMotionListener l) {
		motionListeners.add(l);

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].addSequenceViewerMotionListener(l);
		}
	}

	public void removeSequenceViewerMotionListener(
			SequenceViewerMotionListener l) {
		motionListeners.remove(l);

		for (int i = 0; i < seqPanels.length; i++) {
			seqPanels[i].removeSequenceViewerMotionListener(l);
		}
	}

	public void setFont(Font f) {
		try {
			super.setFont(f);

			for (int i = 0; i < seqPanels.length; i++) {
				seqPanels[i].setFont(f);
				seqPanels[i].resizeAndValidate();
			}
		} catch (NullPointerException e) {
		}
	}

	public void paint(Graphics g) {
		// some wierd sequence graphics exception is going on in biojava
		try {
			super.paint(g);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Caught sequence graphics exception in paint() "
					+ e.toString());
			repaint();
		}
	}

	public void paintComponent(Graphics g) {
		// some wierd sequence graphics exception is going on in biojava-1.3
		try {
			super.paintComponent(g);
		} catch (Exception e) {
			e.printStackTrace();
			System.out
					.println("Caught sequence graphics exception in paintComponent() "
							+ e.toString());
			repaint();
		}
	}
} // class
