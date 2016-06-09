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

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.JPrintPanel;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.MultipleAlignmentJmolDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.multiple.MultipleAlignment;
import org.biojava.nbio.structure.align.multiple.MultipleAlignmentEnsembleImpl;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentTools;
import org.biojava.nbio.structure.align.multiple.util.MultipleAlignmentWriter;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;

/**
 * A JPanel that can display the sequence alignment of a
 * {@link MultipleAlignment} in a nice way and interact with Jmol by
 * selecting the aligned atoms of the  sequence selection.
 * <p>
 * Coloring options include: sequence similarity, by Block or by Structure.
 * Colors are connected with the JmolPanel using the same pattelete.
 * <p>
 * The positions can be selected individually or in ranges and they will be
 * translated to jmol commands.
 *
 * @author Aleix Lafita
 * @since 4.1.0
 *
 */
public class MultipleAligPanel extends JPrintPanel
implements AlignmentPositionListener, WindowListener {

	private static final long serialVersionUID = -6892229111166263764L;

	private MultipleAlignment multAln;
	private List<String> alnSeq; //sequence alignment
	private List<Integer> mapSeqToStruct;  //mapping from sequence to structure

	int size; 			//number of structures
	int length; 		//number of aligned positions in sequence alignment

	private Font seqFont;
	private Font eqFont;

	private AbstractAlignmentJmol jmol;
	private MultipleAligPanelMouseMotionListener mouseMoLi;
	private MultipleAlignmentCoordManager coordManager;

	private BitSet selection;
	private boolean selectionLocked;

	private boolean colorBySimilarity=false;
	private boolean colorByAlignmentBlock=false;

	private static final Color COLOR_EQUAL   = Color.decode("#6A93D4");
	private static final Color COLOR_SIMILAR = Color.decode("#D460CF");

	/**
	 * Default constructor. Empty MultipleAligPanel instance.
	 */
	public MultipleAligPanel(){
		super();
		this.setBackground(Color.white);
		seqFont = new Font("SansSerif",Font.PLAIN,12);
		eqFont = new Font("SansSerif",Font.BOLD,12);

		mouseMoLi = new MultipleAligPanelMouseMotionListener(this);
		this.addMouseMotionListener(mouseMoLi);
		this.addMouseListener(mouseMoLi);
		mouseMoLi.addAligPosListener(this);
		selection = new BitSet();

		multAln = null;
		alnSeq = null;
		mapSeqToStruct = null;
	}

	/**
	 * Constructor using an afpChain and the atom arrays for pairwise
	 * alignments. The AFPChain is converted into a MultipleAlignment.
	 *
	 * @param afpChain
	 * @param ca1
	 * @param ca2
	 * @throws StructureException
	 */
	public MultipleAligPanel(AFPChain afpChain, Atom[] ca1, Atom[] ca2,
			AbstractAlignmentJmol jmol) throws StructureException {

		this();

		String algorithm = afpChain.getAlgorithmName();
		boolean flex = false;
		if (algorithm != null){
			if (algorithm.contains("flexible")) flex = true;
		}

		//Convert the apfChain into a MultipleAlignment object
		MultipleAlignmentEnsembleImpl ensemble =
				new MultipleAlignmentEnsembleImpl(afpChain, ca1, ca2, flex);
		this.multAln = ensemble.getMultipleAlignment(0);

		//Create the sequence alignment and the structure-sequence mapping.
		this.mapSeqToStruct = new ArrayList<Integer>();
		this.alnSeq = MultipleAlignmentTools.getSequenceAlignment(
				this.multAln, this.mapSeqToStruct);

		//Initialize other memeber variables of the panel
		this.size = multAln.size();
		this.length = alnSeq.get(0).length();

		coordManager = new MultipleAlignmentCoordManager(size, length);
		this.jmol = jmol;
	}

	/**
	 * Constructor using a MultipleAlignment.
	 *
	 * @param multAln
	 * @param colors
	 */
	public MultipleAligPanel(MultipleAlignment msa, AbstractAlignmentJmol jm) {
		this();
		this.multAln = msa;

		//Create the sequence alignment and the structure-sequence mapping.
		this.mapSeqToStruct = new ArrayList<Integer>();
		this.alnSeq = MultipleAlignmentTools.getSequenceAlignment(
				this.multAln, this.mapSeqToStruct);

		this.size = multAln.size();
		this.length = this.alnSeq.get(0).length();

		coordManager = new MultipleAlignmentCoordManager(size, length);
		this.jmol = jm;
	}

	public MultipleAlignmentCoordManager getCoordManager() {
		return coordManager;
	}

	public void addAlignmentPositionListener(AlignmentPositionListener li){
		mouseMoLi.addAligPosListener(li);
	}

	public void destroy(){

		multAln = null;
		alnSeq = null;
		mouseMoLi.destroy();
		jmol = null;
		selection = null;
	}

	@Override
	public void paintComponent(Graphics g){

		super.paintComponent(g);

		Graphics2D g2D = (Graphics2D) g;

		g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
				RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);

		int startpos = 0;
		int endpos = length;

		String summary = multAln.toString();
		g2D.drawString(summary, 20, coordManager.getSummaryPos());

		Color significantCol = Color.black;
		//if (multAln.isSignificantResult()) significantCol = Color.green;

		g2D.setPaint(significantCol);
		Rectangle sig = new Rectangle(10,10,10,10);
		g2D.fill(sig);

		for (int i = startpos; i < endpos; i++){

			boolean isGapped = false;
			g2D.setFont(seqFont);

			if (mapSeqToStruct.get(i)!=-1) g2D.setFont(eqFont);
			else isGapped = true;

			//Loop through every structure to get all the points
			List<Point> points = new ArrayList<Point>();
			for (int str=0; str<size; str++) points.add(
					coordManager.getPanelPos(str,i));
			Point p1 = points.get(0);
			Point p2 = points.get(points.size()-1);

			for (int str=0; str<size; str++){

				char c = alnSeq.get(str).charAt(i);
				Color bg = jmol.getColorPalette().getColorPalette(size)[str];

				//Color only if the position is aligned
				if (!isGapped){
					//Color by sequence similarity (equal or similar)
					if (colorBySimilarity){
						boolean equal = true;
						boolean similar = true;
						char c1 = '-';
						for (int st=0; st<size-1; st++){
							if (alnSeq.get(st).charAt(i) != '-') {
								c1 = alnSeq.get(st).charAt(i);
							}
							char c2 = alnSeq.get(st+1).charAt(i);
							//If any position is a gap continue
							if (c1=='-' || c2=='-' ||
									Character.isLowerCase(c1) ||
									Character.isLowerCase(c2)) {
								continue;
							}
							if (equal && c1 == c2)
								continue;
							else equal = false;
							if (AFPAlignmentDisplay.aaScore(c1, c2) > 0)
								continue;
							else similar = false; break;
						}
						if (equal) bg = COLOR_EQUAL;
						else if (similar) bg = COLOR_SIMILAR;
						else bg = Color.LIGHT_GRAY;
					}
					//Color by alignment block the same way as in the Jmol
					else if (colorByAlignmentBlock){
						int blockNr = MultipleAlignmentTools.
								getBlockForSequencePosition(
										multAln,mapSeqToStruct,i);
						bg = jmol.getColorPalette().getColorPalette(
								multAln.getBlocks().size())[blockNr];
					}
					if (isSelected(i)) bg = Color.YELLOW;

					if (Character.isUpperCase(c) && c!='-'){
						g2D.setPaint(bg);
						Rectangle rec = new Rectangle(points.get(str).x-1,
								points.get(str).y-11, (p2.x-p1.x)+12,
								(p2.y-p1.y)/size);
						g2D.fill(rec);
					}
				}

				// draw the AA sequence
				g2D.setColor(Color.black);
				g2D.drawString(String.valueOf(c), points.get(str).x, points.get(str).y);
			}
		}

		int nrLines = (length-1) /
				(MultipleAlignmentCoordManager.DEFAULT_LINE_LENGTH);

		for (int i = 0 ; i < nrLines+1 ; i++){

			// draw legend at i
			for (int str=0; str<size; str++){

				Point p1 = coordManager.getLegendPosition(i,str);

				int aligPos = i *
						MultipleAlignmentCoordManager.DEFAULT_LINE_LENGTH;
				Atom a1 = null;
				while (a1==null &&
						aligPos < Math.min((i+1)*MultipleAlignmentCoordManager.
								DEFAULT_LINE_LENGTH-1,length)){
					a1 = MultipleAlignmentTools.getAtomForSequencePosition(
							multAln, mapSeqToStruct, str, aligPos);
					aligPos++;
				}
				String label1 = JmolTools.getPdbInfo(a1,false);
				g2D.drawString(label1, p1.x,p1.y);

				Point p3 = coordManager.getEndLegendPosition(i,str);

				aligPos = (i*MultipleAlignmentCoordManager.DEFAULT_LINE_LENGTH+
						MultipleAlignmentCoordManager.DEFAULT_LINE_LENGTH - 1);
				if (aligPos > length) aligPos = length-1;
				Atom a3 = null;
				while (a3==null && aligPos > Math.max(i*
						MultipleAlignmentCoordManager.DEFAULT_LINE_LENGTH,0)){
					a3 = MultipleAlignmentTools.getAtomForSequencePosition(
							multAln, mapSeqToStruct, str, aligPos);
					aligPos--;
				}

				String label3 = JmolTools.getPdbInfo(a3,false);

				g2D.drawString(label3, p3.x,p3.y);
			}
		}
	}

	private boolean isSelected(int alignmentPosition) {
		return selection.get(alignmentPosition);
	}

	@Override
	public void mouseOverPosition(AlignedPosition p) {

		if (!selectionLocked) selection.clear();

		selection.set(p.getPos1());
		updateJmolDisplay();
		this.repaint();
	}

	private void updateJmolDisplay() {

		if (jmol == null) return;

		StringBuffer cmd = new StringBuffer("select ");
		int nrSelected = 0;
		for (int i=0; i<length; i++){
			if (selection.get(i)){
				for (int str=0; str<size; str++){
					Atom a = MultipleAlignmentTools.getAtomForSequencePosition(
							multAln, mapSeqToStruct,str,i);
					if (a != null) {
						cmd.append(JmolTools.getPdbInfo(a));
						cmd.append("/"+(str+1)+", ");
					}
				}
				nrSelected++;
			}
		}
		if (nrSelected == 0) cmd.append(" none;");
		else cmd.append(" none; set display selected;");
		//System.out.println(cmd.toString());
		jmol.evalString(cmd.toString());
	}


	@Override
	public void positionSelected(AlignedPosition p) {
		mouseOverPosition(p);
	}

	@Override
	public void rangeSelected(AlignedPosition start, AlignedPosition end) {

		if (!selectionLocked) selection.clear();
		selection.set(start.getPos1(), end.getPos1()+1);
		updateJmolDisplay();
		this.repaint();
	}

	@Override
	public void selectionLocked() {
		selectionLocked = true;
	}

	@Override
	public void selectionUnlocked() {
		selectionLocked = false;
		selection.clear();
		this.repaint();
	}

	@Override
	public void toggleSelection(AlignedPosition p) {
		selection.flip(p.getPos1());
		updateJmolDisplay();
		this.repaint();
	}

	public void setStructureAlignmentJmol(AbstractAlignmentJmol jmol) {
		this.jmol = jmol;

	}

	@Override
	public void windowActivated(WindowEvent e) {}

	@Override
	public void windowClosed(WindowEvent e) {}

	@Override
	public void windowClosing(WindowEvent e) {
		destroy();
	}

	@Override
	public void windowDeactivated(WindowEvent e) {}

	@Override
	public void windowDeiconified(WindowEvent e) {}

	@Override
	public void windowIconified(WindowEvent e) {}

	@Override
	public void windowOpened(WindowEvent e) {}

	@Override
	public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();
		if ( cmd.equals(MenuCreator.PRINT)) {
			super.actionPerformed(e);
		} else if (cmd.equals(MenuCreator.FASTA_FORMAT)){
			String result = MultipleAlignmentWriter.toFASTA(multAln);
			MultipleAlignmentJmolDisplay.showAlignmentImage(multAln, result);
		} else if ( cmd.equals(MenuCreator.PAIRS_ONLY)) {
			String result = MultipleAlignmentWriter.toAlignedResidues(multAln);
			MultipleAlignmentJmolDisplay.showAlignmentImage(multAln, result);
		} else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
			String result = MultipleAlignmentWriter.toFatCat(multAln);
			MultipleAlignmentJmolDisplay.showAlignmentImage(multAln, result);
		} else if (cmd.equals(MenuCreator.SELECT_EQR)){
			selectEQR();
		} else if ( cmd.equals(MenuCreator.SIMILARITY_COLOR)){
			colorBySimilarity(true);
		} else if (cmd.equals(MenuCreator.EQR_COLOR)){
			colorBySimilarity(false);
		} else if ( cmd.equals(MenuCreator.FATCAT_BLOCK)){
			colorByAlignmentBlock();
		} else {
			System.err.println("Unknown command:" + cmd);
		}
	}

	private void selectEQR() {

		selection.clear();

		for (int pos=0; pos<length; pos++){
			if (mapSeqToStruct.get(pos)!=-1) selection.flip(pos);
		}
		mouseMoLi.triggerSelectionLocked(true);
		updateJmolDisplay();
		this.repaint();
	}

	private void colorByAlignmentBlock() {
		colorByAlignmentBlock = true;
		colorBySimilarity = false;
		this.repaint();
	}

	private void colorBySimilarity(boolean flag) {
		this.colorBySimilarity = flag;
		colorByAlignmentBlock = false;
		this.repaint();
	}

	public List<Atom[]> getAtomArrays() {
		return multAln.getAtomArrays();
	}
	public MultipleAlignment getMultipleAlignment(){
		return multAln;
	}
	public List<String> getAlnSequences(){
		return alnSeq;
	}
	public List<Integer> getMapSeqToStruct(){
		return mapSeqToStruct;
	}
}
