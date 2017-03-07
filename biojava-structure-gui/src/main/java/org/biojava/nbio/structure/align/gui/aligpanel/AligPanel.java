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

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.align.gui.DisplayAFP;
import org.biojava.nbio.structure.align.gui.JPrintPanel;
import org.biojava.nbio.structure.align.gui.MenuCreator;
import org.biojava.nbio.structure.align.gui.StructureAlignmentDisplay;
import org.biojava.nbio.structure.align.gui.jmol.AbstractAlignmentJmol;
import org.biojava.nbio.structure.align.gui.jmol.JmolTools;
import org.biojava.nbio.structure.align.gui.jmol.StructureAlignmentJmol;
import org.biojava.nbio.structure.align.model.AFPChain;
import org.biojava.nbio.structure.align.model.AfpChainWriter;
import org.biojava.nbio.structure.align.util.AFPAlignmentDisplay;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.align.util.UserConfiguration;
import org.biojava.nbio.structure.align.webstart.WebStartMain;
import org.biojava.nbio.structure.align.xml.AFPChainXMLParser;
import org.biojava.nbio.structure.gui.events.AlignmentPositionListener;
import org.biojava.nbio.structure.gui.util.AlignedPosition;
import org.biojava.nbio.structure.gui.util.color.ColorUtils;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.BitSet;
import java.util.List;


/** A JPanel that can display an AFPChain in a nice way and interact with Jmol.
 *
 * @author Andreas Prlic
 *
 */
public class AligPanel  extends JPrintPanel implements AlignmentPositionListener, WindowListener{

	/**
	 *
	 */
	private static final long serialVersionUID = -6892229111166263764L;

	private AFPChain afpChain;
	private AFPChainCoordManager coordManager ;
	private Font seqFont;
	private Font eqFont;
	private AbstractAlignmentJmol jmol;
	private AligPanelMouseMotionListener mouseMoLi;

	private BitSet selection;

	private boolean selectionLocked;
	private Atom[] ca1;
	private Atom[] ca2;

	private boolean colorBySimilarity;

	private boolean colorByAlignmentBlock;

	private static final Color COLOR_EQUAL   = Color.decode("#6A93D4");
	private static final Color COLOR_SIMILAR = Color.decode("#D460CF");


	public static void main(String[] args){

		String file = "/Users/ap3/tmp/4hhb.ce";

		try {
			BufferedReader in = new BufferedReader(new FileReader(file));
			StringBuffer xml = new StringBuffer();
			String str;
			while ((str = in.readLine()) != null) {
				xml.append(str);
			}
			in.close();

			AFPChain[] afps = AFPChainXMLParser.parseMultiXML(xml.toString());
			AFPChain afpChain = afps[0];

			UserConfiguration config = WebStartMain.getWebStartConfig();
			AtomCache cache = new AtomCache(config.getPdbFilePath(),config.getCacheFilePath());

			Atom[] ca1 = cache.getAtoms(afpChain.getName1());
			Atom[] ca2 = cache.getAtoms(afpChain.getName2());

			AFPChainXMLParser.rebuildAFPChain(afpChain, ca1, ca2);


			//StructureAlignment algorithm = StructureAlignmentFactory.getAlgorithm(afpChain.getAlgorithmName());
			StructureAlignmentJmol jmol= StructureAlignmentDisplay.display(afpChain, ca1, ca2);

			DisplayAFP.showAlignmentPanel(afpChain, ca1, ca2, jmol);




		} catch (Exception e){
			e.printStackTrace();
		}
	}


	public AligPanel(){
		super();
		this.setBackground(Color.white);
		coordManager = new AFPChainCoordManager();
		seqFont = new Font("SansSerif",Font.PLAIN,12);
		eqFont = new Font("SansSerif",Font.BOLD,12);

		mouseMoLi = new AligPanelMouseMotionListener(this);
		this.addMouseMotionListener(mouseMoLi);
		this.addMouseListener(mouseMoLi);
		mouseMoLi.addAligPosListener(this);

		selection = new BitSet();
		colorBySimilarity = false;
		colorByAlignmentBlock = false;
	}



	public AFPChainCoordManager getCoordManager() {
		return coordManager;
	}


	public void addAlignmentPositionListener(AlignmentPositionListener li){
		mouseMoLi.addAligPosListener(li);
	}

	public void destroy(){

		setAFPChain(null);
		mouseMoLi.destroy();
		jmol = null;
		ca1 = null;
		ca2 = null;
		selection = null;
	}

	public AFPChain getAFPChain(){
		return afpChain;
	}

	public void setAFPChain(AFPChain afpChain) {

		this.afpChain = afpChain;
		coordManager.setAFPChain(afpChain);
		if ( afpChain != null) {
			selection = new BitSet (afpChain.getAlnLength());
			if ( afpChain.getBlockNum() > 1) {
				colorByAlignmentBlock = true;
			}
		}

	}


@Override
public void paintComponent(Graphics g){

		super.paintComponent(g);

		Graphics2D g2D = (Graphics2D) g;


		g2D.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
				RenderingHints.VALUE_TEXT_ANTIALIAS_ON);

		g2D.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);


		// only draw within the ranges of the Clip
		//Rectangle drawHere = g2D.getClipBounds();


		//int startpos = coordManager.getSeqPos(0,drawHere.x);
		//int endpos   = coordManager.getSeqPos(0,drawHere.x+drawHere.width-2);


		char[] seq1 = afpChain.getAlnseq1();
		char[] seq2 = afpChain.getAlnseq2();
		char[] symb = afpChain.getAlnsymb();

		int startpos = 0;
		int endpos = afpChain.getAlnLength();

		String summary = afpChain.toString();
		g2D.drawString(summary, 20, coordManager.getSummaryPos());

		Color significantCol = Color.red;
		if ( afpChain.isSignificantResult())
			significantCol = Color.green;

		g2D.setPaint(significantCol);
		// draw a darker backgroun
		Rectangle sig = new Rectangle(10,10,10,10);
		g2D.fill(sig);
		boolean isFATCAT = false;
		if ( afpChain.getAlgorithmName().startsWith("jFatCat")){
			isFATCAT = true;
		}
		for ( int i = startpos ; ((i <= endpos) && ( i < afpChain.getAlnLength())) ;i++){

			// TODO:
			// color amino acids by hydrophobicity
			char c1 = seq1[i];
			char c2 = seq2[i];
			boolean isGapped = false;
			g2D.setFont(seqFont);

			List<Integer> alignedPos = null;
			if ( colorByAlignmentBlock){
				alignedPos = DisplayAFP.getEQRAlignmentPos(afpChain);
			}
			if ( isFATCAT){
				char s = symb[i];
				if ( s != ' ') {
					isGapped = false;
					g2D.setFont(eqFont);
				}
				else
					isGapped = true;
			} else {
			char s = symb[i];
				if ( c1 !=  '-' && c2 != '-' && s!= ' '){
					// no gap
					g2D.setFont(eqFont);
				} else {
					isGapped = true;

				}
			}

			Point p1 = coordManager.getPanelPos(0,i);
			int xpos1 = p1.x;
			int ypos1 = p1.y;
			Point p2 = coordManager.getPanelPos(1, i);
			int xpos2 = p2.x;
			int ypos2 = p2.y;
			int blockNum = afpChain.getBlockNum();
			if (! isGapped) {
				Color bg = Color.white;
				Color bg2 = Color.white;
				Color end1 = ColorUtils.rotateHue(ColorUtils.orange,  (1.0f  / 24.0f) * blockNum  );
				Color end2 = ColorUtils.rotateHue(ColorUtils.cyan,    (1.0f  / 24.0f) * (blockNum  +1)) ;

				if ( colorByAlignmentBlock) {

					if (! alignedPos.contains(i)){

						// unaligned!
						bg = Color.white;
						bg2 = Color.white;
					} else  {

						int colorPos = 0;
						if (isFATCAT) {
							int block = 0;
							char s = symb[i];
							try {
								block = Integer.parseInt(String.valueOf(s)) - 1;
								bg  = ColorUtils.getIntermediate(ColorUtils.orange, end1, blockNum, block);
								bg2   = ColorUtils.getIntermediate(ColorUtils.cyan, end2, blockNum, block);
								//bg = ColorUtils.rotateHue(ColorUtils.orange,  (1.0f  / 24.0f) * block  );
								//bg2 = ColorUtils.rotateHue(ColorUtils.cyan,  (1.0f  / 16.0f) * block );
							} catch (Exception e){}

							if ( colorPos > ColorUtils.colorWheel.length){
								colorPos = ColorUtils.colorWheel.length % colorPos ;
							}
						} else {
						  colorPos = AFPAlignmentDisplay.getBlockNrForAlignPos(afpChain, i);
						  bg  = ColorUtils.getIntermediate(ColorUtils.orange, end1, blockNum, colorPos);
						  bg2   = ColorUtils.getIntermediate(ColorUtils.cyan, end2, blockNum, colorPos);
						  //bg = ColorUtils.rotateHue(ColorUtils.orange,  (1.0f  / 24.0f) * colorPos );
						  //bg2 = ColorUtils.rotateHue(ColorUtils.cyan,  (1.0f  / 16.0f) * colorPos);
						}

					}
				} else {

					bg = Color.LIGHT_GRAY;
					bg2 = Color.LIGHT_GRAY;
				}

				// draw a darker background
				g2D.setPaint(bg);
				Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+1);
				g2D.fill(rec);
				g2D.setPaint(bg2);
				Rectangle rec2 = new Rectangle(p1.x-1,p1.y+4, (p2.x-p1.x)+12, (p2.y-p1.y)-3);
				g2D.fill(rec2);

				//g2D.setPaint(Color.black);
				//g2D.draw(rec);
			}
			if ( colorBySimilarity){
				if ( c1 == c2){
					Color bg = COLOR_EQUAL;
					g2D.setPaint(bg);
					Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+12);
					g2D.fill(rec);
				} else if (AFPAlignmentDisplay.aaScore(c1, c2) > 0) {
					Color bg = COLOR_SIMILAR;
					g2D.setPaint(bg);
					Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+12);
					g2D.fill(rec);
				}
			}

			//if ( selectionStart != null && selectionEnd != null){
			//	if ( i >= selectionStart.getPos1() && i <= selectionEnd.getPos1()) {

			if ( isSelected(i)){
				// draw selection
				Color bg = Color.YELLOW;
				g2D.setPaint(bg);
				// draw a darker backgroun
				Rectangle rec = new Rectangle(p1.x-1,p1.y-11, (p2.x-p1.x)+12, (p2.y-p1.y)+12);
				g2D.fill(rec);
				//	}
			}

			// draw the AA sequence
			g2D.setColor(Color.black);
			g2D.drawString(String.valueOf(c1), xpos1, ypos1);
			g2D.drawString(String.valueOf(c2), xpos2, ypos2);




			//System.out.println(seq1[i] + " " + xpos1 + " " + ypos1 + " " + seq2[i] + xpos2 + " " + ypos2);
		}

		int nrLines = (afpChain.getAlnLength() -1) / AFPChainCoordManager.DEFAULT_LINE_LENGTH;


		for ( int i = 0 ; i <= nrLines ; i++){

			try {
				// draw legend at i
				Point p1 = coordManager.getLegendPosition(i,0);
				Point p2 = coordManager.getLegendPosition(i,1);


				int aligPos = i * AFPChainCoordManager.DEFAULT_LINE_LENGTH ;
				Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0,aligPos, ca1,false);
				Atom a2 = DisplayAFP.getAtomForAligPos(afpChain, 1,aligPos, ca2,false);
				String label1 = JmolTools.getPdbInfo(a1,false);
				String label2 = JmolTools.getPdbInfo(a2,false);
				g2D.drawString(label1, p1.x,p1.y);
				g2D.drawString(label2, p2.x,p2.y);

				Point p3 = coordManager.getEndLegendPosition(i,0);
				Point p4 = coordManager.getEndLegendPosition(i,1);

				aligPos = i * AFPChainCoordManager.DEFAULT_LINE_LENGTH + AFPChainCoordManager.DEFAULT_LINE_LENGTH -1 ;
				if ( aligPos > afpChain.getAlnLength())
					aligPos = afpChain.getAlnLength() - 1;
				Atom a3 = DisplayAFP.getAtomForAligPos(afpChain, 0,aligPos, ca1,true);
				Atom a4 = DisplayAFP.getAtomForAligPos(afpChain, 1,aligPos, ca2,true);

				String label3 = JmolTools.getPdbInfo(a3,false);
				String label4 = JmolTools.getPdbInfo(a4,false);

				g2D.drawString(label3, p3.x,p3.y);
				g2D.drawString(label4, p4.x,p4.y);


			} catch (StructureException e){
				e.printStackTrace();
			}
		}


	}




	private boolean isSelected(int alignmentPosition) {

		return selection.get(alignmentPosition);

	}


	@Override
public void mouseOverPosition(AlignedPosition p) {
		//System.out.println("AligPanel: mouse over position " + p.getPos1() );

		if ( ! selectionLocked)
			selection.clear();
		selection.set(p.getPos1());

		updateJmolDisplay();

		this.repaint();

	}

	private void updateJmolDisplay() {

		if ( jmol == null)
			return;

		int size = afpChain.getAlnLength();

		StringBuffer cmd = new StringBuffer("select ");

		int nrSelected = 0;
		try {

			for (int i = 0 ; i< size ; i++){
				if ( selection.get(i)){

					Atom a1 = DisplayAFP.getAtomForAligPos(afpChain, 0,i, ca1, false);
					Atom a2 = DisplayAFP.getAtomForAligPos(afpChain, 1,i, ca2, false);

					String select1 = "";

					if ( a1 != null )
						select1 = JmolTools.getPdbInfo(a1);
					String select2 = "" ;
					if ( a2 != null)
						select2 = JmolTools.getPdbInfo(a2);

					// nothing to display
					if ( select1.equals("") && select2.equals(""))
						continue;

					if ( nrSelected > 0)
						cmd.append(", ");

					cmd.append(select1);
					cmd.append("/1, ");
					cmd.append(select2);
					cmd.append("/2");
					nrSelected++;
				}
			}


		} catch (StructureException e){
			e.printStackTrace();
		}
		if ( nrSelected == 0)
			cmd.append(" none;");
		else
			cmd.append("; set display selected;");

		jmol.evalString(cmd.toString());


	}


	@Override
public void positionSelected(AlignedPosition p) {
		mouseOverPosition(p);

	}

	@Override
public void rangeSelected(AlignedPosition start, AlignedPosition end) {
		//System.out.println("AligPanel: range selected " + start.getPos1() + " - " + end.getPos1() + " selectionLockedL " + selectionLocked);
		if ( ! selectionLocked )
			selection.clear();
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
		//System.out.println("AligPanel: toggle selection " + p.getPos1() + " " + selection.get(p.getPos1()));
		updateJmolDisplay();
		this.repaint();

	}



	public void setAlignmentJmol(AbstractAlignmentJmol jmol) {
		this.jmol = jmol;

	}


	@Override
public void windowActivated(WindowEvent e) {

		// TODO Auto-generated method stub

	}


	@Override
public void windowClosed(WindowEvent e) {
		// TODO Auto-generated method stub

	}


	@Override
public void windowClosing(WindowEvent e) {
		destroy();

	}


	@Override
public void windowDeactivated(WindowEvent e) {
		// TODO Auto-generated method stub

	}


	@Override
public void windowDeiconified(WindowEvent e) {
		// TODO Auto-generated method stub

	}


	@Override
public void windowIconified(WindowEvent e) {
		// TODO Auto-generated method stub

	}


	@Override
public void windowOpened(WindowEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
public void actionPerformed(ActionEvent e) {
		String cmd = e.getActionCommand();
		// print is handled by superclass
		if ( cmd.equals(MenuCreator.PRINT)) {
			super.actionPerformed(e);
		} else if (cmd.equals(MenuCreator.TEXT_ONLY)){
			String result = AfpChainWriter.toWebSiteDisplay(afpChain, ca1, ca2);
			DisplayAFP.showAlignmentImage(afpChain, result);
		} else if ( cmd.equals(MenuCreator.PAIRS_ONLY)) {
			String result = AfpChainWriter.toAlignedPairs(afpChain, ca1, ca2) ;
			DisplayAFP.showAlignmentImage(afpChain, result);
		} else if (cmd.equals(MenuCreator.FATCAT_TEXT)){
			String result = afpChain.toFatcat(ca1, ca2);
			result += AFPChain.newline;
			result += afpChain.toRotMat();
			DisplayAFP.showAlignmentImage(afpChain, result);
		} else if ( cmd.equals(MenuCreator.SELECT_EQR)){
			selectEQR();
		} else if ( cmd.equals(MenuCreator.SIMILARITY_COLOR)){
			colorBySimilarity(true);
		} else if ( cmd.equals(MenuCreator.EQR_COLOR)){
			colorBySimilarity(false);
		} else if ( cmd.equals(MenuCreator.FATCAT_BLOCK)){
			colorByAlignmentBlock();
		}
		else {
			System.err.println("Unknown command:" + cmd);
		}

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


	private void selectEQR() {

		selection.clear();

		List<Integer> pos1 = DisplayAFP.getEQRAlignmentPos(afpChain);

		for (int pos : pos1){
			selection.flip(pos);
		}
		mouseMoLi.triggerSelectionLocked(true);
		updateJmolDisplay();
		this.repaint();

	}

	public Atom[] getCa1() {
		return ca1;
	}


	public void setCa1(Atom[] ca1) {
		this.ca1 = ca1;
	}


	public Atom[] getCa2() {
		return ca2;
	}


	public void setCa2(Atom[] ca2) {
		this.ca2 = ca2;
	}


}
