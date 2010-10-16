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
 * Created on Jul 16, 2006
 *
 */
package org.biojava.bio.structure.gui.util;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.text.MessageFormat;
import java.util.List;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.TableCellRenderer;

import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.AtomImpl;
import org.biojava.bio.structure.Calc;
import org.biojava.bio.structure.Chain;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureImpl;
import org.biojava.bio.structure.align.StructurePairAligner;
import org.biojava.bio.structure.align.pairwise.AlternativeAlignment;
import org.biojava.bio.structure.gui.events.JmolAlignedPositionListener;
import org.biojava.bio.structure.jama.Matrix;
import org.biojava.bio.structure.gui.BiojavaJmol;
import org.biojava.bio.structure.gui.ScaleableMatrixPanel;
import org.biojava.bio.structure.gui.SequenceDisplay;


/** a frame showing the alternative alignments, which are the result of a structure superimposition
 * 
 * @author Andreas Prlic
 * @since 1.7
 * @version %I% %G%
 */
public class AlternativeAlignmentFrame 
extends JFrame{

	private static final long serialVersionUID=0l;


	public static Logger logger =  Logger.getLogger("org.biojava");

	private static String[] columnNames = new String[]{"#","eqr","score", "rms", "gaps","cluster", "show distance matrix","show alignment"};

	AlternativeAlignment[] aligs;
	JPanel panel;

	Structure structure1;
	Structure structure2;
	StructurePairAligner structurePairAligner;

	public AlternativeAlignmentFrame(Structure s1, Structure s2) {
		super();
		panel = new JPanel();
		panel.setPreferredSize(new Dimension(800,400));
		this.getContentPane().add(panel);


		structure1  = s1;
		structure2  = s2;
		String pdb1 = s1.getPDBCode();
		String pdb2 = s2.getPDBCode();

		String t = "Alternative Alignments";
		Object[] args = {pdb1,pdb2};

		String title =  MessageFormat.format(t,args);
		this.setTitle(title);
	}

	public void setStructurePairAligner(StructurePairAligner aligner){
		this.structurePairAligner = aligner;
	}
	
	public void setAlternativeAlignments(AlternativeAlignment[] aligs) {
		this.aligs = aligs;
		panel.removeAll();

		//Box vBox = Box.createVerticalBox();
		//panel.add(vBox);

		Object[][] data = getDataFromAligs(aligs);
		JTableDataButtonModel model = new JTableDataButtonModel(data, columnNames);
		JTable table = new JTable(model);
		
		
		TableCellRenderer defaultRenderer = table.getDefaultRenderer(JButton.class);
		
		JButtonTableCellRenderer myRenderer = new JButtonTableCellRenderer(defaultRenderer);
		
		table.setDefaultRenderer(JButton.class, myRenderer);
		
		table.addMouseListener(new JTableMouseButtonListener(table));
		
		JScrollPane scrollPane = new JScrollPane(table);
		scrollPane.setPreferredSize(new Dimension(800,400));
		//vBox.add(e);
		panel.add(scrollPane);
		

	}

	private Object[][] getDataFromAligs(AlternativeAlignment[] aligs){


		Object[][] data = new Object[aligs.length][columnNames.length];

		for ( int i=0;i< aligs.length;i++){
			AlternativeAlignment alig = aligs[i];

			data[i][0] = new Integer(i+1);
			data[i][1] = new Integer(alig.getEqr());
			data[i][2] = new Double(alig.getScore());
			data[i][3] = new Double(alig.getRmsd());
			data[i][4] = new Integer(alig.getGaps());
			data[i][5] = new Integer(alig.getCluster());
			JButton maxb = new JButton("Distance Matrix");
			maxb.addMouseListener(new MatrixMouseListener(this,i));

			data[i][6] = maxb;
			

			//Action action1 = new MyButtonAction(t,this,i);
			JButton but = new JButton("Show in Jmol");
			but.addMouseListener(new MyButtonMouseListener(this,i));
			data[i][7] = but;

			
		}
		return data;
	}

	public void showDistanceMatrix(int position){
		if ( position > aligs.length){
			return;
		}
		AlternativeAlignment alig = aligs[position];
		logger.info("display distance matrix for alternative alignment " + (position +1));
		
		ScaleableMatrixPanel smp = new ScaleableMatrixPanel();
		JFrame frame = new JFrame();
		frame.setTitle("Alt. Alig [" + position+"] - Distance Matrix & path");
		
		frame.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				JFrame f = (JFrame) e.getSource();
				f.setVisible(false);
				f.dispose();
			}

					
			
		});
					
		smp.setMatrix(alig.getDistanceMatrix());		
		smp.setAlternativeAligs(new AlternativeAlignment[]{alig});
		
		frame.getContentPane().add(smp);

		frame.pack();
		frame.setVisible(true);

	}

	public void showAlternative(int position){
		if ( position > aligs.length){
			return;
		}
		AlternativeAlignment alig = aligs[position];
		logger.info("display alternative alignment " + (position +1));

		// create the structure alignment object and tell the listeners ...


//		Matrix m1 = Matrix.identity(3,3);
		Matrix m2 = alig.getRotationMatrix();

		String pdb1 = structure1.getPDBCode();
		String pdb2 = structure2.getPDBCode();


		Atom shift1 = new AtomImpl();
		shift1.setCoords(new double[]{0,0,1});
		Atom shift2 = alig.getShift();

		Structure s3 = (Structure)structure2.clone();

		Calc.rotate(s3,m2);
		Calc.shift(s3,shift2);

		BiojavaJmol jmol = new BiojavaJmol();
		jmol.setTitle(pdb1 + " vs. " + pdb2);

		Structure n = new StructureImpl();
		
		List<Chain> chains1 = structure1.getChains();
		
		n.addModel(chains1);
		
		List<Chain> chains3 = s3.getChains();
		n.addModel(chains3);
	
		n.setNmr(true);
		
		jmol.setStructure(n);
		String[] cmds = createRasmolScripts(alig);
		jmol.evalString("model 0 ; select * ; wireframe off ; spacefill off; backbone 0.3;");
		jmol.evalString(cmds[0]);
		jmol.evalString(cmds[1]);		
		
		JFrame frame = new JFrame("Sequences for AlternativeAlignment ["+position+"]");
		
		SequenceDisplay seqdisp;
		seqdisp =  new SequenceDisplay(structurePairAligner);
		seqdisp.setStructure1(structure1);
		seqdisp.setStructure2(structure2);
	
		seqdisp.setAlternativeAlignment(alig);
		
		frame.getContentPane().add(seqdisp);
		
		frame.pack();
		frame.setVisible(true);
		frame.addWindowListener(new WindowAdapter(){
			public void windowClosing(WindowEvent e){
				JFrame f = (JFrame) e.getSource();
				f.setVisible(false);
				f.dispose();
			}

					
			
		});
		
		seqdisp.updateDisplay();
		
		JmolAlignedPositionListener jmolBridge = new JmolAlignedPositionListener(jmol,structurePairAligner);
		jmolBridge.setStructure1(structure1);		
		jmolBridge.setStructure2(s3);
		
		seqdisp.addAlignmentPositionListener(jmolBridge);
		
	}

	
	
	private String[] createRasmolScripts(AlternativeAlignment alig){
		String[] scripts = new String[2];

		Color col1 = Color.red;
		Color col2 = Color.blue;

		Color chaincol1 = new Color(col1.getRed()/2,col1.getGreen()/2,col1.getBlue()/2);
		Color chaincol2 = new Color(col2.getRed()/2,col2.getGreen()/2,col2.getBlue()/2);

		
		
		String cmd1 = "";
		String cmd2 = "";

		cmd1 += "select */"+1+"; ";
		cmd1 += " color [" +chaincol1.getRed()+","+chaincol1.getGreen() +","+chaincol1.getBlue() +"];";

		cmd2 += "select */"+2+"; ";
		cmd2 += " color [" +chaincol2.getRed()+","+chaincol2.getGreen() +","+chaincol2.getBlue() +"];";

		cmd1 += "select ";
		cmd2 += "select ";

		String[] pdb1s = alig.getPDBresnum1();
		String[] pdb2s = alig.getPDBresnum2();


		for ( int i =0 ; i< pdb1s.length;i++){

			String p1 = pdb1s[i];
			String p2 = pdb2s[i];

			cmd1 += p1 +"/1";
			cmd2 += p2 +"/2";

			if ( i <= pdb1s.length -2){
				cmd1 += ",";
				cmd2 += ",";
			}
		}

		cmd1 += "; color [" +col1.getRed()+","+col1.getGreen() +","+col1.getBlue() +"];";
		cmd1 += " backbone 0.6;";   

		cmd2 += "; color [" +col2.getRed()+","+col2.getGreen() +","+col2.getBlue() +"];";
		cmd2 += " backbone 0.6;";   

		//System.out.println(cmd1);
		scripts[0] = cmd1;
		scripts[1] = cmd2;

		return scripts;
	}




}

class MyButtonMouseListener implements MouseListener{
	AlternativeAlignmentFrame parent;
	int pos;
	public MyButtonMouseListener(AlternativeAlignmentFrame parent, int position){

		this.parent = parent;
		this.pos = position;
	}



	public void mouseClicked(MouseEvent arg0) {


	}

	public void mousePressed(MouseEvent arg0) {

	}

	public void mouseReleased(MouseEvent arg0) {     
		parent.showAlternative(pos);

	}

	public void mouseEntered(MouseEvent arg0) {

	}

	public void mouseExited(MouseEvent arg0) {

	}

}

class MatrixMouseListener implements MouseListener{
	AlternativeAlignmentFrame parent;
	int pos;
	public MatrixMouseListener( AlternativeAlignmentFrame parent, int position){

		this.parent = parent;
		this.pos = position;
	}

	public void mouseClicked(MouseEvent arg0) {}
	public void mousePressed(MouseEvent arg0) {}

	public void mouseReleased(MouseEvent arg0) {     
		parent.showDistanceMatrix(pos);

	}


	public void mouseEntered(MouseEvent arg0) { }

	public void mouseExited(MouseEvent arg0) {}

}

